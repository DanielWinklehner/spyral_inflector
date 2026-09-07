"""Track an RFQ bunch through the STEP geometry, without or with space charge.

Vacuum field: the BEM solution of the STEP assembly (ef_itp_<tag>.pickle from
bem_reload) or a superposition of basis solves. Space charge: PyPATools' PyAMG Poisson
solver on a Cartesian grid spanning the assembly with every STEP conductor grounded; the
solve is repeated every few steps with the live particles as the charge distribution, one
RFQ bunch (I / f_RF) tracked in isolation. Both fields superpose exactly (CompositeField).

Transmission counts particles that cross the design exit plane outward near the design
exit point AND clear the housing exit opening on the way out. Transmitted particles coast
on for `coast` steps (trajectory plots, the asymptotic state 120 steps = ~31 mm past the
crossing, and the hand-off plane crossing for the openPMD files).

    python -m spyral_inflector.tracking.bunch --tag final --reload-dir OUT --step-dir STEPS --particles core.txt \\
        --bfield B.pickle --phi 90 --n 43969 [--sc --h 0.0015 --current-ma 8] [--save-openpmd handoff.h5 --save-mode both]
"""
import argparse
import json
import logging
import os
import time

import numpy as np

from PyPATools.field import Field, CompositeField
from PyPATools.pusher import Pusher
from PyPATools.trackers import Tracker

from .deck import (load_particles, orient_beam, load_step_assembly, mesh_assembly, drop_electrodes, load_bfield,
                   load_state, superpose_basis)
from .hooks import (ElectrodeCollision, ExitPlane, PlaneCrossing, TrajectoryRecorder, SnapshotRecorder, Envelope,
                    SpaceCharge, continue_design, _PD)
from .handoff import save_handoff_openpmd, save_snapshot_openpmd


def track_bunch(reload_dir, tag, step_dir, particles, bfield, out_dir=None, out_tag=None, n=10000, phi=0.0, swap_xy=False,
                sc=False, current_ma=8.0, rf_mhz=32.8, h=2.0e-3, pad=0.01, xy_max=None, resolve_every=8, tol=1e-5, gpu=True,
                nsteps=1900, dt=1.0e-10, coast=300, asym_steps=120, post_exit_steps=100, record=1000, exclude=(),
                superpose=None, vscale=1.0, basis_dir=None, unit=3500.0,
                save_openpmd=None, save_mode="plane", handoff_distance=0.030, phase_reference="mean",
                seed=20260905, reference=None, plot=True, log=None):
    """Track n particles of the RFQ file through the geometry; returns the summary dict
    and writes bunch_<out_tag>.json / .npz / .png into out_dir (default: reload_dir).

    tag names ef_itp_<tag>.pickle and si_state_<tag>.pickle in reload_dir. With
    superpose=(q1, a1, q2, a2) the vacuum field is instead summed from the basis solves
    of basis_dir (default reload_dir) and tag only names the design state. sc=True adds
    PyAMG space charge (cell size h, current_ma at rf_mhz). save_openpmd writes the
    hand-off file(s) at handoff_distance of design path past the exit (save_mode plane,
    lab6d or both). reference: a vacuum run's bunch json for the comparison plot.
    """
    log = log or (lambda m: print(m, flush=True))
    out_dir = out_dir or reload_dir
    out_tag = out_tag or ("{}_sc{:g}mA".format(tag, current_ma) if sc else tag)
    os.makedirs(out_dir, exist_ok=True)
    efield_fn = os.path.join(reload_dir, "ef_itp_{}.pickle".format(tag))
    state_fn = os.path.join(reload_dir, "si_state_{}.pickle".format(tag))

    log("=" * 78)
    log("BUNCH TRACKING {}, tag '{}' -> '{}'".format(
        "WITH SPACE CHARGE ({:g} mA at {:g} MHz)".format(current_ma, rf_mhz) if sc else "(no space charge)", tag, out_tag))
    log("=" * 78)
    r0, v0, ion, raw = load_particles(particles, n, seed=seed)
    r0, v0 = orient_beam(r0, v0, phi, swap_xy)
    log("  particles      : {:,d} (file has {:,d} unlost) from {}{}{}".format(
        len(r0), len(raw), os.path.basename(particles), ", x and y swapped" if swap_xy else "",
        ", rotated {:+.1f} deg".format(phi) if phi else ""))
    log("  B-field        : {}".format(os.path.basename(bfield)))
    log("  mean energy    : {:.6f} MeV, start centroid {} m, z range {:+.4f}..{:+.4f} m".format(
        raw[:, 8].mean(), np.round(r0.mean(axis=0), 5), r0[:, 2].min(), r0[:, 2].max()))

    assembly = load_step_assembly(step_dir)
    if exclude:
        drop_electrodes(assembly, exclude)
        log("  EXCLUDED from collisions and the Poisson boundary: {}".format(sorted(exclude)))
    n_tri = mesh_assembly(assembly)
    log("  electrodes     : {} ({:,d} triangles)".format(len(assembly.electrodes), n_tri))

    b_field = load_bfield(bfield)
    if superpose:
        q1v, a1, q2v, a2 = superpose
        bdir = basis_dir or reload_dir
        e_vac = superpose_basis(bdir, q1v, a1, q2v, a2, vscale=vscale, unit=unit)
        efield_fn = "superposed from {}: spiral x{:.4f}, q1 {:+.0f} V at {:.0f} deg, q2 {:+.0f} V at {:.0f} deg".format(
            bdir, vscale, q1v, a1, q2v, a2)
    else:
        e_vac = Field.from_file(efield_fn)
    state = load_state(state_fn)
    trj, vdes = state["trj_design"], state["v_design"]
    r_exit = float(np.linalg.norm(trj[-1][:2]))
    log("  E-field        : {}".format(efield_fn))
    log("  voltages       : {}".format({k: round(v) for k, v in state["electrode_voltages"].items()}))
    log("  design exit    : r = {:.4f} m, point {} m".format(r_exit, np.round(trj[-1], 4)))

    # ------------------------------------------------------------------ space charge
    e_total = CompositeField([e_vac, Field.zero(dim=3)], [1.0, 1.0])
    interactions = []
    sc_obj = None
    sc_info = {}
    if sc:
        from PyPATools.poisson_amg import PyAMGSolverConfig, PyAMGPoissonSolver

        bbox = np.array(assembly.bounding_box(), dtype=float)   # xmin ymin zmin xmax ymax zmax
        lo = bbox[:3] - pad
        hi = bbox[3:] + pad
        lo[2] = min(lo[2], r0[:, 2].min() - 0.005)                # the bunch must start inside the box
        if xy_max is not None:
            lo[:2] = np.maximum(lo[:2], -xy_max)
            hi[:2] = np.minimum(hi[:2], xy_max)
        extent = hi - lo
        origin = 0.5 * (hi + lo)
        cells = tuple(int(np.ceil(e / h)) for e in extent)
        log("  assembly bbox  : {} .. {} m".format(np.round(bbox[:3], 4), np.round(bbox[3:], 4)))
        log("  Poisson box    : {} .. {} m, cells {} (h = {} mm), {:,d} DOFs".format(
            np.round(lo, 4), np.round(hi, 4), cells, np.round(1e3 * extent / np.array(cells), 2), int(np.prod(cells))))
        config = PyAMGSolverConfig(domain_extent=tuple(extent), mesh_cells=cells, domain_origin=tuple(origin), use_gpu=gpu, solver_tol=tol)
        t0 = time.time()
        solver = PyAMGPoissonSolver(config, assembly)
        t_setup = time.time() - t0
        n_cond = int((solver.cell_type == 2).sum())
        n_bnd = int((solver.cell_type == 1).sum())
        log("  Poisson setup  : {:.1f} s; conductor cells {:,d}, boundary cells {:,d}, active DOFs {:,d}".format(
            t_setup, n_cond, n_bnd, solver.n_active_dofs))
        q_bunch = current_ma * 1e-3 / (rf_mhz * 1e6)
        charges = np.full(len(r0), q_bunch / len(r0))
        log("  bunch charge   : {:.4e} C = {:.4e} C per macro-particle".format(q_bunch, charges[0]))
        sc_obj = SpaceCharge(solver, e_total, charges, resolve_every)
        sc_obj.solve_now(r0, np.ones(len(r0), dtype=bool), 0, verbose=True)   # field for the first steps
        interactions.append(sc_obj)
        sc_info = {"current_mA": current_ma, "rf_MHz": rf_mhz, "bunch_charge_C": q_bunch, "h_m": h, "cells": list(cells),
                   "box_lo_m": lo.tolist(), "box_hi_m": hi.tolist(), "n_conductor_cells": n_cond, "n_boundary_cells": n_bnd,
                   "active_dofs": int(solver.n_active_dofs), "setup_time_s": t_setup, "resolve_every": resolve_every,
                   "tol": tol, "gpu": bool(solver.use_gpu)}

    # ------------------------------------------------------------------ track
    collision = ElectrodeCollision(assembly)
    exit_plane = ExitPlane(r_exit, point=trj[-1], normal=vdes[-1], coast_steps=coast, asym_steps=asym_steps)
    collision.skip = exit_plane
    collision.post_exit_electrodes = [(i, e) for i, e in enumerate(assembly.electrodes.values()) if "housing" in e.name.lower()]
    collision.post_exit_steps = post_exit_steps
    ho_point, ho_normal = continue_design(ion, e_vac, b_field, trj[-1], vdes[-1], dt, handoff_distance)
    log("  hand-off plane : {:.0f} mm of design path past the exit, origin {} m, normal {}".format(
        1e3 * handoff_distance, np.round(ho_point, 4), np.round(ho_normal, 4)))
    handoff = PlaneCrossing(point=ho_point, normal=ho_normal, gate=exit_plane, dt=dt)
    snapshot = SnapshotRecorder(range(1000, nsteps, 20)) if save_mode in ("lab6d", "both") else None
    if sc_obj:
        sc_obj.exclude = exit_plane
    envelope = Envelope()
    trajectories = TrajectoryRecorder(len(r0), n_record=record)
    tracker = Tracker(Pusher(ion, algorithm="rk4_rel"), e_total, b_field, interactions=interactions,
                      terminators=[exit_plane, collision, handoff],
                      recorders=[envelope, trajectories] + ([snapshot] if snapshot is not None else []))
    log("  tracking {:,d} steps of {:.1e} s{} ...".format(nsteps, dt, ", Poisson solve every {} steps".format(resolve_every) if sc else ""))
    t0 = time.time()
    result = tracker.run(_PD(r0, v0), dt, nsteps, show_progress=False, sync_back=False)
    wall = time.time() - t0
    log("  done in {:.1f} s{}".format(wall, " ({} Poisson solves, {:.1f} s)".format(sc_obj.n_solves, sc_obj.solve_time) if sc_obj else ""))

    # ------------------------------------------------------------------ diagnostics
    index_map = {i: uuid for i, uuid in enumerate(assembly.electrodes.keys())}
    losses = collision.losses_by_electrode(index_map)
    post_exit = collision.post_exit_hit if collision.post_exit_hit is not None else np.zeros(len(r0), dtype=bool)
    crossed_plane = exit_plane.crossed
    crossed = crossed_plane & ~post_exit            # through the housing exit opening
    if post_exit.any():
        losses["Housing_exit"] = int(post_exit.sum())
    log("  housing exit opening: {:,d} of the {:,d} particles that crossed the exit plane hit it ({:.2f} % of the bunch)".format(
        int(post_exit.sum()), int(crossed_plane.sum()), 100.0 * post_exit.sum() / len(r0)))
    n_lost = int((collision.hit_electrode >= 0).sum())
    n_neither = int(len(r0) - crossed.sum() - n_lost)
    z_exit = exit_plane.state[crossed, 2]
    env = np.array(envelope.rows)
    hit_z = collision.hit_point[:, 2]
    loss_z = {name: hit_z[(collision.hit_electrode == i)] for i, name in ((i, assembly.electrodes[u].name) for i, u in index_map.items())}

    summary = {
        "tag": tag, "out_tag": out_tag, "space_charge": bool(sc), "n_particles": int(len(r0)),
        "n_transmitted": int(crossed.sum()), "transmission": float(crossed.sum()) / len(r0),
        "n_lost": n_lost, "loss_fraction": n_lost / len(r0), "n_still_inside_at_end": n_neither,
        "losses_by_electrode": losses,
        "loss_z_mean_by_electrode": {k: float(np.mean(v)) for k, v in loss_z.items() if len(v)},
        "z_exit_mean_m": float(np.mean(z_exit)) if z_exit.size else None,
        "z_exit_rms_m": float(np.std(z_exit)) if z_exit.size else None,
        "z_exit_p5_p95_m": [float(np.percentile(z_exit, 5)), float(np.percentile(z_exit, 95))] if z_exit.size else None,
        "mean_exit_step": float(np.mean(exit_plane.step[crossed])) if crossed.any() else None,
        "exit_offset_mean_m": float(np.nanmean(exit_plane.offset[crossed])) if crossed.any() else None,
        "voltages": state["electrode_voltages"], "efield": efield_fn, "particles": particles, "step_dir": step_dir,
        "phi_deg": phi, "swap_xy": bool(swap_xy), "nsteps": nsteps, "dt": dt, "wall_time_s": wall,
    }
    if superpose:
        summary["superpose"] = {"q1": q1v, "alpha1": a1, "q2": q2v, "alpha2": a2, "vscale": vscale, "unit": unit, "basis_dir": basis_dir or reload_dir}
    if sc_obj:
        sc_info.update({"n_solves": sc_obj.n_solves, "solve_time_s": sc_obj.solve_time,
                        "phi_min_V": float(min(r[2] for r in sc_obj.log)), "phi_max_V": float(max(r[3] for r in sc_obj.log)),
                        "e_sc_max_on_beam_V_per_m": float(max(r[4] for r in sc_obj.log))})
        summary["sc"] = sc_info

    ref_fn = reference or os.path.join(reload_dir, "bunch_{}.json".format(tag))
    ref = json.load(open(ref_fn)) if (ref_fn != os.path.join(out_dir, "bunch_{}.json".format(out_tag)) and os.path.exists(ref_fn)) else None
    log("")
    log("  TRANSMITTED     : {:,d} / {:,d}  ({:.2f} %){}".format(
        summary["n_transmitted"], len(r0), 100 * summary["transmission"],
        "   [reference run: {:.2f} %]".format(100 * ref["transmission"]) if ref else ""))
    log("  intercepted     : {:,d}  ({:.2f} %)".format(n_lost, 100 * summary["loss_fraction"]))
    log("  neither         : {:,d}  (still inside when tracking ended)".format(n_neither))
    if z_exit.size:
        log("  exit z          : mean {:+.3f} mm, rms {:.3f} mm, 5-95 % {:+.2f}..{:+.2f} mm".format(
            1e3 * summary["z_exit_mean_m"], 1e3 * summary["z_exit_rms_m"], 1e3 * summary["z_exit_p5_p95_m"][0], 1e3 * summary["z_exit_p5_p95_m"][1]))
    log("  losses by electrode (count, %, mean z of the hit{}):".format(", reference %" if ref else ""))
    for name, count in sorted(losses.items(), key=lambda kv: -kv[1]):
        line = "      {:<22s} {:>7,d}  ({:5.2f} %)   z = {:+.3f} m".format(
            name, count, 100.0 * count / len(r0), summary["loss_z_mean_by_electrode"].get(name, float("nan")))
        if ref:
            line += "   [{:5.2f} %]".format(100.0 * ref["losses_by_electrode"].get(name, 0) / ref["n_particles"])
        log(line)
    if sc_obj:
        log("  space charge    : {} solves, {:.1f} s total ({:.2f} s each); phi {:+.1f}..{:+.1f} V; |E_sc| max on the beam {:.2e} V/m".format(
            sc_obj.n_solves, sc_obj.solve_time, sc_obj.solve_time / max(1, sc_obj.n_solves),
            sc_info["phi_min_V"], sc_info["phi_max_V"], sc_info["e_sc_max_on_beam_V_per_m"]))
    log("=" * 78)

    if save_openpmd:
        extra = {"source_tag": tag, "step_dir": step_dir, "transmission_through_housing": float(crossed.mean()),
                 "n_tracked": int(len(r0)), "space_charge": bool(sc), "beam_rotation_deg": float(phi), "swap_xy": bool(swap_xy)}
        if save_mode in ("plane", "both"):
            info = save_handoff_openpmd(save_openpmd, handoff, crossed, ion, rf_mhz * 1e6, current_ma, trj, vdes,
                                        handoff_distance, phase_reference, extra_meta=extra)
            summary["handoff_openpmd"] = dict(info, path=save_openpmd, handoff_distance_m=handoff_distance)
            log("  hand-off plane {:.0f} mm past the exit: {:,d} particles written to {} (phase rms {:.1f} deg, u/v rms {:.2f}/{:.2f} mm)".format(
                1e3 * handoff_distance, info["n"], save_openpmd, info["phase_rms_deg"], info["u_rms_mm"], info["v_rms_mm"]))
        if save_mode in ("lab6d", "both") and snapshot is not None and snapshot.snapshots:
            sel6 = crossed & handoff.crossed
            k_mean = float(np.mean(handoff.step[sel6])) if sel6.any() else nsteps - 1
            k = min(snapshot.snapshots, key=lambda kk: abs(kk - k_mean))
            p6 = save_openpmd if save_mode == "lab6d" else os.path.splitext(save_openpmd)[0] + "_lab6d.h5"
            info6 = save_snapshot_openpmd(p6, snapshot.snapshots[k], crossed, ion, rf_mhz * 1e6, current_ma, extra_meta=dict(extra, snapshot_step=int(k)))
            summary["lab6d_openpmd"] = dict(info6, path=p6, step=int(k))
            log("  6-D snapshot at step {} ({:.1f} ns): {:,d} particles written to {}".format(k, 1e9 * info6["t_s"], info6["n"], p6))

    with open(os.path.join(out_dir, "bunch_{}.json".format(out_tag)), "w") as fh:
        json.dump(summary, fh, indent=2)
    np.savez_compressed(os.path.join(out_dir, "bunch_{}.npz".format(out_tag)),
                        envelope=env, z_exit=z_exit, exit_state=exit_plane.state[crossed], exit_step=exit_plane.step, crossed=crossed,
                        hit_electrode=collision.hit_electrode, hit_step=collision.hit_step, hit_point=collision.hit_point, active=result.active,
                        electrode_names=np.array([assembly.electrodes[u].name for u in index_map.values()]),
                        r0=r0, v0=v0, sc_log=np.array(sc_obj.log) if sc_obj else np.zeros((0, 7)),
                        traj=trajectories.data, traj_idx=trajectories.idx, traj_every=trajectories.every,
                        asym_state=exit_plane.asym_state, asym_steps=exit_plane.asym_steps, dt=dt,
                        handoff_state=(handoff.state if handoff.state is not None else np.full((len(r0), 6), np.nan)),
                        handoff_time=(handoff.time if handoff.time is not None else np.full(len(r0), np.nan)),
                        post_exit_hit=post_exit,
                        post_exit_point=(collision.post_exit_point if collision.post_exit_point is not None else np.full((len(r0), 3), np.nan)))
    if plot:
        _plot_bunch(os.path.join(out_dir, "bunch_{}.png".format(out_tag)), summary, env, z_exit, losses, hit_z, collision, sc_obj, state, trj,
                    ref, os.path.splitext(ref_fn)[0] + ".npz" if ref else None, len(r0), out_tag, sc, current_ma, rf_mhz, h)
    log("wrote bunch_{}.json/.npz{} in {}".format(out_tag, "/.png" if plot else "", out_dir))
    return summary


def _plot_bunch(out_png, summary, env, z_exit, losses, hit_z, collision, sc_obj, state, trj, ref, ref_npz, n0, out_tag, sc, current_ma, rf_mhz, h):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref_env = np.load(ref_npz)["envelope"] if ref_npz and os.path.exists(ref_npz) else None
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    ax = axs[0, 0]
    zc = env[:, 5] * 1e3
    ax.plot(zc, 1e3 * env[:, 6], color="C0", label="rms x")
    ax.plot(zc, 1e3 * env[:, 7], color="C1", label="rms y")
    if ref_env is not None:
        ax.plot(ref_env[:, 5] * 1e3, 1e3 * ref_env[:, 6], color="C0", ls=":", label="rms x, reference")
        ax.plot(ref_env[:, 5] * 1e3, 1e3 * ref_env[:, 7], color="C1", ls=":", label="rms y, reference")
    ax.axvline(1e3 * trj[0][2], color="k", lw=0.8, ls=":", label="inflector entrance")
    ax.set_xlabel("mean z of live particles [mm]")
    ax.set_ylabel("[mm]")
    ax.set_title("Beam envelope along the path (live particles)")
    ax.legend(fontsize=8)
    ax2 = ax.twinx()
    ax2.plot(zc, env[:, 2], color="0.4", lw=0.8)
    if ref_env is not None:
        ax2.plot(ref_env[:, 5] * 1e3, ref_env[:, 2], color="0.4", lw=0.8, ls=":")
    ax2.set_ylabel("live particles", color="0.4")

    ax = axs[0, 1]
    names = sorted(set(losses) | (set(ref["losses_by_electrode"]) if ref else set()), key=lambda nm: -losses.get(nm, 0))
    xs = np.arange(len(names))
    ax.bar(xs - 0.2, [100.0 * losses.get(nm, 0) / n0 for nm in names], width=0.4, color="#d62728", label="space charge" if sc else "this run")
    if ref:
        ax.bar(xs + 0.2, [100.0 * ref["losses_by_electrode"].get(nm, 0) / ref["n_particles"] for nm in names], width=0.4, color="#7f7f7f", label="reference")
    ax.set_xticks(xs)
    ax.set_xticklabels(names, rotation=60, fontsize=8)
    ax.set_ylabel("lost [% of bunch]")
    ax.set_title("Transmission {:.1f} %{}: losses by electrode".format(
        100 * summary["transmission"], " (reference {:.1f} %)".format(100 * ref["transmission"]) if ref else ""))
    ax.legend(fontsize=8)

    ax = axs[0, 2]
    if z_exit.size:
        ax.hist(1e3 * z_exit, bins=60, color="#1f77b4", alpha=0.8, label="this run")
    if ref_npz and os.path.exists(ref_npz):
        ax.hist(1e3 * np.load(ref_npz)["z_exit"], bins=60, histtype="step", color="k", label="reference")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("z at the exit plane [mm]")
    ax.set_ylabel("particles")
    ax.set_title("Vertical position at the design exit plane")
    ax.legend(fontsize=8)

    ax = axs[1, 0]
    hits = collision.hit_electrode >= 0
    if hits.any():
        ax.hist(1e3 * hit_z[hits], bins=80, color="#7f7f7f")
    ax.set_xlabel("z of the interception [mm]")
    ax.set_ylabel("particles lost")
    ax.set_title("Where particles are lost")

    ax = axs[1, 1]
    if sc_obj and sc_obj.log:
        lg = np.array(sc_obj.log)
        ax.plot(lg[:, 0], lg[:, 4], color="C3", label="|E_sc| max on the beam")
        ax.plot(lg[:, 0], lg[:, 5], color="C3", ls="--", label="|E_sc| median on the beam")
        ax.set_yscale("log")
        ax.set_xlabel("step")
        ax.set_ylabel("V/m")
        ax.legend(fontsize=8, loc="upper right")
        ax3 = ax.twinx()
        ax3.plot(lg[:, 0], lg[:, 3] - lg[:, 2], color="0.4", lw=0.8)
        ax3.set_ylabel("space-charge potential span [V]", color="0.4")
        ax.set_title("Space-charge field seen by the bunch")
        ax = axs[1, 2]
        ax.plot(lg[:, 0], lg[:, 6], ".", ms=3)
        ax.set_xlabel("step")
        ax.set_ylabel("Poisson solve wall time [s]")
        ax.set_title("{} solves, {:.1f} s total".format(sc_obj.n_solves, sc_obj.solve_time))
    else:
        ax.axis("off")
        axs[1, 2].axis("off")
    fig.suptitle("bunch tracking '{}': {:,d} particles, {}quad voltages {}".format(
        out_tag, n0, "{:g} mA at {:g} MHz, h = {:g} mm, ".format(current_ma, rf_mhz, 1e3 * h) if sc else "",
        {k: int(v) for k, v in state["electrode_voltages"].items() if k.startswith("D")}))
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description="track an RFQ bunch through the STEP geometry, with or without space charge")
    p.add_argument("--reload-dir", required=True, help="folder with ef_itp_<tag>.pickle and si_state_<tag>.pickle")
    p.add_argument("--tag", default="final")
    p.add_argument("--step-dir", required=True)
    p.add_argument("--particles", required=True, help="RFQ output file (TenThousandRFQParticles.txt columns)")
    p.add_argument("--bfield", required=True)
    p.add_argument("--out-dir", default=None, help="default: --reload-dir")
    p.add_argument("--out-tag", default=None)
    p.add_argument("--n", type=int, default=10000)
    p.add_argument("--phi", type=float, default=0.0, help="rotate the input beam about z [deg]")
    p.add_argument("--swap-xy", action="store_true")
    p.add_argument("--sc", action="store_true", help="PyAMG space charge")
    p.add_argument("--current-ma", type=float, default=8.0)
    p.add_argument("--rf-mhz", type=float, default=32.8)
    p.add_argument("--h", type=float, default=2.0e-3, help="Poisson cell size [m]")
    p.add_argument("--pad", type=float, default=0.01)
    p.add_argument("--xy-max", type=float, default=None)
    p.add_argument("--resolve-every", type=int, default=8)
    p.add_argument("--tol", type=float, default=1e-5)
    p.add_argument("--cpu", action="store_true", help="no GPU in the Poisson solver")
    p.add_argument("--nsteps", type=int, default=1900)
    p.add_argument("--dt", type=float, default=1.0e-10)
    p.add_argument("--coast", type=int, default=300)
    p.add_argument("--asym-steps", type=int, default=120)
    p.add_argument("--post-exit-steps", type=int, default=100)
    p.add_argument("--record", type=int, default=1000)
    p.add_argument("--exclude", default="", help="comma-separated electrode names to drop")
    p.add_argument("--superpose", type=float, nargs=4, default=None, metavar=("Q1", "A1", "Q2", "A2"))
    p.add_argument("--vscale", type=float, default=1.0)
    p.add_argument("--basis-dir", default=None)
    p.add_argument("--unit", type=float, default=3500.0)
    p.add_argument("--save-openpmd", default=None)
    p.add_argument("--save-mode", choices=["plane", "lab6d", "both"], default="plane")
    p.add_argument("--handoff-distance", type=float, default=0.030)
    p.add_argument("--phase-reference", choices=["mean", "median"], default="mean")
    p.add_argument("--seed", type=int, default=20260905)
    p.add_argument("--reference", default=None, help="bunch json of a vacuum run for the comparison plot")
    p.add_argument("--no-plot", action="store_true")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.WARNING)
    return track_bunch(a.reload_dir, a.tag, a.step_dir, a.particles, a.bfield, out_dir=a.out_dir, out_tag=a.out_tag, n=a.n,
                       phi=a.phi, swap_xy=a.swap_xy, sc=a.sc, current_ma=a.current_ma, rf_mhz=a.rf_mhz, h=a.h, pad=a.pad,
                       xy_max=a.xy_max, resolve_every=a.resolve_every, tol=a.tol, gpu=not a.cpu, nsteps=a.nsteps, dt=a.dt,
                       coast=a.coast, asym_steps=a.asym_steps, post_exit_steps=a.post_exit_steps, record=a.record,
                       exclude=[s.strip() for s in a.exclude.split(",") if s.strip()], superpose=a.superpose, vscale=a.vscale,
                       basis_dir=a.basis_dir, unit=a.unit, save_openpmd=a.save_openpmd, save_mode=a.save_mode,
                       handoff_distance=a.handoff_distance, phase_reference=a.phase_reference, seed=a.seed,
                       reference=a.reference, plot=not a.no_plot)


if __name__ == "__main__":
    main()
