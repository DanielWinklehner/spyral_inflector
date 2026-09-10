"""Direct BEM solve of an exported STEP assembly at given voltages (quads optionally
rotated about z), the interpolated E-field on a Cartesian grid, and a check of the
reloaded geometry against the optimizer: surface-mesh distances and the design particle
tracked in the re-solved field versus the optimizer's own track.

Writes into out_dir: ef_itp_<tag>.pickle (the field), si_state_<tag>.pickle (design orbit
and electrode voltages, what the bunch tracking needs), reload_<tag>.json and
reload_<tag>_tracks.npz.

    python -m spyral_inflector.tracking.bem_reload --step-dir STEPS --voltages voltages.csv --state state.pickle \\
        --out-dir OUT --tag final --bfield B.pickle --quad-voltages 6575 -8450 --rotate-quads 16 24 --spiral-voltage 11594.6
"""
import argparse
import json
import os
import pickle
import time

import numpy as np

from .deck import (SPECIES, Z_START, read_voltages, set_quad_voltages, set_spiral_voltage, load_state,
                   load_step_assembly, rotate_quads, rotate_assembly, mesh_assembly)

DEFAULT_BOX = (-0.10, 0.10, -0.10, 0.10, -0.29, 0.06)


def _path_length(r):
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(r, axis=0), axis=1))])


def _exit_state(state, r, v):
    p_out, n_out = state["trj_design"][-1], state["v_design"][-1] / np.linalg.norm(state["v_design"][-1])
    d = (r - p_out) @ n_out
    i = int(np.where(d >= 0.0)[0][0])
    ang = np.degrees(np.arcsin(v[i, 2] / np.linalg.norm(v[i])))
    return i, r[i].copy(), float(ang)


def _outside(r, v, s, i_exit):
    """Mean vertical angle and height 2 gaps .. 2 gaps + 25 mm past the exit (as in the optimizer)."""
    w = (s >= s[i_exit] + 0.038) & (s <= s[i_exit] + 0.063)
    vn = np.linalg.norm(v[w], axis=1)
    return float(np.mean(np.degrees(np.arcsin(v[w, 2] / vn)))), float(np.mean(r[w, 2]))


def surface_field_stats(n_fun_coeff, mesh, domain_names):
    """Field on the conductor surfaces from the BEM solution, per electrode.

    The Neumann trace of the potential on the piecewise-constant space is the normal
    derivative per triangle, i.e. the surface field |E| (the tangential field vanishes on
    a conductor). Corner and edge triangles carry mesh-dependent spikes, so besides the
    maximum the area-weighted 99th and 99.9th percentiles and the mean over the 1 % of the
    surface with the highest field are returned; those are the numbers to compare between
    geometries. Values in kV/cm."""
    verts = np.asarray(mesh["verts"], dtype=float)
    elems = np.asarray(mesh["elems"], dtype=int)
    domns = np.asarray(mesh["domns"], dtype=int).ravel()
    if verts.shape[0] == 3 and verts.shape[1] != 3:
        verts = verts.T
    if elems.shape[0] == 3 and elems.shape[1] != 3:
        elems = elems.T
    e_n = np.abs(np.asarray(n_fun_coeff, dtype=float)) * 1e-5          # V/m -> kV/cm
    tri = verts[elems]
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    out = {}
    for dom in np.unique(domns):
        m = domns == dom
        e, a = e_n[m], area[m]
        order = np.argsort(e)
        e, a = e[order], a[order]
        cum = np.cumsum(a) / a.sum()
        p99 = float(e[np.searchsorted(cum, 0.99)]) if len(e) else float("nan")
        p999 = float(e[np.searchsorted(cum, 0.999)]) if len(e) else float("nan")
        top = cum >= 0.99
        top_mean = float(np.average(e[top], weights=a[top])) if top.any() else float("nan")
        out[domain_names.get(int(dom), "domain_{}".format(dom))] = {
            "max_kv_cm": float(e.max()) if len(e) else float("nan"), "p999_kv_cm": p999, "p99_kv_cm": p99,
            "top1pct_mean_kv_cm": top_mean, "mean_kv_cm": float(np.average(e, weights=a)) if len(e) else float("nan"),
            "n_elements": int(m.sum()), "area_cm2": float(1e4 * a.sum())}
    return out


def _compare_mesh(assembly, state, log):
    """Distance of every reloaded STEP vertex to the surface the optimizer solved on."""
    from ..optimization import _surface_distance
    ov = np.asarray(state["mesh"]["verts"], float).T
    oe = np.asarray(state["mesh"]["elems"], int).T
    od = np.asarray(state["mesh"]["domns"], int)
    by_name = {nm: oe[od == int(dom)] for dom, nm in state["mesh"]["names"].items()}
    log("mesh comparison (reloaded STEP vertices -> optimizer's surface):")
    log("   {:<20s} {:>8s} {:>8s} {:>8s}".format("electrode", "mean[mm]", "p99[mm]", "max[mm]"))
    out = {}
    for e in assembly.electrodes.values():
        tris = by_name.get(e.name)
        if tris is None or len(tris) == 0:
            log("   {:<20s} no counterpart in the optimizer mesh".format(e.name))
            continue
        pts = np.asarray(e._gmsh_msh["vertices"], float)[np.unique(np.asarray(e._gmsh_msh["elements"]))]
        pts = pts[:: max(1, len(pts) // 4000)]
        d = _surface_distance(pts, ov, tris, k=16)
        out[e.name] = {"mean_mm": 1e3 * float(d.mean()), "p99_mm": 1e3 * float(np.percentile(d, 99)), "max_mm": 1e3 * float(d.max())}
        log("   {:<20s} {:8.3f} {:8.3f} {:8.3f}".format(e.name, out[e.name]["mean_mm"], out[e.name]["p99_mm"], out[e.name]["max_mm"]))
    return out


def _quad_field_check(efield, state, volts, log):
    """Transverse gradient at the centre of each quad, fitted over +-8 mm (a hyperbolic
    quad at aperture a gives G = 2 V / a^2)."""
    out = {}
    dz = state.get("shift_lab", (0.0, 0.0, 0.0))[2]
    for label, zc, name in (("quad1", -0.27 + 0.0225 + dz, "D0"), ("quad2", -0.19 + 0.0225 + dz, "D4")):
        xs = np.linspace(-8e-3, 8e-3, 9)
        ex = efield(np.column_stack([xs, np.zeros_like(xs), np.full_like(xs, zc)]))[:, 0]
        ey = efield(np.column_stack([np.zeros_like(xs), xs, np.full_like(xs, zc)]))[:, 1]
        gx, gy = np.polyfit(xs, ex, 1)[0], np.polyfit(xs, ey, 1)[0]
        e_axis = np.linalg.norm(efield(np.array([[0.0, 0.0, zc]]))[0])
        out[label] = {"z_m": zc, "dEx_dx_V_per_m2": float(gx), "dEy_dy_V_per_m2": float(gy), "E_on_axis_V_per_m": float(e_axis),
                      "G_hyperbolic_2V_over_a2": float(2.0 * abs(volts.get(name, 0.0)) / 0.013 ** 2)}
        log("{}: z = {:+.4f} m  dEx/dx = {:+.3e}  dEy/dy = {:+.3e} V/m^2, |E| on axis {:.2e} V/m".format(label, zc, gx, gy, e_axis))
    return out


def _rotz(deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rotate_state(state, deg):
    """The geometry build's state with every trajectory, velocity and mesh vertex rotated
    about the z axis by deg (counter-clockwise viewed from +z)."""
    R = _rotz(deg)
    out = dict(state)
    for key in ("trj_design", "v_design", "trj_design_full", "v_design_full", "track_r", "track_v"):
        if key in out and out[key] is not None:
            out[key] = np.asarray(out[key], dtype=float) @ R.T
    for key in ("shift", "shift_lab"):
        if key in out and out[key] is not None:
            out[key] = R @ np.asarray(out[key], dtype=float)
    if "mesh" in out and out["mesh"] is not None:
        m = dict(out["mesh"])
        v = np.asarray(m["verts"], dtype=float)
        m["verts"] = (R @ v) if v.shape[0] == 3 else (v @ R.T)
        out["mesh"] = m
    out["rotation_deg"] = float(out.get("rotation_deg", 0.0) + deg)
    return out


def solve_step_assembly(step_dir, voltages, state, out_dir, tag, bfield, res=0.0025, h=0.005,
                        quad_voltages=None, spiral_voltage=None, rotate=None, box=DEFAULT_BOX, domain_decomp=(4, 4, 4),
                        test_particle=True, nsteps_axis=2400, energy_mev=0.069337, rotate_all=0.0, log=print):
    """Solve the STEP assembly with bempp and save the field and the design state.

    voltages: dict or voltages.csv of the geometry build; quad_voltages=(q1, q2) and
    spiral_voltage override it (D0,D1 = +q1, D2,D3 = -q1, D4,D5 = +q2, D6,D7 = -q2; anode +V,
    cathode -V); rotate=(a1, a2) rotates the quads about z [deg]; rotate_all rotates the
    WHOLE assembly (every electrode, and the design orbit of the state) about z by that
    angle first, i.e. a rigid rotation of the inflector system in the magnet. state: the
    geometry build's state.pickle (design orbit). Returns the summary of reload_<tag>.json.
    """
    from PyPATools.particles import ParticleDistribution
    from PyPATools.species import IonSpecies
    from ..spyral_inflector import SpiralInflector

    os.makedirs(out_dir, exist_ok=True)
    state = load_state(state)
    if rotate_all:
        state = rotate_state(state, rotate_all)
    volts = read_voltages(voltages)
    if quad_voltages is not None:
        set_quad_voltages(volts, *quad_voltages)
    if spiral_voltage is not None:
        set_spiral_voltage(volts, spiral_voltage)

    t0 = time.time()
    assembly = load_step_assembly(step_dir, volts, name="{} from STEP".format(tag))
    log("loaded {} electrodes from {} in {:.1f} s".format(len(assembly.electrodes), step_dir, time.time() - t0))
    if rotate_all:
        rotate_assembly(assembly, rotate_all)
        log("   whole assembly rotated about z by {} deg".format(rotate_all))
    if rotate is not None and any(a != 0.0 for a in rotate):
        rotate_quads(assembly, *rotate)
        log("   quadrupoles rotated about z by {} / {} deg".format(*rotate))
    for e in assembly.electrodes.values():
        log("   {:<20s} {:+9.1f} V".format(e.name, e.voltage))
    t0 = time.time()
    mesh_assembly(assembly, h)
    mesh = assembly.get_bempp_mesh(brep_h=h)
    log("meshed: {} triangles in {:.1f} s".format(mesh["elems"].shape[1], time.time() - t0))
    mesh_cmp = _compare_mesh(assembly, state, log) if "mesh" in state else {}

    ion = ParticleDistribution(species=IonSpecies(SPECIES))
    ion.set_mean_energy_z_mev(energy_mev)
    si = SpiralInflector(ion=ion, method="numerical", solver="bempp",
                         volt=12000.0, gap=0.019, tilt=31.0, dx=0.01, sigma=0.0022, vee_shape="parabolic", ns=100,
                         aspect_ratio=2.4, rotation=0.0, debug=False, gammaAng=5.0, anglingAng=11.0)
    si.load_bfield(bfield=bfield)
    si.initialize()
    si.set_parameter(key="h", value=h)
    si.numerical_variables["objects"] = assembly
    si.numerical_variables["full mesh"] = {"verts": mesh["verts"], "elems": mesh["elems"], "domns": mesh["domns"]}
    t0 = time.time()
    si.solve()
    t_solve = time.time() - t0
    b = box
    t0 = time.time()
    si.calculate_potential(limits=((b[0], b[1]), (b[2], b[3]), (b[4], b[5])), res=res, domain_decomp=tuple(domain_decomp), overlap=0)
    si.calculate_efield()
    t_pot = time.time() - t0
    log("BEM solve {:.1f} s, potential+field {:.1f} s".format(t_solve, t_pot))

    efield = si.numerical_variables["ef_itp"]
    ef_fn = os.path.join(out_dir, "ef_itp_{}.pickle".format(tag))
    efield.save(ef_fn)
    electrode_voltages = {e.name: e.voltage for e in assembly.electrodes.values()}
    surface = surface_field_stats(si.numerical_variables["n_fun_coeff"], mesh, {int(e.bempp_domain): e.name for e in assembly.electrodes.values()})
    for name, st in sorted(surface.items(), key=lambda kv: -kv[1]["p99_kv_cm"])[:5]:
        log("surface field {:<20s} p99 {:.1f} kV/cm, top-1%-area mean {:.1f}, p99.9 {:.1f}, max {:.1f} (mesh h {:.0f} mm, {} elements)".format(
            name, st["p99_kv_cm"], st["top1pct_mean_kv_cm"], st["p999_kv_cm"], st["max_kv_cm"], 1e3 * h, st["n_elements"]))
    with open(os.path.join(out_dir, "si_state_{}.pickle".format(tag)), "wb") as fh:
        pickle.dump({"trj_design": state["trj_design"], "v_design": state["v_design"], "voltage": state.get("voltage"),
                     "electrode_voltages": electrode_voltages, "rotation_deg": float(rotate_all or 0.0)}, fh)
    log("wrote {} and si_state_{}.pickle".format(ef_fn, tag))
    quad_check = _quad_field_check(efield, state, volts, log)
    summary = {"tag": tag, "step_dir": step_dir, "voltages": electrode_voltages, "rotate_quads": list(rotate) if rotate else [0.0, 0.0],
               "rotate_all_deg": float(rotate_all or 0.0),
               "n_triangles": int(mesh["elems"].shape[1]), "res_m": res, "h_m": h, "solve_s": t_solve, "potential_s": t_pot,
               "mesh_comparison_mm": mesh_cmp, "quad_field_check": quad_check, "surface_field": surface}
    if not test_particle or "track_r" not in state:
        with open(os.path.join(out_dir, "reload_{}.json".format(tag)), "w") as fh:
            json.dump(summary, fh, indent=2)
        log("wrote reload_{}.json (no test particle)".format(tag))
        return summary

    # the design particle in the re-solved field vs the optimizer's track (same start, same dt)
    z_start, v0, dt, nsteps = state["z_start"], state["v0"], state["dt"], state["nsteps"]
    r_ref, v_ref = state["track_r"], state["track_v"]
    r_new, v_new = si.fast_track(r_start=np.array([0.0, 0.0, z_start]), v_start=np.array([0.0, 0.0, v0]), nsteps=nsteps, dt=dt)
    s_ref, s_new = _path_length(r_ref), _path_length(r_new)
    n = min(len(s_ref), len(s_new))
    dev = np.linalg.norm(r_new[:n] - r_ref[:n], axis=1)
    i_ref, p_ref, a_ref = _exit_state(state, r_ref, v_ref)
    i_new, p_new, a_new = _exit_state(state, r_new, v_new)
    ang_ref, z_ref = _outside(r_ref, v_ref, s_ref, i_ref)
    ang_new, z_new = _outside(r_new, v_new, s_new, i_new)
    log("test particle, optimizer track vs STEP-reloaded field:")
    log("   exit point [mm]      ref {}  new {}  |d| = {:.3f} mm".format(np.round(1e3 * p_ref, 2), np.round(1e3 * p_new, 2), 1e3 * np.linalg.norm(p_new - p_ref)))
    log("   exit angle at plane  ref {:+.3f}  new {:+.3f} deg".format(a_ref, a_new))
    log("   outside angle        ref {:+.3f}  new {:+.3f} deg".format(ang_ref, ang_new))
    log("   outside z            ref {:+.3f}  new {:+.3f} mm".format(1e3 * z_ref, 1e3 * z_new))
    log("   max |dr| up to exit  {:.3f} mm, over the whole track {:.3f} mm".format(1e3 * dev[:min(i_ref, i_new) + 1].max(), 1e3 * dev.max()))
    # an on-axis particle from the bunch start, through the quadrupoles
    r_ax, v_ax = si.fast_track(r_start=np.array([0.0, 0.0, Z_START]), v_start=np.array([0.0, 0.0, v0]), nsteps=nsteps_axis, dt=dt)
    i_ax, p_ax, a_ax = _exit_state(state, r_ax, v_ax)
    ang_ax, z_ax = _outside(r_ax, v_ax, _path_length(r_ax), i_ax)
    in_quads = (r_ax[:, 2] > -0.28) & (r_ax[:, 2] < -0.13)
    off_quads = np.linalg.norm(r_ax[in_quads, :2], axis=1)
    log("on-axis particle from z = {:+.3f} m: exit point {} mm, angle at plane {:+.3f} deg, outside angle {:+.3f} deg, "
        "outside z {:+.3f} mm; max transverse offset in the quads {:.4f} mm".format(
            Z_START, np.round(1e3 * p_ax, 2), a_ax, ang_ax, 1e3 * z_ax, 1e3 * off_quads.max()))
    summary["test_particle"] = {"exit_point_ref_mm": (1e3 * p_ref).tolist(), "exit_point_new_mm": (1e3 * p_new).tolist(),
                                "exit_point_diff_mm": 1e3 * float(np.linalg.norm(p_new - p_ref)),
                                "exit_angle_ref_deg": a_ref, "exit_angle_new_deg": a_new,
                                "outside_angle_ref_deg": ang_ref, "outside_angle_new_deg": ang_new,
                                "outside_z_ref_mm": 1e3 * z_ref, "outside_z_new_mm": 1e3 * z_new,
                                "max_dev_to_exit_mm": 1e3 * float(dev[:min(i_ref, i_new) + 1].max()), "max_dev_total_mm": 1e3 * float(dev.max())}
    summary["on_axis_from_bunch_start"] = {"exit_point_mm": (1e3 * p_ax).tolist(), "exit_angle_deg": a_ax, "outside_angle_deg": ang_ax,
                                           "outside_z_mm": 1e3 * z_ax, "max_offset_in_quads_mm": 1e3 * float(off_quads.max())}
    with open(os.path.join(out_dir, "reload_{}.json".format(tag)), "w") as fh:
        json.dump(summary, fh, indent=2)
    np.savez_compressed(os.path.join(out_dir, "reload_{}_tracks.npz".format(tag)),
                        r_ref=r_ref, v_ref=v_ref, r_new=r_new, v_new=v_new, r_axis=r_ax, v_axis=v_ax,
                        trj_design=state["trj_design"], trj_design_full=state.get("trj_design_full", state["trj_design"]),
                        mesh_verts=np.asarray(mesh["verts"]), mesh_elems=np.asarray(mesh["elems"]), mesh_domns=np.asarray(mesh["domns"]),
                        mesh_names=np.array([e.name for e in assembly.electrodes.values()]),
                        mesh_domain_ids=np.array([e.bempp_domain for e in assembly.electrodes.values()]))
    log("wrote reload_{}.json and reload_{}_tracks.npz".format(tag, tag))
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description="BEM solve of a STEP assembly, field pickle and design-state for the bunch tracking")
    p.add_argument("--step-dir", required=True)
    p.add_argument("--voltages", required=True, help="voltages.csv of the geometry build")
    p.add_argument("--state", required=True, help="state.pickle of the geometry build (design orbit)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--bfield", required=True)
    p.add_argument("--tag", default="final")
    p.add_argument("--res", type=float, default=0.0025, help="potential grid resolution [m]")
    p.add_argument("--h", type=float, default=0.005, help="surface mesh size [m]")
    p.add_argument("--quad-voltages", type=float, nargs=2, default=None, metavar=("Q1", "Q2"))
    p.add_argument("--spiral-voltage", type=float, default=None, help="anode +V, cathode -V; 0 for a quad-only basis field")
    p.add_argument("--rotate-quads", type=float, nargs=2, default=None, metavar=("A1", "A2"), help="45 deg = skew basis")
    p.add_argument("--rotate-all", type=float, default=0.0, help="rigid rotation of the whole assembly and the design orbit about z [deg]")
    p.add_argument("--no-test", action="store_true", help="skip the test particle (basis fields without the spiral field)")
    p.add_argument("--box", type=float, nargs=6, default=list(DEFAULT_BOX), metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"))
    p.add_argument("--domain-decomp", type=int, nargs=3, default=[4, 4, 4])
    p.add_argument("--nsteps-axis", type=int, default=2400)
    p.add_argument("--energy-mev", type=float, default=0.069337)
    a = p.parse_args(argv)
    return solve_step_assembly(a.step_dir, a.voltages, a.state, a.out_dir, a.tag, a.bfield, res=a.res, h=a.h,
                               quad_voltages=a.quad_voltages, spiral_voltage=a.spiral_voltage, rotate=a.rotate_quads,
                               box=tuple(a.box), domain_decomp=tuple(a.domain_decomp), test_particle=not a.no_test,
                               nsteps_axis=a.nsteps_axis, energy_mev=a.energy_mev, rotate_all=a.rotate_all,
                               log=lambda m: print(m, flush=True))


if __name__ == "__main__":
    main()
