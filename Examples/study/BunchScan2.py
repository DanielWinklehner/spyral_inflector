"""Generalized bunch scan by field superposition: quadrupole voltages, quadrupole
rotation angles, spiral voltage scale and a rotation of the input beam about z.

Basis fields (written by TrackFromStep.py, all on the same grid):
    spiral  : spiral electrodes at the optimized voltage, quads at 0 V   (tag quads_off)
    q1n/q2n : one quad at +-UNIT volts, everything else 0 V              (q1unit, q2unit)
    q1s/q2s : the same quad rotated by 45 deg about z                    (q1skew, q2skew)
A quad rotated by alpha is cos(2 alpha) * normal + sin(2 alpha) * skew, exact for the
quadrupole term of the field (the 12-pole and fringe terms rotate differently, so a
chosen design point should be confirmed with a direct solve, TrackFromStep.py
--rotate-quads).

The beam rotation phi turns x, y, x', y' about the z axis (a rotation of the RFQ
relative to the inflector). Positive angles are counter-clockwise seen along +z.

    python BunchScan2.py --particles ... --swap-xy --phi 0 30 60 90 120 150 --q1 -9000 9000 7 --q2 -9000 9000 7
    python BunchScan2.py ... --alpha1 0 10 20 --alpha2 0 10 20 30 --points 3000 -3000
"""
import argparse
import itertools
import json
import os
import pickle
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import track_inflector as ti  # noqa: E402
from spyral_inflector.tracking.deck import rotate_assembly, rotate_quads, shift_assembly  # noqa: E402
from PyPATools.field import Field  # noqa: E402
from PyPATools.pusher import Pusher  # noqa: E402
from PyPATools.trackers import Tracker  # noqa: E402

DECK = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

parser = argparse.ArgumentParser()
parser.add_argument("--reload-dir", default=os.path.join(DECK, "Results", "reload"))
parser.add_argument("--out-dir", default=None, help="default: the reload dir")
parser.add_argument("--step-dir", default=os.path.join(DECK, "Geometry", "final_steps"))
parser.add_argument("--bfield", default=None, help="B-field pickle (default: track_inflector.BFIELD)")
parser.add_argument("--basis", nargs="+", default=["quads_off", "q1unit", "q1skew", "q2unit", "q2skew"], metavar="NAME")
parser.add_argument("--unit", type=float, default=3500.0, help="quad voltage of the basis fields [V]")
parser.add_argument("--q1", type=float, nargs=3, default=[-9000, 9000, 7], metavar=("MIN", "MAX", "N"))
parser.add_argument("--q2", type=float, nargs=3, default=[-9000, 9000, 7], metavar=("MIN", "MAX", "N"))
parser.add_argument("--q3", type=float, nargs=3, default=None, metavar=("MIN", "MAX", "N"), help="third quad (needs 7 basis names)")
parser.add_argument("--alpha3", type=float, nargs="+", default=[0.0], help="quad 3 rotation angles [deg]")
parser.add_argument("--points", type=float, nargs="+", default=None, help="explicit q1 q2 pairs instead of the grid")
parser.add_argument("--alpha1", type=float, nargs="+", default=[0.0], help="quad 1 rotation angles [deg]")
parser.add_argument("--alpha2", type=float, nargs="+", default=[0.0], help="quad 2 rotation angles [deg]")
parser.add_argument("--phi", type=float, nargs="+", default=[0.0], help="beam rotation angles [deg]")
parser.add_argument("--vscale", type=float, nargs="+", default=[1.0], help="spiral voltage scale factors")
parser.add_argument("--particles", default=None)
parser.add_argument("--swap-xy", action="store_true")
parser.add_argument("--n", type=int, default=2000)
parser.add_argument("--nsteps", type=int, default=1750)
parser.add_argument("--post-exit-steps", type=int, default=100,
                    help="steps after the exit-plane crossing during which particles are still checked against the housing "
                         "(its exit opening clips the beam); transmission counts particles through the housing. 0 = exit plane only")
parser.add_argument("--coast", type=int, default=300, help="steps a transmitted particle is followed past the exit plane")
parser.add_argument("--asym-steps", type=int, default=210, help="steps after the crossing for the asymptotic state (~55 mm)")
parser.add_argument("--dt", type=float, default=1.0e-10)
parser.add_argument("--seed", type=int, default=20260906)
parser.add_argument("--tag", default="scan2")
parser.add_argument("--rank", choices=("transmission", "vert"), default="transmission",
                    help="best point: highest transmission, or (vert) the smallest rms vertical angle vz/v_long of the transmitted bunch at the asymptotic state (vfom, mrad) among the points within --rank-tol of the best transmission")
parser.add_argument("--rank-tol", type=float, default=0.03, help="transmission tolerance for --rank vert (fraction)")
args = parser.parse_args()

out_dir = args.out_dir or args.reload_dir
os.makedirs(out_dir, exist_ok=True)
ti.STEP_DIR = args.step_dir
if args.particles:
    ti.PARTICLES = args.particles
if args.bfield:
    ti.BFIELD = args.bfield
need_skew = any(a != 0.0 for a in args.alpha1 + args.alpha2 + args.alpha3)
TRIPLET = args.q3 is not None
if TRIPLET and len(args.basis) < 7:
    raise SystemExit("--q3 needs 7 basis names: spiral q1 q1skew q2 q2skew q3 q3skew")

# ------------------------------------------------------------------ inputs
keys = ("spiral", "q1n", "q1s", "q2n", "q2s", "q3n", "q3s")[:len(args.basis)]
vals = {}
grid = None
for key, tag in zip(keys, args.basis):
    if key in ("q1s", "q2s", "q3s") and not need_skew:
        continue
    fn = os.path.join(args.reload_dir, "ef_itp_{}.pickle".format(tag))
    f = Field.from_file(fn)
    if grid is None:
        grid = f.grid
    vals[key] = f.grid_values
    for c in "xyz":
        if vals[key][c].shape != vals["spiral"][c].shape:
            raise RuntimeError("basis field {} is not on the same grid".format(tag))
    print("  basis {:7s} <- {}".format(key, os.path.basename(fn)), flush=True)

with open(os.path.join(args.reload_dir, "si_state_{}.pickle".format(args.basis[0])), "rb") as fh:
    state = pickle.load(fh)
trj, vdes = state["trj_design"], state["v_design"]
r_exit = float(np.linalg.norm(trj[-1][:2]))

r0, v0, ion, raw = ti.load_particles(args.n, seed=args.seed)
if args.swap_xy:
    r0[:, [0, 1]] = r0[:, [1, 0]]
    v0[:, [0, 1]] = v0[:, [1, 0]]
print("scan '{}': {:,d} particles per point from {}{}, B-field {}".format(
    args.tag, len(r0), os.path.basename(ti.PARTICLES), " (x<->y)" if args.swap_xy else "",
    os.path.basename(ti.BFIELD)), flush=True)

assembly = ti.load_assembly()
for uuid in [u for u, e in assembly.electrodes.items() if "assembly" in e.name.lower()]:
    assembly.electrodes.pop(uuid)
# a rigidly rotated system (bem_reload --rotate-all) records its angle in the design state;
# the collision geometry has to turn with the field.
rot_all = float(state.get("rotation_deg", 0.0) or 0.0)
if rot_all:
    rotate_assembly(assembly, rot_all)
    print("  assembly rotated {:+.1f} deg about z (rigid system rotation)".format(rot_all), flush=True)
shift_all = float(state.get("shift_all_m", 0.0) or 0.0)
if shift_all:
    shift_assembly(assembly, shift_all)
    print("  assembly shifted {:+.2f} mm along z".format(1e3 * shift_all), flush=True)
# the quads' collision geometry at the orientation of the basis solve (the scan's alpha1/2/3
# are extra rotations of the FIELD by superposition; the pole geometry stays at the basis orientation)
base_quad_rot = [float(a) for a in (state.get("quad_rotation_deg") or [])]
if any(a != 0.0 for a in base_quad_rot):
    rotate_quads(assembly, *base_quad_rot)
    print("  quad collision geometry rotated by {} deg (basis orientation)".format(" / ".join("{:g}".format(a) for a in base_quad_rot)), flush=True)
for e in assembly.electrodes.values():
    if e.generate_mesh() != 0:
        raise RuntimeError("failed to mesh {}".format(e.name))
index_map = {i: uuid for i, uuid in enumerate(assembly.electrodes.keys())}
housing_electrodes = [(i, e) for i, e in enumerate(assembly.electrodes.values()) if "housing" in e.name.lower()]
from spyral_inflector.tracking.deck import load_bfield as _load_bfield  # noqa: E402
bfield = _load_bfield(ti.BFIELD)


def quad_field(c, key_n, key_s, volt, alpha_deg):
    w = volt / args.unit
    if alpha_deg == 0.0 or key_s not in vals:
        return w * vals[key_n][c]
    a = np.radians(alpha_deg)
    return w * (np.cos(2 * a) * vals[key_n][c] + np.sin(2 * a) * vals[key_s][c])


def combined_field(q1, q2, a1, a2, s, q3=None, a3=0.0):
    values = {c: s * vals["spiral"][c] + quad_field(c, "q1n", "q1s", q1, a1) + quad_field(c, "q2n", "q2s", q2, a2)
              + (quad_field(c, "q3n", "q3s", q3, a3) if q3 is not None else 0.0)
              for c in "xyz"}
    return Field.from_arrays(grid=grid, values=values, dim=3, units="m",
                             label="spiral x{:.3f} q1={:.0f}@{:.0f} q2={:.0f}@{:.0f}".format(s, q1, a1, q2, a2)
                             + (" q3={:.0f}@{:.0f}".format(q3, a3) if q3 is not None else ""))


def rotated_beam(phi_deg):
    if phi_deg == 0.0:
        return r0.copy(), v0.copy()
    c, s = np.cos(np.radians(phi_deg)), np.sin(np.radians(phi_deg))
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return r0 @ R.T, v0 @ R.T


class _PD(object):
    pass


def run_point(q1, q2, a1, a2, s, phi, q3=None, a3=0.0):
    efield = combined_field(q1, q2, a1, a2, s, q3, a3)
    rr, vv = rotated_beam(phi)
    collision = ti.ElectrodeCollision(assembly)
    exit_plane = ti.ExitPlane(r_exit, point=trj[-1], normal=vdes[-1], coast_steps=max(args.coast, args.post_exit_steps),
                              asym_steps=min(args.asym_steps, max(args.coast, args.post_exit_steps)))
    collision.skip = exit_plane
    collision.post_exit_electrodes = housing_electrodes
    collision.post_exit_steps = args.post_exit_steps
    pd = _PD()
    pd.x_vec, pd.v_vec = rr, vv
    pd.alive = np.ones(len(rr), dtype=bool)
    pd.set_p_from_v_vec = lambda v: None
    t0 = time.time()
    Tracker(Pusher(ion, algorithm="rk4_rel"), efield, bfield,
            terminators=[exit_plane, collision]).run(pd, args.dt, args.nsteps, show_progress=False, sync_back=False)
    post = collision.post_exit_hit if collision.post_exit_hit is not None else np.zeros(len(rr), dtype=bool)
    crossed = exit_plane.crossed & ~post          # through the housing exit opening
    z_exit = exit_plane.state[crossed, 2]
    losses = collision.losses_by_electrode(index_map)
    if post.any():
        losses["Housing_exit"] = int(post.sum())
    spiral = sum(v for k, v in losses.items() if k in ("SI_Anode", "SI_Cathode"))
    quads = sum(v for k, v in losses.items() if k in ("D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7"))
    za = zr = aa = ar = vfom = zp_mean = zp_rms = None
    z_env = None
    if exit_plane.asym_state is not None:
        ok = crossed & np.all(np.isfinite(exit_plane.asym_state), axis=1)
        if ok.any():
            st_a = exit_plane.asym_state[ok]
            vnorm = np.linalg.norm(st_a[:, 3:], axis=1)
            va = np.degrees(np.arcsin(st_a[:, 5] / vnorm))
            za, zr = 1e3 * float(st_a[:, 2].mean()), 1e3 * float(st_a[:, 2].std())
            aa, ar = float(va.mean()), float(va.std())
            # vertical envelope past the exit, ballistic from the asymptotic state (taken
            # asym_steps after the crossing, i.e. s_a = asym_steps * dt * |v| along the path)
            s_a = exit_plane.asym_steps * args.dt * vnorm
            slope = st_a[:, 5] / vnorm
            env_s = np.array([0.0, 0.05, 0.10, 0.15, 0.20])
            z_env = [1e3 * float(np.std(st_a[:, 2] + (s - s_a) * slope)) for s in env_s]
            # the figure of merit: rms vertical angle vz / v_longitudinal (in-plane speed) [mrad]
            zp = st_a[:, 5] / np.hypot(st_a[:, 3], st_a[:, 4])
            zp_mean, zp_rms = 1e3 * float(zp.mean()), 1e3 * float(zp.std())
            vfom = zp_rms
    return {"q1": q1, "q2": q2, "alpha1": a1, "alpha2": a2, "vscale": s, "phi": phi, "q3": q3, "alpha3": a3,
            "transmission": float(crossed.sum()) / len(rr), "n_transmitted": int(crossed.sum()),
            "lost_spiral": spiral / len(rr), "lost_quads": quads / len(rr),
            "lost_apertures": (len(rr) - int(crossed.sum()) - spiral - quads) / len(rr),
            "z_exit_mean_mm": 1e3 * float(np.mean(z_exit)) if z_exit.size else None,
            "z_exit_rms_mm": 1e3 * float(np.std(z_exit)) if z_exit.size else None,
            "z_asym_mean_mm": za, "z_asym_rms_mm": zr, "vert_angle_asym_mean_deg": aa, "vert_angle_asym_rms_deg": ar,
            "z_env_s_mm": [0, 50, 100, 150, 200], "z_env_rms_mm": z_env,
            "z_env_mean_mm": None if z_env is None else float(np.mean(z_env[:4])),
            "zp_asym_mean_mrad": zp_mean, "zp_asym_rms_mrad": zp_rms,
            "vfom": vfom, "vfom_unit": "mrad (rms vz/v_long at the asymptotic state)",
            "losses": losses, "wall_s": time.time() - t0}


def pick_best(rows):
    """The best point of the scan so far: highest transmission, or with --rank vert the
    smallest vertical envelope among the points within --rank-tol of the best transmission."""
    rows = [r for r in rows if r["transmission"] > 0]
    if not rows:
        return None
    top = max(r["transmission"] for r in rows)
    if args.rank == "vert":
        cands = [r for r in rows if r["transmission"] >= top - args.rank_tol and r["vfom"] is not None]
        if cands:
            return min(cands, key=lambda r: (r["vfom"], -r["transmission"]))
    return max(rows, key=lambda r: r["transmission"])


if args.points:
    step = 3 if TRIPLET else 2
    pts = [tuple(args.points[i:i + step]) for i in range(0, len(args.points) - step + 1, step)]
else:
    q1s = np.linspace(args.q1[0], args.q1[1], int(args.q1[2]))
    q2s = np.linspace(args.q2[0], args.q2[1], int(args.q2[2]))
    if TRIPLET:
        q3s = np.linspace(args.q3[0], args.q3[1], int(args.q3[2]))
        pts = [(float(a), float(b), float(c)) for a in q1s for b in q2s for c in q3s]
    else:
        pts = [(float(a), float(b)) for a in q1s for b in q2s]

combos = list(itertools.product(args.phi, args.alpha1, args.alpha2, args.vscale, args.alpha3))
print("{} field/beam combinations x {} quad points = {} runs".format(len(combos), len(pts), len(combos) * len(pts)), flush=True)
print("{:>6s} {:>6s} {:>6s} {:>6s} {:>7s} {:>7s} {:>8s} {:>7s} {:>7s} {:>7s} {:>8s} {:>7s} {:>5s}   top losses".format(
    "phi", "alpha1", "alpha2", "Vsc", "q1[V]", "q2[V]", "transm.", "spiral", "quads", "apert.", "z_exit", "zp_rms", "s"))
results = []
best = None
for phi, a1, a2, s, a3 in combos:
    for pt in pts:
        q1, q2 = pt[0], pt[1]
        q3 = pt[2] if TRIPLET else None
        res = run_point(q1, q2, a1, a2, s, phi, q3, a3)
        results.append(res)
        top = sorted(res["losses"].items(), key=lambda kv: -kv[1])[:3]
        print("{:6.1f} {:6.1f} {:6.1f} {:6.3f} {:7.0f} {:7.0f} {:7.1f} % {:6.1f}% {:6.1f}% {:6.1f}% {:8s} {:7s} {:5.0f}   {}{}".format(
            phi, a1, a2, s, q1, q2, 100 * res["transmission"], 100 * res["lost_spiral"], 100 * res["lost_quads"],
            100 * res["lost_apertures"],
            "-" if res["z_exit_mean_mm"] is None else "{:+.1f} mm".format(res["z_exit_mean_mm"]),
            "-" if res["vfom"] is None else "{:.1f}mrad".format(res["vfom"]),
            res["wall_s"], ", ".join("{} {}".format(n, c) for n, c in top),
            "  | q3 {:+.0f}@{:g}".format(q3, a3) if TRIPLET else ""), flush=True)
        best = pick_best(results)
        with open(os.path.join(out_dir, "scan2_{}.json".format(args.tag)), "w") as fh:
            json.dump({"tag": args.tag, "n": len(r0), "unit_V": args.unit, "particles": ti.PARTICLES,
                       "swap_xy": bool(args.swap_xy), "bfield": ti.BFIELD, "basis": args.basis,
                       "results": results, "best": best}, fh, indent=2)

print("best ({}): phi {:.1f}, alpha {:.1f}/{:.1f}, spiral x{:.3f}, q1 = {:.0f} V, q2 = {:.0f} V{}: {:.1f} % "
      "(spiral {:.1f} %, quads {:.1f} %, apertures {:.1f} %), vertical angle {} mrad rms".format(
          args.rank, best["phi"], best["alpha1"], best["alpha2"], best["vscale"], best["q1"], best["q2"],
          ", q3 = {:.0f} V @ {:g}".format(best["q3"], best.get("alpha3", 0.0)) if TRIPLET else "",
          100 * best["transmission"], 100 * best["lost_spiral"], 100 * best["lost_quads"],
          100 * best["lost_apertures"], "-" if best["vfom"] is None else "{:.1f}".format(best["vfom"])), flush=True)

# best point as a regular field + state for BunchTrack.py (beam rotation is not part of the field:
# pass --phi to BunchTrack separately)
combined_field(best["q1"], best["q2"], best["alpha1"], best["alpha2"], best["vscale"], best.get("q3"), best.get("alpha3", 0.0)).save(
    os.path.join(out_dir, "ef_itp_{}_best.pickle".format(args.tag)))
volts = dict(state["electrode_voltages"])
for name, sgn in (("D0", 1), ("D1", 1), ("D2", -1), ("D3", -1)):
    volts[name] = sgn * best["q1"]
for name, sgn in (("D4", 1), ("D5", 1), ("D6", -1), ("D7", -1)):
    volts[name] = sgn * best["q2"]
if TRIPLET:
    for name, sgn in (("D8", 1), ("D9", 1), ("D10", -1), ("D11", -1)):
        volts[name] = sgn * best["q3"]
for name in ("SI_Anode", "SI_Cathode"):
    volts[name] = best["vscale"] * volts[name]
with open(os.path.join(out_dir, "si_state_{}_best.pickle".format(args.tag)), "wb") as fh:
    pickle.dump({"trj_design": trj, "v_design": vdes, "voltage": best["vscale"] * state["voltage"],
                 "electrode_voltages": volts,
                 "quad_rotation_deg": [b + a for b, a in zip(base_quad_rot + [0.0] * 3, [best["alpha1"], best["alpha2"], best.get("alpha3", 0.0)])][:max(len(base_quad_rot), 3 if TRIPLET else 2)],
                 "beam_rotation_deg": best["phi"], "rotation_deg": rot_all, "shift_all_m": shift_all}, fh)
print("wrote scan2_{0}.json, ef_itp_{0}_best.pickle, si_state_{0}_best.pickle".format(args.tag), flush=True)

# ------------------------------------------------------------------ plot: best transmission per combination
fig, ax = plt.subplots(figsize=(9, 5))
by_combo = {}
for r in results:
    k = (r["phi"], r["alpha1"], r["alpha2"], r["vscale"])
    if k not in by_combo or r["transmission"] > by_combo[k]["transmission"]:
        by_combo[k] = r
if len(args.phi) > 1:
    for (a1, a2, s) in sorted(set((k[1], k[2], k[3]) for k in by_combo)):
        ph = [k[0] for k in sorted(by_combo) if k[1:] == (a1, a2, s)]
        tr = [100 * by_combo[(p, a1, a2, s)]["transmission"] for p in ph]
        ax.plot(ph, tr, "o-", label="alpha {:.0f}/{:.0f}, V x{:.2f}".format(a1, a2, s))
    ax.set_xlabel("beam rotation phi [deg]")
else:
    labels = ["a{:.0f}/{:.0f}\nV{:.2f}".format(k[1], k[2], k[3]) for k in sorted(by_combo)]
    ax.bar(range(len(labels)), [100 * by_combo[k]["transmission"] for k in sorted(by_combo)])
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("best transmission over the quad grid [%]")
ax.set_title("scan '{}': {:,d} particles per point, best {:.1f} %".format(args.tag, len(r0), 100 * best["transmission"]))
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, "scan2_{}.png".format(args.tag)), dpi=140)
print("wrote scan2_{}.png".format(args.tag), flush=True)
