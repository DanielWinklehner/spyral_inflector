"""Validate a new B-field pickle against the previous one: agreement in the overlap
region, the on-axis profile down the bore, the design particle's exit point, and a
2,000-particle bunch through the optimum E-field on both maps.

    python validate_bfield.py --new ../Fields/HCHC-60_CentralBField_z-40to5cm_1mm.pickle
"""
import argparse
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import track_inflector as ti  # noqa: E402
from PyPATools.field import Field  # noqa: E402
from PyPATools.pusher import Pusher  # noqa: E402
from PyPATools.trackers import Tracker  # noqa: E402

DECK = os.environ.get("SI_DECK") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
parser = argparse.ArgumentParser()
parser.add_argument("--new", required=True)
parser.add_argument("--old", default=ti.BFIELD)
parser.add_argument("--tag", default="fine_best")
parser.add_argument("--reload-dir", default=os.path.join(DECK, "Results", "reload"))
parser.add_argument("--particles", default=None)
parser.add_argument("--swap-xy", action="store_true")
parser.add_argument("--n", type=int, default=2000)
args = parser.parse_args()
if args.particles:
    ti.PARTICLES = args.particles

new = Field.from_file(args.new)
old = Field.from_file(args.old)
gn, vn = new.grid, new.grid_values
go, vo = old.grid, old.grid_values
print("new: x {:+.3f}..{:+.3f} y {:+.3f}..{:+.3f} z {:+.3f}..{:+.3f} ({} x {} x {})".format(
    gn["x"][0], gn["x"][-1], gn["y"][0], gn["y"][-1], gn["z"][0], gn["z"][-1], len(gn["x"]), len(gn["y"]), len(gn["z"])))
print("old: x {:+.3f}..{:+.3f} y {:+.3f}..{:+.3f} z {:+.3f}..{:+.3f} ({} x {} x {})".format(
    go["x"][0], go["x"][-1], go["y"][0], go["y"][-1], go["z"][0], go["z"][-1], len(go["x"]), len(go["y"]), len(go["z"])))

# overlap on common grid points (both 1 mm grids containing the origin)
def common(a, b):
    ia = {round(v, 6): i for i, v in enumerate(a)}
    ib = {round(v, 6): i for i, v in enumerate(b)}
    keys = sorted(set(ia) & set(ib))
    return np.array([ia[k] for k in keys]), np.array([ib[k] for k in keys]), np.array(keys)

ixn, ixo, xs = common(gn["x"], go["x"])
iyn, iyo, ys = common(gn["y"], go["y"])
izn, izo, zs = common(gn["z"], go["z"])
print("overlap: {} x {} x {} points, z {:+.3f}..{:+.3f}".format(len(xs), len(ys), len(zs), zs[0], zs[-1]))
X, Y = np.meshgrid(xs, ys, indexing="ij")
R = np.hypot(X, Y)
for comp in "xyz":
    a = vn[comp][np.ix_(ixn, iyn, izn)]
    b = vo[comp][np.ix_(ixo, iyo, izo)]
    d = np.abs(a - b)
    for rmax in (0.03, 0.10, 0.15):
        m = np.broadcast_to((R < rmax)[:, :, None], d.shape)
        dd = d[m]
        print("  B{} |new - old| inside r < {:.0f} cm: max {:.3f} mT, 99.99 % {:.4f} mT, median {:.5f} mT".format(
            comp, 100 * rmax, 1e3 * dd.max(), 1e3 * np.percentile(dd, 99.99), 1e3 * np.median(dd)))

print("\non-axis profile of the new map (deck frame, beam travels +z):")
zz = np.arange(gn["z"][0], gn["z"][-1] + 1e-6, 0.025)
b = new(np.column_stack([np.zeros_like(zz), np.zeros_like(zz), zz]))
for z_, bz in zip(zz, b[:, 2]):
    print("    z = {:+.3f} m: Bz = {:+.5f} T".format(z_, bz))
# radial field at r = 10 mm (the solenoidal focusing in the bore)
zz2 = np.arange(-0.40, 0.0, 0.05)
br = new(np.column_stack([np.full_like(zz2, 0.01), np.zeros_like(zz2), zz2]))[:, 0]
print("  Bx at x = 10 mm: " + ", ".join("z {:+.2f}: {:+.2f} mT".format(z_, 1e3 * b_) for z_, b_ in zip(zz2, br)))

# tracking: design particle and a bunch, old vs new map
efield = Field.from_file(os.path.join(args.reload_dir, "ef_itp_{}.pickle".format(args.tag)))
with open(os.path.join(args.reload_dir, "si_state_{}.pickle".format(args.tag)), "rb") as fh:
    state = pickle.load(fh)
trj, vdes = state["trj_design"], state["v_design"]
r_exit = float(np.linalg.norm(trj[-1][:2]))
assembly = ti.load_assembly()
for uuid in [u for u, e in assembly.electrodes.items() if "assembly" in e.name.lower()]:
    assembly.electrodes.pop(uuid)
for e in assembly.electrodes.values():
    e.generate_mesh()
index_map = {i: uuid for i, uuid in enumerate(assembly.electrodes.keys())}
ion = ti.IonSpecies(ti.SPECIES)
r0b, v0b, _, _ = ti.load_particles(args.n)
if args.swap_xy:
    r0b[:, [0, 1]] = r0b[:, [1, 0]]
    v0b[:, [0, 1]] = v0b[:, [1, 0]]
gam = 1 + 1e-3 * 69.34 / ion.mass_mev
r0s = np.array([[0.0, 0.0, -0.130], [0.0, 0.0, -0.288]])
v0s = np.array([[0.0, 0.0, np.sqrt(1 - 1 / gam ** 2) * ti.CLIGHT]] * 2)


class _PD(object):
    pass


def run(bf, r0, v0, nsteps):
    pd = _PD()
    pd.x_vec, pd.v_vec = r0.copy(), v0.copy()
    pd.alive = np.ones(len(r0), dtype=bool)
    pd.set_p_from_v_vec = lambda v: None
    col = ti.ElectrodeCollision(assembly)
    ep = ti.ExitPlane(r_exit, point=trj[-1], normal=vdes[-1])
    Tracker(Pusher(ion, algorithm="rk4_rel"), efield, bf, terminators=[ep, col]).run(
        pd, 1e-10, nsteps, show_progress=False, sync_back=False)
    return ep, col


res = {}
for name, bf in (("old", ti.with_fast_interpolator(old, "old")), ("new", ti.with_fast_interpolator(new, "new"))):
    t0 = time.time()
    ep, _ = run(bf, r0s, v0s, 2200)
    ep2, col2 = run(bf, r0b, v0b, 1800)
    losses = col2.losses_by_electrode(index_map)
    res[name] = (ep.state[:, :3].copy(), ep.crossed.copy(), ep2.crossed.mean(), ep2.state[ep2.crossed, 2], losses)
    print("\n{} map ({:.0f} s): on-axis 69.3 keV from z = -130 mm -> exit {} ; from z = -288 mm (through the quads) -> {}".format(
        name, time.time() - t0,
        np.round(1e3 * ep.state[0, :3], 2) if ep.crossed[0] else "lost",
        np.round(1e3 * ep.state[1, :3], 2) if ep.crossed[1] else "lost"))
    print("   bunch of {}: transmission {:.2f} %, exit z mean {:+.2f} mm; top losses {}".format(
        len(r0b), 100 * res[name][2], 1e3 * res[name][3].mean(),
        sorted(losses.items(), key=lambda kv: -kv[1])[:4]))
for i, lab in enumerate(("from -130 mm", "from -288 mm")):
    if res["old"][1][i] and res["new"][1][i]:
        print("exit-point shift new vs old ({}): {:.3f} mm".format(lab, 1e3 * np.linalg.norm(res["new"][0][i] - res["old"][0][i])))
