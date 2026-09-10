"""Does the rotated inflector housing sit inside the pole cutouts of a field map?

The map tells where the iron is (|B| ~ 2 T inside, < 1.2 T in the gap and the cutouts).
Every vertex of the housing, the cylinder and the electrodes is rotated by the run's
solved R and the map is sampled there; vertices in iron are reported with their location.
Also prints the cutout half-width along +-x and the pole-hole radius at a few heights.

    python HousingFitCheck.py --run Results/final/final1 [--iron 1.4]
"""
import argparse
import json
import os

import numpy as np

from spyral_inflector.tracking.deck import load_step_assembly, mesh_assembly, load_bfield
from spyral_inflector.tracking.exit_plane import rotz

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)

p = argparse.ArgumentParser()
p.add_argument("--run", required=True, help="Results/final/<name>")
p.add_argument("--bfield", default=None, help="default: the run's phase-1 pickle")
p.add_argument("--iron", type=float, default=1.4, help="|B| above this [T] counts as iron")
a = p.parse_args()

run = a.run if os.path.isabs(a.run) else os.path.join(DECK, a.run)
spec = json.load(open(os.path.join(run, "exit_plane.json")))
R = spec["rotation_deg"]
bf_path = a.bfield or json.load(open(os.path.join(run, "phase_1.done")))["info"]["pickle"]
bf = load_bfield(bf_path)                       # deck frame (mirrored in memory if the map is Baseline)
Rz = rotz(R)

assy = load_step_assembly(os.path.join(run, "steps"))
mesh_assembly(assy)
print("=" * 96)
print("HOUSING FIT CHECK  run {}  R = {:+.3f} deg  map {}  (iron: |B| > {:.1f} T)".format(os.path.basename(run), R, os.path.basename(bf_path), a.iron))
print("=" * 96)
print("{:20s} {:>8s} {:>8s}   {}".format("electrode", "vertices", "in iron", "where (BCS frame: r mm, azimuth deg, z mm)"))
worst = []
for e in assy.electrodes.values():
    V = np.asarray(e._gmsh_msh["vertices"], float) @ Rz.T          # rotated into the deployed orientation, deck frame
    B = np.linalg.norm(bf(V), axis=1)
    bad = B > a.iron
    where = ""
    if bad.any():
        Vb = V[bad]
        r = np.hypot(Vb[:, 0], Vb[:, 1])
        az = np.degrees(np.arctan2(Vb[:, 1], Vb[:, 0]))
        z = -Vb[:, 2]                                               # BCS z (mirror of the deck)
        where = "r {:.0f}..{:.0f}, az {:+.0f}..{:+.0f}, BCS z {:+.0f}..{:+.0f}".format(r.min() * 1e3, r.max() * 1e3, az.min(), az.max(), z.min() * 1e3, z.max() * 1e3)
        worst.append((e.name, int(bad.sum())))
    print("{:20s} {:8d} {:8d}   {}".format(e.name, len(V), int(bad.sum()), where))

# the cutout as the map shows it: along +x at a few BCS heights, the largest |y| still free of iron at r = 60..100 mm
print()
print("cutout along +x as the map shows it (free of iron): half-width in y at r = 60, 80, 100 mm")
for zb in (0.0, 0.02, 0.04, 0.05, 0.06, 0.07):
    row = []
    for r in (0.06, 0.08, 0.10):
        ys = np.linspace(0, 0.14, 141)
        pts = np.column_stack([np.full_like(ys, r), ys, np.full_like(ys, -zb)])   # deck z = -BCS z
        B = np.linalg.norm(bf(pts), axis=1)
        k = np.argmax(B > a.iron) if (B > a.iron).any() else len(ys)
        row.append("{:5.0f} mm".format(1e3 * (ys[k - 1] if 0 < k < len(ys) else (ys[-1] if k == len(ys) else 0.0))))
    print("   BCS z = {:+.0f} mm: {}".format(1e3 * zb, "   ".join(row)))
print()
print("pole hole radius above the cutout (iron-free r along +-y) at BCS z = 65, 80, 100 mm:")
for zb in (0.065, 0.08, 0.10):
    rs = np.linspace(0, 0.14, 281)
    pts = np.column_stack([np.zeros_like(rs), rs, np.full_like(rs, -zb)])
    B = np.linalg.norm(bf(pts), axis=1)
    k = np.argmax(B > a.iron) if (B > a.iron).any() else len(rs)
    print("   BCS z = {:+.0f} mm: r_free = {:.0f} mm".format(1e3 * zb, 1e3 * (rs[k - 1] if 0 < k < len(rs) else rs[-1])))
print("=" * 96)
print("VERDICT: " + ("no housing/electrode vertex in iron" if not worst else "IN IRON: " + ", ".join("{} ({})".format(n, c) for n, c in worst)))
