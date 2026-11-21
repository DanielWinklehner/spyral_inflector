"""Where to centre the housing EXIT opening on the beam (Daniel, 2026-10-20: "use the handoff distribution to find the
center of the aperture").

The housing cuts its exit opening along the exit direction n (horizontal part of the design velocity at the electrode
exit), centred at c0 = trj[-1] + 0.99 * bottom_distance * (vx, vy, 0)/|v_exit| -- at the orbit's exit HEIGHT, with no
climb (SIHousing.create_geo_str; the script checks c0 against the opening the build recorded). This takes every particle of a run's hand-off (the electrode-exit plane, tracking frame),
carries it on a straight line to the mid-plane of the exit plate (exit_plane.json: outer face minus half the plate
thickness, ~4.5 mm on), and returns the beam centroid there as the offset of the opening's centre from c0:
du along the plate (t = z x n) and dw vertically -- the knob exit_circle = [radius, du, dw] of build_geometry.
Both are invariant under the rigid rotation R, so the tracking-frame numbers are the deck-frame knob values.

    python exit_hole_center.py <run>   (uses <run>/sc_level/bunch_sc_all_level.npz, basis/si_state_spiral.pickle,
                                        exit_plane.json, steps/housing_outline.txt, geometry/summary.json)
"""
import argparse
import json
import os

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("run")
p.add_argument("--npz", default=None, help="default: <run>/sc_level/bunch_sc_all_level.npz")
p.add_argument("--out", default=None, help="json with the result (default: <run>/exit_hole_center.json)")
a = p.parse_args()

from spyral_inflector.tracking.deck import load_state  # noqa: E402
from spyral_inflector.tracking.export import HOUSING_OUTLINE_NAME, read_housing_outline  # noqa: E402

run = os.path.abspath(a.run)
npz = a.npz or os.path.join(run, "sc_level", "bunch_sc_all_level.npz")
st = load_state(os.path.join(run, "basis", "si_state_spiral.pickle"))
ep = json.load(open(os.path.join(run, "exit_plane.json")))
summ = json.load(open(os.path.join(run, "geometry", "summary.json")))
R = float(st.get("rotation_deg", 0.0) or 0.0)
if abs(R - float(ep["rotation_deg"])) > 1e-6:
    raise SystemExit("basis state rotation {} != exit_plane.json {}".format(R, ep["rotation_deg"]))

trj, vdes = np.asarray(st["trj_design"], float), np.asarray(st["v_design"], float)
v3 = vdes[-1] / np.linalg.norm(vdes[-1])
n2 = v3[:2] / np.linalg.norm(v3[:2])
n3 = np.array([n2[0], n2[1], 0.0])
t3 = np.array([-n2[1], n2[0], 0.0])
z3 = np.array([0.0, 0.0, 1.0])
b_gap = float(summ["knobs"]["bottom_distance"])
c0 = trj[-1] + np.array([v3[0], v3[1], 0.0]) * b_gap * 0.99    # the housing keeps z = trj[-1][2] (no climb)

# check: c0 must be the centre of the opening the build recorded (deck frame, unrotated -> rotate by R)
ol = read_housing_outline(os.path.join(run, "steps", HOUSING_OUTLINE_NAME))
c, s = np.cos(np.radians(R)), np.sin(np.radians(R))
ap = np.asarray(ol["aperture"], float)
ap_c = ap.mean(axis=0)
ap_c = np.array([c * ap_c[0] - s * ap_c[1], s * ap_c[0] + c * ap_c[1], ap_c[2]])
print("default opening centre c0 (tracking frame) {} mm; recorded opening centre {} mm; difference {:.4f} mm".format(
    np.round(1e3 * c0, 3), np.round(1e3 * ap_c, 3), 1e3 * np.linalg.norm(c0 - ap_c)))

# plate mid-plane (tracking frame = the spec's 'deck' entry, which is rotated by R)
face = 1e-3 * np.asarray(ep["deck"]["point_mm"], float)
N = np.asarray(ep["deck"]["normal"], float)
thick = 1e-3 * float(ep["plate_thickness_mm"])
mid = face - 0.5 * thick * N
print("plate: outer face {:.2f} mm past the hand-off plane along n, thickness {:.2f} mm, normal vs n {:.3f} deg".format(
    1e3 * (face - trj[-1]) @ n3, 1e3 * thick, np.degrees(np.arccos(min(1.0, abs(N @ n3))))))

d = np.load(npz, allow_pickle=True)
hs = d["handoff_state"]
ok = d["crossed"] & np.all(np.isfinite(hs), axis=1)
r, v = hs[ok, :3], hs[ok, 3:]
tail = d["is_tail"][ok] if "is_tail" in d.files else np.zeros(len(r), bool)
# straight line to the plate mid-plane
sgn = (r - mid) @ N
r_mid = r - v * (sgn / (v @ N))[:, None]


def cen(rr):
    rel = rr - c0
    return float(np.mean(rel @ t3)), float(np.mean(rel @ z3))


u_h, w_h = cen(r)
u_m, w_m = cen(r_mid)
vz = v[:, 2] / np.hypot(v[:, 0], v[:, 1])
vu = (v @ t3) / (v @ n3)
des_z = v3[2] / np.linalg.norm(v3[:2])
print("hand-off: {:,d} particles ({:,d} tail) of {}".format(int(ok.sum()), int(tail.sum()), os.path.basename(npz)))
print("  centroid vs the default centre line, at the hand-off plane: du {:+.3f} mm, dw {:+.3f} mm (deck z)".format(1e3 * u_h, 1e3 * w_h))
print("  centroid angles vs the hole axis n: horizontal {:+.2f} mrad, vertical {:+.2f} mrad (design orbit {:+.2f} mrad)".format(
    1e3 * vu.mean(), 1e3 * vz.mean(), 1e3 * des_z))
print("  at the plate mid-plane: du {:+.3f} mm, dw {:+.3f} mm (deck z) = {:+.3f} mm in the Baseline frame".format(
    1e3 * u_m, 1e3 * w_m, -1e3 * w_m))
rel = r_mid - c0
rad = np.hypot(rel @ t3 - u_m, rel @ z3 - w_m)
for rr in (8e-3, 10e-3, 12e-3):
    print("  share of the hand-off particles outside r = {:.0f} mm of that centre at the plate: {:.2f} %".format(1e3 * rr, 100 * (rad > rr).mean()))
out = {"run": run, "npz": npz, "R_deg": R, "c0_m": c0.tolist(), "du_m": u_m, "dw_m": w_m, "du_handoff_m": u_h, "dw_handoff_m": w_h,
       "n_particles": int(ok.sum()), "n_tail": int(tail.sum()), "plate_mid_m": mid.tolist(),
       "frac_outside": {"8mm": float((rad > 8e-3).mean()), "10mm": float((rad > 10e-3).mean()), "12mm": float((rad > 12e-3).mean())}}
with open(a.out or os.path.join(run, "exit_hole_center.json"), "w") as fh:
    json.dump(out, fh, indent=2)
