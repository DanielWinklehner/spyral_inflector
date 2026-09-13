"""The design exit radius of every FullOpt evaluation, against k' and the objective.

The exit radius is not part of the objective: the optimizer chases the vertical angle and the
radius follows k' (about -0.6 mm per degree on the doublet). The central region is re-optimized
for whatever radius comes out, so this table is the record of what was traded, and it works on a
run in progress as well as a finished one.

    python FullOptRadius.py --run <deck>\\Results\\final\\vfocus1_level

Writes <run>/full_opt_radius.md and prints the table.
"""
import argparse
import json
import math
import os

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run", required=True)
p.add_argument("--out", default=None, help="default: <run>/full_opt_radius.md")
a = p.parse_args()

RUN = os.path.abspath(a.run)
OPT = os.path.join(RUN, "full_opt")
evals = json.load(open(os.path.join(OPT, "evals.json")))
NAMES = ["tilt", "angling", "gamma", "sigma_mm", "t_ent", "t_exit", "dz_mm", "q1", "q2", "gap_mm"]

rows = []
for e in evals:
    g = os.path.join(OPT, e["tag"], "geometry", "summary.json")
    if not os.path.exists(g):
        continue
    pt = json.load(open(g))["optimizer"].get("exit_point_mm")
    if not pt:
        continue
    m = e["measure"] or {}
    rows.append({"tag": e["tag"], "tilt": e["x"][0], "gap": e["x"][9], "r": math.hypot(pt[0], pt[1]), "z": pt[2],
                 "angle": 1e3 * m.get("zp_rms", float("nan")), "T": 100 * m.get("transmission", float("nan")),
                 "f": e["objective"], "x": e["x"]})
if not rows:
    raise SystemExit("no evaluation with a geometry summary in " + OPT)

best = min(rows, key=lambda r: r["f"])
seed = rows[0]
lo, hi = min(rows, key=lambda r: r["tilt"]), max(rows, key=lambda r: r["tilt"])
slope = (hi["r"] - lo["r"]) / (hi["tilt"] - lo["tilt"]) if hi["tilt"] != lo["tilt"] else float("nan")

L = ["# Design exit radius across the FullOpt evaluations: {}".format(os.path.basename(RUN)), "",
     "The radius is not in the objective; it follows k'. Over the range probed so far "
     "({:.2f} to {:.2f} deg) the slope is **{:+.2f} mm per degree of k'**. Seed: k' {:.2f} deg, r = {:.2f} mm. "
     "Best objective so far ({}): k' {:.2f} deg, **r = {:.2f} mm**, exit height {:.2f} mm, {:.1f} mrad at {:.1f} %."
     .format(lo["tilt"], hi["tilt"], slope, seed["tilt"], seed["r"], best["tag"], best["tilt"], best["r"], best["z"],
             best["angle"], best["T"]),
     "", "| tag | k' [deg] | gap [mm] | exit radius [mm] | exit height [mm] | rms vert. angle [mrad] | T [%] | objective |",
     "|---|---|---|---|---|---|---|---|"]
for r in rows:
    mark = " **" if r is best else ""
    L.append("| {} | {:.2f} | {:.2f} | {:.2f}{} | {:.2f} | {:.1f} | {:.1f} | {:.2f} |".format(
        r["tag"], r["tilt"], r["gap"], r["r"], mark, r["z"], r["angle"], r["T"], r["f"]))
L += ["", "Reference points: Jarrett's deck 71.0 mm, final1 71.8 mm, vfocus1_level (the seed of this run) {:.1f} mm."
      .format(seed["r"])]
out = a.out or os.path.join(RUN, "full_opt_radius.md")
with open(out, "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
print("\n".join(L[4:]))
print("\nwrote", out)
