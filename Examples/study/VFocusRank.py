"""Re-rank the finished points of a VFocusScan campaign from their recorded quad grids, for
any transmission tolerance, without re-running anything: for each shape point the best
transmission on its grid, the smallest rms vertical angle within the tolerance (and the
quad setting that gives it), and the angle at the transmission optimum.

    python VFocusRank.py --name vfocus1 --tol 0.03 0.05 0.08
"""
import argparse
import glob
import json
import os

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--name", default="vfocus1")
p.add_argument("--tol", type=float, nargs="+", default=[0.03, 0.05, 0.08], help="transmission tolerances (fraction)")
p.add_argument("--write", action="store_true", help="write rerank.md next to the scan report")
p.add_argument("--write-winner", default=None, help="write the winner's knobs (at the last tolerance) to this json")
a = p.parse_args()

rows = []
for fn in sorted(glob.glob(os.path.join(DECK, "Results", "vfocus", a.name, "points", "*", "basis", "scan2_scan.json"))):
    name = os.path.basename(os.path.dirname(os.path.dirname(fn)))
    grid = [r for r in json.load(open(fn))["results"] if r.get("vfom") is not None]
    if not grid:
        continue
    t_max = max(r["transmission"] for r in grid)
    t_opt = max(grid, key=lambda r: r["transmission"])
    row = {"name": name, "t_max": t_max, "angle_at_t_opt": t_opt["vfom"], "q_t_opt": (t_opt["q1"], t_opt["q2"])}
    for tol in a.tol:
        cands = [r for r in grid if r["transmission"] >= t_max - tol]
        b = min(cands, key=lambda r: (r["vfom"], -r["transmission"]))
        row[tol] = (b["vfom"], b["transmission"], b["q1"], b["q2"])
    rows.append(row)

hdr = "| point | best T | angle at T-opt [mrad] | " + " | ".join("tol {:.0f} pts: angle / T / q1,q2".format(100 * t) for t in a.tol) + " |"
L = ["# {}: re-ranking by the rms vertical angle for several transmission tolerances".format(a.name), "",
     "Per shape point: the best transmission on its 4x4 quad grid, the angle at that setting, and for each tolerance "
     "the smallest angle among the settings within the tolerance of the best transmission.", "", hdr,
     "|---|---|---|" + "---|" * len(a.tol)]
for r in sorted(rows, key=lambda r: r[a.tol[-1]][0]):
    L.append("| {} | {:.1f} % | {:.0f} | ".format(r["name"], 100 * r["t_max"], r["angle_at_t_opt"]) +
             " | ".join("{:.0f} / {:.1f} % / {:+.0f},{:+.0f}".format(r[t][0], 100 * r[t][1], r[t][2], r[t][3]) for t in a.tol) + " |")
for t in a.tol:
    best = min(rows, key=lambda r: (r[t][0], -r[t][1]))
    L.append("")
    L.append("tolerance {:.0f} points: winner {} at {:.0f} mrad, {:.1f} % (q {:+.0f}/{:+.0f})".format(
        100 * t, best["name"], best[t][0], 100 * best[t][1], best[t][2], best[t][3]))
print("\n".join(L))
if a.write:
    with open(os.path.join(DECK, "Results", "vfocus", a.name, "rerank.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
if a.write_winner and rows:
    t = a.tol[-1]
    best = min(rows, key=lambda r: (r[t][0], -r[t][1]))
    s_mm, ang, gam = best["name"].split("_")            # s1.2_a11_g5
    knobs = {"sigma": 1e-3 * float(s_mm[1:]), "angling": float(ang[1:]), "gamma": float(gam[1:]),
             "_winner": best["name"], "_tolerance": t, "_angle_mrad": best[t][0], "_transmission": best[t][1],
             "_q1": best[t][2], "_q2": best[t][3]}
    with open(a.write_winner, "w") as fh:
        json.dump(knobs, fh, indent=2)
    print("winner knobs ->", a.write_winner)
