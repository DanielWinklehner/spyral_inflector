"""Refine the quad setting of the best shapes of a VFocusScan campaign with a finer q2 grid.

The coarse 4x4 grids (1000 V steps) locate the region of small rms vertical angle but the
minimum is narrow in q2: some shapes land on it, others between grid points. For the
--top shapes (ranked from the coarse grids at --tol), the basis fields are reused and a
finer grid (q1 --q1, q2 --q2, --n particles) is tracked by superposition; each shape's
setting is then the smallest rms vertical angle within --tol of its best transmission over
coarse + fine grids, and the shapes are ranked the same way.

    python VFocusRefine.py --name vfocus1 --top 4 --tol 0.05
Writes rerank_refined.md, refined.json and winner_knobs_refined.json in the campaign folder.
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
PY = sys.executable

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--name", default="vfocus1")
p.add_argument("--top", type=int, default=4)
p.add_argument("--tol", type=float, default=0.05)
p.add_argument("--q1", type=float, nargs=3, default=[6075, 7075, 3], metavar=("MIN", "MAX", "N"))
p.add_argument("--q2", type=float, nargs=3, default=[-7200, -5700, 7], metavar=("MIN", "MAX", "N"))
p.add_argument("--n", type=int, default=1500)
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle"))
a = p.parse_args()

OUT = os.path.join(DECK, "Results", "vfocus", a.name)
LOG = open(os.path.join(OUT, "log_refine.txt"), "a", encoding="utf-8")


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def pick(grid, tol):
    grid = [r for r in grid if r.get("vfom") is not None]
    t_max = max(r["transmission"] for r in grid)
    cands = [r for r in grid if r["transmission"] >= t_max - tol]
    return min(cands, key=lambda r: (r["vfom"], -r["transmission"])), t_max


# ---------------------------------------------------------------- the coarse ranking
shapes = []
for fn in sorted(glob.glob(os.path.join(OUT, "points", "*", "basis", "scan2_scan.json"))):
    pdir = os.path.dirname(os.path.dirname(fn))
    grid = json.load(open(fn))["results"]
    best, t_max = pick(grid, a.tol)
    shapes.append({"name": os.path.basename(pdir), "dir": pdir, "coarse": grid, "coarse_best": best, "coarse_t_max": t_max})
shapes.sort(key=lambda s: (s["coarse_best"]["vfom"], -s["coarse_best"]["transmission"]))
top = shapes[:a.top]
stamp("REFINE '{}': {} shapes, refining the top {}: {}".format(
    a.name, len(shapes), len(top), ", ".join("{} ({:.0f} mrad, {:.1f} %)".format(s["name"], s["coarse_best"]["vfom"], 100 * s["coarse_best"]["transmission"]) for s in top)))

# ---------------------------------------------------------------- the fine grids
for s in top:
    pdir, basis = s["dir"], os.path.join(s["dir"], "basis")
    fine_fn = os.path.join(basis, "scan2_fine.json")
    if not os.path.exists(fine_fn):
        R = json.load(open(os.path.join(pdir, "exit_plane.json")))["rotation_deg"]
        stamp("  {}: fine grid {}x{} at {} particles (R {:+.3f})".format(s["name"], int(a.q1[2]), int(a.q2[2]), a.n, R))
        cmd = [PY, "-u", os.path.join(SCRIPTS, "BunchScan2.py"), "--reload-dir", basis, "--out-dir", basis,
               "--step-dir", os.path.join(pdir, "steps"), "--bfield", a.bfield, "--basis", "spiral", "q1", "q1", "q2", "q2",
               "--phi", str(90.0 + R), "--particles", a.particles, "--n", str(a.n),
               "--q1", str(a.q1[0]), str(a.q1[1]), str(int(a.q1[2])), "--q2", str(a.q2[0]), str(a.q2[1]), str(int(a.q2[2])),
               "--tag", "fine", "--rank", "vert", "--rank-tol", str(a.tol)]
        with open(os.path.join(pdir, "log_fine.txt"), "w", encoding="utf-8") as fh:
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS)
        if r.returncode != 0:
            stamp("  {}: fine grid FAILED (exit {})".format(s["name"], r.returncode))
            continue
    fine = json.load(open(fine_fn))["results"]
    s["fine"] = fine
    s["best"], s["t_max"] = pick(s["coarse"] + fine, a.tol)
    s["fine_best"], s["fine_t_max"] = pick(fine, a.tol)
    stamp("  {}: refined {:.1f} mrad at {:.1f} % (q {:+.0f}/{:+.0f}); coarse pick was {:.1f} mrad at {:.1f} %".format(
        s["name"], s["best"]["vfom"], 100 * s["best"]["transmission"], s["best"]["q1"], s["best"]["q2"],
        s["coarse_best"]["vfom"], 100 * s["coarse_best"]["transmission"]))

done = [s for s in top if "best" in s]
done.sort(key=lambda s: (s["best"]["vfom"], -s["best"]["transmission"]))
L = ["# {}: refinement of the top {} shapes (q2 in {:.0f} V steps, {} particles, tolerance {:.0f} points)".format(
        a.name, len(top), (a.q2[1] - a.q2[0]) / max(1, int(a.q2[2]) - 1), a.n, 100 * a.tol), "",
     "| rank | shape | coarse pick: angle / T | refined: angle / T / q1,q2 | best T on all grids |", "|---|---|---|---|---|"]
for i, s in enumerate(done):
    L.append("| {} | {} | {:.1f} mrad / {:.1f} % | **{:.1f} mrad** / {:.1f} % / {:+.0f},{:+.0f} | {:.1f} % |".format(
        i + 1, s["name"], s["coarse_best"]["vfom"], 100 * s["coarse_best"]["transmission"], s["best"]["vfom"],
        100 * s["best"]["transmission"], s["best"]["q1"], s["best"]["q2"], 100 * s["t_max"]))
if done:
    w = done[0]
    L += ["", "Winner: **{}** at {:.1f} mrad rms, {:.1f} %, q1 {:+.0f} / q2 {:+.0f} V.".format(
        w["name"], w["best"]["vfom"], 100 * w["best"]["transmission"], w["best"]["q1"], w["best"]["q2"])]
    s_mm, ang, gam = w["name"].split("_")
    with open(os.path.join(OUT, "winner_knobs_refined.json"), "w") as fh:
        json.dump({"sigma": 1e-3 * float(s_mm[1:]), "angling": float(ang[1:]), "gamma": float(gam[1:]),
                   "_winner": w["name"], "_tolerance": a.tol, "_angle_mrad": w["best"]["vfom"], "_transmission": w["best"]["transmission"],
                   "_q1": w["best"]["q1"], "_q2": w["best"]["q2"]}, fh, indent=2)
with open(os.path.join(OUT, "rerank_refined.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
with open(os.path.join(OUT, "refined.json"), "w") as fh:
    json.dump([{k: v for k, v in s.items() if k not in ("coarse", "fine")} for s in done], fh, indent=2)
print("\n".join(L))
stamp("refine done: winner {}".format(done[0]["name"] if done else "-"))
