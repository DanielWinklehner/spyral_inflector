"""Doublet quad-rotation scan by superposition on a final-push run (follow-up 3 of the vfocus
campaign): solves the two skew basis fields of the run (each quad rotated by a further 45 deg
on the run's own basis grid) if they are missing, then scans extra rotations alpha1 x alpha2
of the two quads at a fixed voltage pair on a bunch, ranked by the rms vertical angle vz/v_long
at 55 mm past the exit within --rank-tol of the best transmission. The (0, 0) point is the run's
own setting. Cheap test of the "phase-space orientation" hypothesis before believing a triplet.

    python DoubletRotScan.py --run <deck>\\Results\\final\\vfocus1_final --q 6575 -6700
    python DoubletRotScan.py --run ... --alpha1 -30 -20 -10 0 10 20 30 --alpha2 -30 -20 -10 0 10 20 30

Writes <run>/basis/ef_itp_q1skew.pickle, ef_itp_q2skew.pickle (bem_reload), <run>/basis/scan2_rot.json
(BunchScan2) and <run>/rot_scan.md.
"""
import argparse
import json
import os
import subprocess
import sys
import time

DECK = os.environ.get("SI_DECK", r"D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron\Spiral_inflector")
SCRIPTS = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run", required=True, help="final-push run folder (geometry/, steps/, basis/, exit_plane.json)")
p.add_argument("--q", type=float, nargs=2, required=True, metavar=("Q1", "Q2"), help="fixed quad voltages [V]")
p.add_argument("--alpha1", type=float, nargs="+", default=[-30, -20, -10, 0, 10, 20, 30], help="extra rotations of quad 1 [deg]")
p.add_argument("--alpha2", type=float, nargs="+", default=[-30, -20, -10, 0, 10, 20, 30], help="extra rotations of quad 2 [deg]")
p.add_argument("--quads", type=float, nargs=2, default=[16.0, 24.0], help="basis rotation of the quads on top of R [deg]")
p.add_argument("--res", type=float, default=0.0025, help="basis-field resolution [m] (the run's)")
p.add_argument("--n", type=int, default=1000)
p.add_argument("--rank-tol", type=float, default=0.05)
p.add_argument("--tag", default="rot")
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle"))
a = p.parse_args()

RUN = os.path.abspath(a.run)
GEO, STEPS, BASIS = (os.path.join(RUN, d) for d in ("geometry", "steps", "basis"))
LOG = open(os.path.join(RUN, "log_rot_scan.txt"), "a", encoding="utf-8")


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def run(args, logname):
    with open(os.path.join(RUN, logname), "w", encoding="utf-8") as fh:
        rc = subprocess.call([PY, "-u"] + [str(x) for x in args], cwd=SCRIPTS, stdout=fh, stderr=subprocess.STDOUT)
    if rc != 0:
        raise SystemExit("FAILED (exit {}): {} -> see {}".format(rc, args[:3], logname))


sys.path.insert(0, os.path.dirname(os.path.dirname(SCRIPTS)))
from spyral_inflector.tracking.deck import particle_rows  # noqa: E402

with open(os.path.join(RUN, "exit_plane.json")) as fh:
    R = json.load(fh)["rotation_deg"]
ENERGY = float(particle_rows(a.particles, core=True)[:, 8].mean())
a1, a2 = a.quads
stamp("DOUBLET ROTATION SCAN of {}: q1 {:+.0f} / q2 {:+.0f} V, basis rotation [{:g}, {:g}] deg on top of R {:+.3f} deg; "
      "alpha1 {} x alpha2 {} at {} particles, res {:.1f} mm".format(os.path.basename(RUN), a.q[0], a.q[1], a1, a2, R,
                                                                     a.alpha1, a.alpha2, a.n, 1e3 * a.res))

# ---------------------------------------------------------------- 1. the skew bases (the normal ones exist)
common = ["--step-dir", STEPS, "--voltages", os.path.join(GEO, "voltages.csv"), "--state", os.path.join(GEO, "state.pickle"),
          "--out-dir", BASIS, "--res", a.res, "--bfield", a.bfield, "--rotate-all", R, "--energy-mev", ENERGY]
for tag in ("spiral", "q1", "q2"):
    if not os.path.exists(os.path.join(BASIS, "ef_itp_{}.pickle".format(tag))):
        raise SystemExit("missing basis {} in {}".format(tag, BASIS))
for tag, extra in (("q1skew", ["--spiral-voltage", 0, "--quad-voltages", 3500, 0, "--rotate-quads", a1 + 45.0, a2, "--no-test"]),
                   ("q2skew", ["--spiral-voltage", 0, "--quad-voltages", 0, 3500, "--rotate-quads", a1, a2 + 45.0, "--no-test"])):
    if not os.path.exists(os.path.join(BASIS, "ef_itp_{}.pickle".format(tag))):
        t0 = time.time()
        run(["-m", "spyral_inflector.tracking.bem_reload", "--tag", tag] + extra + common, "log_basis_{}.txt".format(tag))
        stamp("  basis {} ({:.0f} s)".format(tag, time.time() - t0))
    else:
        stamp("  basis {} exists".format(tag))

# ---------------------------------------------------------------- 2. the rotation grid at the fixed voltages
scan_json = os.path.join(BASIS, "scan2_{}.json".format(a.tag))
if not os.path.exists(scan_json):
    stamp("scanning {} x {} rotations = {} points".format(len(a.alpha1), len(a.alpha2), len(a.alpha1) * len(a.alpha2)))
    run([os.path.join(SCRIPTS, "BunchScan2.py"), "--reload-dir", BASIS, "--out-dir", BASIS, "--step-dir", STEPS, "--bfield", a.bfield,
         "--basis", "spiral", "q1", "q1skew", "q2", "q2skew", "--phi", 90.0 + R, "--particles", a.particles, "--n", a.n,
         "--points", a.q[0], a.q[1], "--alpha1"] + list(a.alpha1) + ["--alpha2"] + list(a.alpha2) +
        ["--tag", a.tag, "--rank", "vert", "--rank-tol", a.rank_tol], "log_scan_{}.txt".format(a.tag))
with open(scan_json) as fh:
    rows = [r for r in json.load(fh)["results"] if r.get("vfom") is not None]
t_max = max(r["transmission"] for r in rows)
ranked = sorted([r for r in rows if r["transmission"] >= t_max - a.rank_tol], key=lambda r: (r["vfom"], -r["transmission"]))
best = ranked[0]
ref = next((r for r in rows if r["alpha1"] == 0.0 and r["alpha2"] == 0.0), None)
by_t = max(rows, key=lambda r: r["transmission"])

# ---------------------------------------------------------------- 3. report
L = ["# Doublet quad-rotation scan of {}".format(os.path.basename(RUN)), "",
     "Quads fixed at q1 {:+.0f} / q2 {:+.0f} V, basis rotation [{:g}, {:g}] deg on top of R = {:+.3f} deg; extra rotations "
     "alpha1 x alpha2 by superposition of the normal and skew basis fields ({:.1f} mm grid), {} particles per point, "
     "vacuum. Ranking: rms vertical angle vz/v_long at 55 mm past the exit within {:.0f} points of the best transmission. "
     "The collision geometry of the quads stays at the basis orientation.".format(a.q[0], a.q[1], a1, a2, R, 1e3 * a.res, a.n,
                                                                                 100 * a.rank_tol), ""]
if ref is not None:
    L.append("Reference (0, 0): {:.1f} mrad rms at {:.1f} %.".format(ref["vfom"], 100 * ref["transmission"]))
L.append("Best transmission: {:.1f} % at ({:+.0f}, {:+.0f}) deg, {:.1f} mrad.".format(100 * by_t["transmission"], by_t["alpha1"], by_t["alpha2"], by_t["vfom"]))
L.append("**Vertically ranked best: {:.1f} mrad rms at {:.1f} %, alpha1 {:+.0f} / alpha2 {:+.0f} deg.**".format(
    best["vfom"], 100 * best["transmission"], best["alpha1"], best["alpha2"]))
L += ["", "| alpha1 [deg] | alpha2 [deg] | transmission | z' rms [mrad] | z rms at 55 mm [mm] | spiral loss |", "|---|---|---|---|---|---|"]
for r in sorted(rows, key=lambda r: (r["alpha1"], r["alpha2"])):
    mark = " **" if r is best else ""
    L.append("| {:+.0f} | {:+.0f} | {:.1f} % | {:.1f}{} | {:.2f} | {:.1f} % |".format(
        r["alpha1"], r["alpha2"], 100 * r["transmission"], r["vfom"], mark, r.get("z_asym_rms_mm", float("nan")), 100 * r["lost_spiral"]))
with open(os.path.join(RUN, "rot_scan.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
stamp("done: best {:.1f} mrad at {:.1f} % (alpha {:+.0f}/{:+.0f}); reference {}; report -> {}".format(
    best["vfom"], 100 * best["transmission"], best["alpha1"], best["alpha2"],
    "-" if ref is None else "{:.1f} mrad at {:.1f} %".format(ref["vfom"], 100 * ref["transmission"]), os.path.join(RUN, "rot_scan.md")))
