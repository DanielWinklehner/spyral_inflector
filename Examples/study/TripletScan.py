"""Quadrupole TRIPLET in the injection line, tuned for a small rms vertical angle at the exit.

Builds the inflector of a finished final-push run (its knobs and rotated field; shape
overrides optional) with THREE electrostatic quadrupoles in the same 140 mm between the RFQ
exit plane and the entrance aperture plate (shared grounded plates; the last quad ends one
plate gap before the entrance plate), optimizes the design orbit (exit truncation and dz
free), solves the rotation R from the exit-plate rule, solves seven basis fields of the
rigidly rotated system (spiral, q1, q1 skew, q2, q2 skew, q3, q3 skew -- the skew ones let
the third quad's rotation be scanned by superposition) and scans the three voltages (and
optionally extra rotations of quad 3) on a bunch, ranked by the rms vertical angle vz/v_long
at 55 mm past the exit within --rank-tol of the best transmission.

    python TripletScan.py --name triplet1
    python TripletScan.py --name triplet1 --sigma 0.0004 --angling 15 --gamma 5 --lens 0.035 0.040 0.035

Run from this folder with PYTHONPATH set to the triplet worktree and SI_DECK to the deck
(RunTriplet.ps1 does that). Outputs: <deck>/Results/triplet/<name>/ (geometry, steps, basis,
scan2_*.json, report.md).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
PY = sys.executable
sys.path.insert(0, SCRIPTS)

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--name", required=True)
p.add_argument("--base-run", default=os.path.join(DECK, "Results", "final", "final1"))
p.add_argument("--sigma", type=float, default=None, help="V-depth [m] (default: the base run's)")
p.add_argument("--angling", type=float, default=None)
p.add_argument("--gamma", type=float, default=None)
p.add_argument("--lens", type=float, nargs=3, default=[0.035, 0.040, 0.035], help="quad lengths [m]; the last is re-derived to end one plate gap before the entrance plate")
p.add_argument("--quads", type=float, nargs=3, default=[16.0, 24.0, 24.0], help="quad rotations on top of R [deg] (basis orientation)")
p.add_argument("--alpha3", type=float, nargs="+", default=[0.0, 20.0, 40.0], help="extra rotations of quad 3 scanned by superposition [deg]")
p.add_argument("--q1", type=float, nargs=3, default=[4000, 10000, 4], metavar=("MIN", "MAX", "N"))
p.add_argument("--q2", type=float, nargs=3, default=[-11000, -5000, 4], metavar=("MIN", "MAX", "N"))
p.add_argument("--q3", type=float, nargs=3, default=[2000, 8000, 4], metavar=("MIN", "MAX", "N"))
p.add_argument("--n-scan", type=int, default=1000)
p.add_argument("--rank-tol", type=float, default=0.05)
p.add_argument("--fix-entrance", type=float, default=0.34)
p.add_argument("--maxiter", type=int, default=4)
p.add_argument("--geo-res", type=float, default=0.005)
p.add_argument("--res", type=float, default=0.005)
p.add_argument("--h", type=float, default=0.005)
p.add_argument("--gap-azimuth", type=float, default=34.0)
p.add_argument("--half-gap", type=float, default=0.005)
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle"))
p.add_argument("--tag", default="scan")
args = p.parse_args()

OUT = os.path.join(DECK, "Results", "triplet", args.name)
GEO, STEPS, BASIS = (os.path.join(OUT, d) for d in ("geometry", "steps", "basis"))
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")
T0 = time.time()


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def run(cmd, log_path):
    with open(log_path, "w", encoding="utf-8") as fh:
        r = subprocess.run([PY, "-u"] + [str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS)
    if r.returncode != 0:
        raise SystemExit("{} failed (exit {}), see {}".format(os.path.basename(str(cmd[0])), r.returncode, log_path))


def jload(fn):
    with open(fn) as fh:
        return json.load(fh)


from spyral_inflector.tracking.deck import particle_rows  # noqa: E402
from spyral_inflector.tracking.geometry_point import build_geometry  # noqa: E402
from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402

knobs = dict(jload(os.path.join(args.base_run, "geometry", "summary.json"))["knobs"])
for k in ("sigma", "angling", "gamma"):
    if getattr(args, k) is not None:
        knobs[k] = getattr(args, k)
knobs.update(n_quads=3, quad_len=args.lens[0], quad_len2=args.lens[1], quad_len3=args.lens[2], shared_plates=True, q2_exit_plate=False)
ROT_BF = os.path.join(args.base_run, "bfield_rotated.pickle")
ENERGY = float(particle_rows(args.particles, core=True)[:, 8].mean())
TRUNC = (args.fix_entrance, None)
stamp("TRIPLET '{}': shape sigma {:.1f} mm / angling {:g} / gamma {:g}, quad lengths {} mm, rotations {} deg; base {}".format(
    args.name, 1e3 * knobs["sigma"], knobs["angling"], knobs["gamma"], [1e3 * x for x in args.lens], args.quads, os.path.basename(args.base_run)))

# ---------------------------------------------------------------- 1. geometry + design orbit
if not os.path.exists(os.path.join(GEO, "summary.json")):
    stamp("building the triplet geometry and the design orbit (exit truncation + dz free)")
    if os.path.isdir(GEO):
        shutil.rmtree(GEO)
    build_geometry(GEO, STEPS, ROT_BF, ENERGY, knobs=knobs, fix_truncations=TRUNC, maxiter=args.maxiter, res=args.geo_res, h=args.h, log=stamp)
opt = jload(os.path.join(GEO, "summary.json"))["optimizer"]
names = sorted(f for f in os.listdir(STEPS) if f.endswith(".step"))
stamp("  electrodes: {}".format(", ".join(n.split("_", 1)[1][:-5] for n in names)))
spec = exit_plane_spec(STEPS, os.path.join(GEO, "state.pickle"), args.gap_azimuth, args.half_gap,
                       out_json=os.path.join(OUT, "exit_plane.json"), log=stamp)
R = spec["rotation_deg"]
stamp("  R = {:+.3f} deg; V {:.0f} V, dz {:+.2f} mm, truncations {}, residuals angle {:+.3f} deg, z {:+.2f} mm ({})".format(
    R, opt["voltage_V"], opt["dz_mm"], opt.get("truncations_deg"), opt["residual_final"]["angle_deg"], opt["residual_final"]["z_offset_mm"], opt["status"][:40]))

# ---------------------------------------------------------------- 2. seven basis fields
os.makedirs(BASIS, exist_ok=True)
a1, a2, a3 = args.quads
common = ["--step-dir", STEPS, "--voltages", os.path.join(GEO, "voltages.csv"), "--state", os.path.join(GEO, "state.pickle"),
          "--out-dir", BASIS, "--res", args.res, "--h", args.h, "--bfield", args.bfield, "--rotate-all", R, "--energy-mev", ENERGY]
bases = (("spiral", ["--quad-voltages", 0, 0, 0, "--rotate-quads", a1, a2, a3]),
         ("q1", ["--spiral-voltage", 0, "--quad-voltages", 3500, 0, 0, "--rotate-quads", a1, a2, a3, "--no-test"]),
         ("q1skew", ["--spiral-voltage", 0, "--quad-voltages", 3500, 0, 0, "--rotate-quads", a1 + 45, a2, a3, "--no-test"]),
         ("q2", ["--spiral-voltage", 0, "--quad-voltages", 0, 3500, 0, "--rotate-quads", a1, a2, a3, "--no-test"]),
         ("q2skew", ["--spiral-voltage", 0, "--quad-voltages", 0, 3500, 0, "--rotate-quads", a1, a2 + 45, a3, "--no-test"]),
         ("q3", ["--spiral-voltage", 0, "--quad-voltages", 0, 0, 3500, "--rotate-quads", a1, a2, a3, "--no-test"]),
         ("q3skew", ["--spiral-voltage", 0, "--quad-voltages", 0, 0, 3500, "--rotate-quads", a1, a2, a3 + 45, "--no-test"]))
for tag, extra in bases:
    if not os.path.exists(os.path.join(BASIS, "ef_itp_{}.pickle".format(tag))):
        t0 = time.time()
        run(["-m", "spyral_inflector.tracking.bem_reload", "--tag", tag] + extra + common, os.path.join(OUT, "log_basis_{}.txt".format(tag)))
        stamp("  basis {} ({:.0f} s)".format(tag, time.time() - t0))
rl = jload(os.path.join(BASIS, "reload_spiral.json"))
qc = rl.get("quad_field_check") or {}
stamp("  quad field check: " + "; ".join("{} z {:+.3f} m".format(k, v["z_m"]) for k, v in qc.items()))

# ---------------------------------------------------------------- 3. the voltage (and quad-3 rotation) grid
scan_json = os.path.join(BASIS, "scan2_{}.json".format(args.tag))
if not os.path.exists(scan_json):
    n_pts = int(args.q1[2]) * int(args.q2[2]) * int(args.q3[2]) * len(args.alpha3)
    stamp("scanning {} x {} x {} voltages x {} quad-3 rotations = {} points at {} particles (vertical ranking, tol {:.0%})".format(
        int(args.q1[2]), int(args.q2[2]), int(args.q3[2]), len(args.alpha3), n_pts, args.n_scan, args.rank_tol))
    run([os.path.join(SCRIPTS, "BunchScan2.py"), "--reload-dir", BASIS, "--out-dir", BASIS, "--step-dir", STEPS, "--bfield", args.bfield,
         "--basis", "spiral", "q1", "q1skew", "q2", "q2skew", "q3", "q3skew", "--phi", 90.0 + R, "--particles", args.particles,
         "--n", args.n_scan, "--q1", args.q1[0], args.q1[1], int(args.q1[2]), "--q2", args.q2[0], args.q2[1], int(args.q2[2]),
         "--q3", args.q3[0], args.q3[1], int(args.q3[2]), "--alpha3"] + list(args.alpha3) +
        ["--tag", args.tag, "--rank", "vert", "--rank-tol", args.rank_tol], os.path.join(OUT, "log_scan_{}.txt".format(args.tag)))
sc = jload(scan_json)
rows = [r for r in sc["results"] if r.get("vfom") is not None]
t_max = max(r["transmission"] for r in rows)
cands = sorted([r for r in rows if r["transmission"] >= t_max - args.rank_tol], key=lambda r: (r["vfom"], -r["transmission"]))
best = cands[0]
by_t = max(rows, key=lambda r: r["transmission"])

# ---------------------------------------------------------------- report
L = ["# Triplet scan {}".format(args.name), "",
     "Shape sigma {:.1f} mm, angling {:g}, gamma {:g}; quads {} mm long at rotations {} deg on top of R = {:+.3f} deg; design orbit "
     "V {:.0f} V, dz {:+.2f} mm, exit truncation {}. Grid {}x{}x{} voltages x quad-3 rotations {} at {} particles, {:.0f} mm fields; "
     "ranking: rms vertical angle vz/v_long at 55 mm past the exit within {:.0f} points of the best transmission.".format(
         1e3 * knobs["sigma"], knobs["angling"], knobs["gamma"], [1e3 * x for x in args.lens], args.quads, R, opt["voltage_V"], opt["dz_mm"],
         "{:.2f} deg".format(opt["truncations_deg"][1]) if opt.get("truncations_deg") else "-", int(args.q1[2]), int(args.q2[2]), int(args.q3[2]),
         args.alpha3, args.n_scan, 1e3 * args.res, 100 * args.rank_tol), "",
     "Best transmission on the grid: {:.1f} % at q1 {:+.0f} / q2 {:+.0f} / q3 {:+.0f} V (+{:g} deg on quad 3), {:.1f} mrad.".format(
         100 * by_t["transmission"], by_t["q1"], by_t["q2"], by_t["q3"], by_t.get("alpha3", 0), by_t["vfom"]), "",
     "**Vertically ranked best: {:.1f} mrad rms at {:.1f} %**, q1 {:+.0f} / q2 {:+.0f} / q3 {:+.0f} V (+{:g} deg on quad 3); "
     "z rms {:.1f} mm at 55 mm, envelope {:.1f} mm.".format(100 * 0 + best["vfom"], 100 * best["transmission"], best["q1"], best["q2"], best["q3"],
                                                            best.get("alpha3", 0), best["z_asym_rms_mm"], best.get("z_env_mean_mm") or float("nan")), "",
     "| rank | q1 | q2 | q3 | +alpha3 | transmission | z' rms [mrad] | z rms [mm] | spiral loss |", "|---|---|---|---|---|---|---|---|---|"]
for i, r in enumerate(cands[:15]):
    L.append("| {} | {:+.0f} | {:+.0f} | {:+.0f} | {:g} | {:.1f} % | {:.1f} | {:.1f} | {:.1f} % |".format(
        i + 1, r["q1"], r["q2"], r["q3"], r.get("alpha3", 0), 100 * r["transmission"], r["vfom"], r["z_asym_rms_mm"], 100 * r["lost_spiral"]))
L += ["", "Doublet reference: take it from Results/vfocus/vfocus1/report.md (5 mm screen) or a final-push run's own quad scan "
      "at the SAME basis-field resolution and particle count as this scan. Do NOT compare a 2.5 mm triplet number with a "
      "5 mm doublet number: halving the basis grid is worth about 13 mrad and 9 transmission points on this geometry "
      "(Results/triplet/triplet2wide vs triplet2wide_push, identical everything else).", "",
      "Wall time {:.0f} min.".format((time.time() - T0) / 60)]
with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
with open(os.path.join(OUT, "summary.json"), "w") as fh:
    json.dump({"args": vars(args), "knobs": knobs, "R_deg": R, "optimizer": opt, "best": best, "best_transmission": by_t,
               "n_points": len(rows)}, fh, indent=2)
stamp("done: best {:.1f} mrad at {:.1f} % (q {:+.0f}/{:+.0f}/{:+.0f}, +{:g} deg); report -> {}".format(
    best["vfom"], 100 * best["transmission"], best["q1"], best["q2"], best["q3"], best.get("alpha3", 0), os.path.join(OUT, "report.md")))
