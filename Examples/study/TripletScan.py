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
p.add_argument("--aspect", type=float, default=None,
               help="electrode width / gap (default: the base run's, 2.4); narrower electrodes open the anode-cathode edge gap")
p.add_argument("--lens", type=float, nargs=3, default=[0.035, 0.040, 0.035], help="quad lengths [m]; the last is re-derived to end one plate gap before the entrance plate")
p.add_argument("--quads", type=float, nargs=3, default=[16.0, 24.0, 24.0], help="quad rotations on top of R [deg] (basis orientation)")
p.add_argument("--plate-gap", type=float, default=None, help="axial gap between a quad pole end and its aperture plate [m] (default: the base run, 1 mm; Daniel 2026-09-14: 3 mm)")
p.add_argument("--plate-thickness", type=float, default=None, help="grounded aperture plate thickness [m] (default: the base run, 5 mm; Daniel 2026-09-14: 2 mm)")
p.add_argument("--quad-z1", type=float, default=None,
               help="quad 1 start z [m, deck frame, before the design-orbit shift] (default: the base run's); the plate in front of it "
                    "spans quad_z1 - plate_gap - plate_thickness .. quad_z1 - plate_gap")
p.add_argument("--entrance-hole", type=float, default=None,
               help="hole RADIUS of the first plate [m] (HCHC-1.5: the RFQ-cyclotron vacuum sleeve, 46 mm diameter, doubles as it)")
p.add_argument("--aper-hole", type=float, default=None, help="hole RADIUS of the plates between the quads [m] (default: the base run's)")
p.add_argument("--entrance-rect", type=float, nargs=3, default=None, metavar=("HALF_U", "HALF_W", "ANGLE_DEG"),
               help="RECTANGULAR first-plate hole (the vacuum-sleeve stand-in) instead of the round --entrance-hole: "
                    "half-widths [m] and the turn of its u axis from the deck x axis [deg]; to line it up with quad 1 "
                    "give quad 1's total rotation (the plate turns with R, not with the quad)")
p.add_argument("--slot-width", type=float, default=None,
               help="inflector entrance slot WIDTH (across the electrode gap) [m] (default: the base run's)")
p.add_argument("--slot-length", type=float, default=None, help="inflector entrance slot LENGTH [m] (default: the base run's)")
p.add_argument("--skew-bases", action="store_true",
               help="solve the three skew bases even when no quad-3 rotation is scanned here (a later rotation scan of all quads needs them)")
p.add_argument("--alpha3", type=float, nargs="+", default=[0.0, 20.0, 40.0], help="extra rotations of quad 3 scanned by superposition [deg]")
p.add_argument("--q1", type=float, nargs=3, default=[4000, 10000, 4], metavar=("MIN", "MAX", "N"))
p.add_argument("--q2", type=float, nargs=3, default=[-11000, -5000, 4], metavar=("MIN", "MAX", "N"))
p.add_argument("--q3", type=float, nargs=3, default=[2000, 8000, 4], metavar=("MIN", "MAX", "N"))
p.add_argument("--n-scan", type=int, default=1000)
p.add_argument("--rank-tol", type=float, default=0.05)
p.add_argument("--rank-key", choices=("vfom_h", "vfom"), default="vfom_h",
               help="what the vertical ranking minimizes: vfom_h = rms vertical angle at the HAND-OFF plane (the cyclotron's injection "
                    "cut; default since 2026-09-14), vfom = the value 55 mm past the exit (the old metric, which hides that cut)")
p.add_argument("--fix-entrance", type=float, default=0.34)
p.add_argument("--maxiter", type=int, default=4)
p.add_argument("--geo-res", type=float, default=0.005)
p.add_argument("--res", type=float, default=0.005)
p.add_argument("--h", type=float, default=0.005)
p.add_argument("--gap-azimuth", type=float, default=34.0)
p.add_argument("--rotation-rule", choices=("offset", "parallel"), default="offset",
               help="offset: the plate's outer face lands half a gap before the gap's radial plane. "
                    "parallel: the plate is turned parallel to that plane instead (it then stands off it).")
p.add_argument("--half-gap", type=float, default=0.005)
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle"))
p.add_argument("--tag", default="scan")
args = p.parse_args()

RESULTS = os.environ.get("SI_RESULTS") or os.path.join(DECK, "Results")   # 2026-09-15: new runs on a local SSD (E:), not Dropbox
OUT = os.path.join(RESULTS, "triplet", args.name)
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


from spyral_inflector.tracking.deck import particle_rows, Z_START  # noqa: E402
from spyral_inflector.tracking.geometry_point import build_geometry  # noqa: E402
from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402

knobs = dict(jload(os.path.join(args.base_run, "geometry", "summary.json"))["knobs"])
for k in ("sigma", "angling", "gamma", "aspect"):
    if getattr(args, k) is not None:
        knobs[k] = getattr(args, k)
knobs.update(n_quads=3, quad_len=args.lens[0], quad_len2=args.lens[1], quad_len3=args.lens[2], shared_plates=True, q2_exit_plate=False)
if args.plate_gap is not None:
    knobs["plate_gap"] = args.plate_gap
if args.plate_thickness is not None:
    knobs["plate_thickness"] = args.plate_thickness
if args.quad_z1 is not None:
    knobs["quad_z1"] = args.quad_z1
if args.entrance_hole is not None:
    knobs["entrance_hole"] = args.entrance_hole
if args.aper_hole is not None:
    knobs["aper_hole"] = args.aper_hole
if args.entrance_rect is not None:
    knobs["entrance_rect"] = list(args.entrance_rect)
if args.slot_width is not None:
    knobs["slot_width"] = args.slot_width
if args.slot_length is not None:
    knobs["slot_length"] = args.slot_length
ROT_BF = os.path.join(args.base_run, "bfield_rotated.pickle")
ENERGY = float(particle_rows(args.particles, core=True)[:, 8].mean())
TRUNC = (args.fix_entrance, None)
stamp("TRIPLET '{}': shape sigma {:.1f} mm / angling {:g} / gamma {:g}, quad lengths {} mm, rotations {} deg; base {}".format(
    args.name, 1e3 * knobs["sigma"], knobs["angling"], knobs["gamma"], [1e3 * x for x in args.lens], args.quads, os.path.basename(args.base_run)))
stamp("  quad 1 at z {:+.2f} mm (before the design-orbit shift), plates {:.1f} mm thick with {:.1f} mm gaps; hole radii: first plate {:.1f} mm, "
      "between the quads {:.1f} mm; pole tips at {:.1f} mm; bunch starts at z {:+.1f} mm".format(
          1e3 * knobs["quad_z1"], 1e3 * knobs["plate_thickness"], 1e3 * knobs["plate_gap"],
          1e3 * (knobs.get("entrance_hole") or knobs["aper_hole"]), 1e3 * knobs["aper_hole"], 1e3 * knobs["quad_bore"], 1e3 * Z_START))

# ---------------------------------------------------------------- 1. geometry + design orbit
if not os.path.exists(os.path.join(GEO, "summary.json")):
    stamp("building the triplet geometry and the design orbit (exit truncation + dz free)")
    if knobs.get("entrance_rect"):
        stamp("  first plate: RECTANGULAR hole, half-widths {:.1f} x {:.1f} mm, u axis at {:+.2f} deg (deck frame, before R)".format(
            1e3 * knobs["entrance_rect"][0], 1e3 * knobs["entrance_rect"][1], knobs["entrance_rect"][2]))
    if os.path.isdir(GEO):
        shutil.rmtree(GEO)
    build_geometry(GEO, STEPS, ROT_BF, ENERGY, knobs=knobs, fix_truncations=TRUNC, maxiter=args.maxiter, res=args.geo_res, h=args.h, log=stamp)
opt = jload(os.path.join(GEO, "summary.json"))["optimizer"]
from spyral_inflector.tracking.deck import electrode_cad_files, cad_name  # noqa: E402
names = electrode_cad_files(STEPS)
stamp("  electrodes: {}".format(", ".join(cad_name(n) for n in names)))
spec = exit_plane_spec(STEPS, os.path.join(GEO, "state.pickle"), args.gap_azimuth, args.half_gap,
                       out_json=os.path.join(OUT, "exit_plane.json"), log=stamp, rule=args.rotation_rule)
R = spec["rotation_deg"]
stamp("  R = {:+.3f} deg; V {:.0f} V, dz {:+.2f} mm, truncations {}, residuals angle {:+.3f} deg, z {:+.2f} mm ({})".format(
    R, opt["voltage_V"], opt["dz_mm"], opt.get("truncations_deg"), opt["residual_final"]["angle_deg"], opt["residual_final"]["z_offset_mm"], opt["status"][:40]))

# ---------------------------------------------------------------- 2. the basis fields
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
# The skew bases only exist so a quad's ROTATION can be scanned by superposition (a quad turned
# by alpha is cos(2 alpha) * normal + sin(2 alpha) * skew). With every scanned rotation at zero
# BunchScan2 never opens them -- it still wants the seven names, but skips loading the skew files
# -- so solving them is three full-assembly BEM solves, about twelve minutes, for nothing.
if all(x == 0.0 for x in args.alpha3) and not args.skew_bases:
    bases = tuple(b for b in bases if not b[0].endswith("skew"))
    stamp("  no quad rotation is scanned, so the three skew bases are not solved")
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
        ["--tag", args.tag, "--rank", ("vert_h" if args.rank_key == "vfom_h" else "vert"), "--rank-tol", args.rank_tol], os.path.join(OUT, "log_scan_{}.txt".format(args.tag)))
sc = jload(scan_json)
KEY = args.rank_key                     # "vfom_h": rms vertical angle at the hand-off plane; "vfom": at 55 mm (the old metric)
rows = [r for r in sc["results"] if r.get(KEY) is not None]
t_max = max(r["transmission"] for r in rows)
cands = sorted([r for r in rows if r["transmission"] >= t_max - args.rank_tol], key=lambda r: (r[KEY], -r["transmission"]))
best = cands[0]
by_t = max(rows, key=lambda r: r["transmission"])

# ---------------------------------------------------------------- report
L = ["# Triplet scan {}".format(args.name), "",
     "Shape sigma {:.1f} mm, angling {:g}, gamma {:g}; quads {} mm long at rotations {} deg on top of R = {:+.3f} deg; design orbit "
     "V {:.0f} V, dz {:+.2f} mm, exit truncation {}. Grid {}x{}x{} voltages x quad-3 rotations {} at {} particles, {:.0f} mm fields; "
     "ranking: rms vertical angle vz/v_long {} within {:.0f} points of the best transmission.".format(
         1e3 * knobs["sigma"], knobs["angling"], knobs["gamma"], [1e3 * x for x in args.lens], args.quads, R, opt["voltage_V"], opt["dz_mm"],
         "{:.2f} deg".format(opt["truncations_deg"][1]) if opt.get("truncations_deg") else "-", int(args.q1[2]), int(args.q2[2]), int(args.q3[2]),
         args.alpha3, args.n_scan, 1e3 * args.res,
         "AT THE HAND-OFF PLANE (the cyclotron's injection cut, cliff near 40 mrad)" if KEY == "vfom_h" else "at 55 mm past the exit", 100 * args.rank_tol), "",
     "Best transmission on the grid: {:.1f} % at q1 {:+.0f} / q2 {:+.0f} / q3 {:+.0f} V (+{:g} deg on quad 3), {:.1f} mrad at the hand-off / {:.1f} at 55 mm.".format(
         100 * by_t["transmission"], by_t["q1"], by_t["q2"], by_t["q3"], by_t.get("alpha3", 0), by_t.get("vfom_h") or float("nan"), by_t["vfom"]), "",
     "**Vertically ranked best ({}): {:.1f} mrad rms at the hand-off, {:.1f} at 55 mm, at {:.1f} %**, q1 {:+.0f} / q2 {:+.0f} / q3 {:+.0f} V (+{:g} deg on quad 3); "
     "z rms {:.1f} mm at 55 mm, envelope {:.1f} mm.".format(KEY, best.get("vfom_h") or float("nan"), best["vfom"], 100 * best["transmission"], best["q1"], best["q2"], best["q3"],
                                                            best.get("alpha3", 0), best["z_asym_rms_mm"], best.get("z_env_mean_mm") or float("nan")), "",
     "| rank | q1 | q2 | q3 | +alpha3 | transmission | z' rms hand-off [mrad] | z' rms 55 mm | z rms [mm] | spiral loss |", "|---|---|---|---|---|---|---|---|---|---|"]
for i, r in enumerate(cands[:15]):
    L.append("| {} | {:+.0f} | {:+.0f} | {:+.0f} | {:g} | {:.1f} % | {:.1f} | {:.1f} | {:.1f} | {:.1f} % |".format(
        i + 1, r["q1"], r["q2"], r["q3"], r.get("alpha3", 0), 100 * r["transmission"], r.get("vfom_h") or float("nan"), r["vfom"], r["z_asym_rms_mm"], 100 * r["lost_spiral"]))
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
stamp("done: best {:.1f} mrad ({}) at {:.1f} % (q {:+.0f}/{:+.0f}/{:+.0f}, +{:g} deg); report -> {}".format(
    best[KEY], KEY, 100 * best["transmission"], best["q1"], best["q2"], best["q3"], best.get("alpha3", 0), os.path.join(OUT, "report.md")))
