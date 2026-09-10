"""The full spiral-inflector re-calculation for a new field map, one command, resumable.

    python RunFinalPush.py --name final1 --comsol "Fields/<export>.comsol"
    python RunFinalPush.py --name rehearsal --bfield Fields/HCHC-60_CentralBField_z-40to5cm_1mm.pickle --quick

Phases (each writes <out>/phase_N.done and is skipped on a re-run):
  1  convert the COMSOL text export (Baseline frame, as-is) to a pickle, fix the seams, validate
  2  build the geometry and optimize the design orbit (voltage, dz) at R = 0
  3  measure the exit plate's outer face, solve the rotation R0 (5 mm along the normal to the 34 deg gap)
  4  rotate the field by R0 (the field a system rotated by R0 sees)
  5  re-optimize the geometry in that field, re-solve R (iterate while |dR| > 0.2 deg)
  6  three basis solves of the rigidly rotated system at R (res 2.5 mm)
  7  quad retune: coarse grid, then a fine grid around the best
  8  vacuum run, every particle of the .dst (core as a bunch, tail injected at its arrival time)
  9  the same with space charge at 8 mA -- the final run
 10  exports in the Baseline frame, mm: per-electrode STEP + assembly, E-field, hand-offs, exit-plane spec
 11  report.md
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

# the housing3 geometry (Results/housing/housing3/geometry/summary.json), aperture radius 47 mm is the library default
KNOBS = dict(volt=12000.0, gap=0.019, tilt=31.0, dx=0.01, sigma=0.0004, aspect=2.4, gamma=11.0, angling=7.0,
             quad_bore=0.018, quad_z1=-0.264, quad_z2=-0.194, quad_len=0.055, quad_len2=0.06, shared_plates=True,
             aper_hole=0.0175, entrance_hole=None, plate_gap=0.005, plate_thickness=0.005,
             slot_width=0.019, slot_length=0.04, exit_opening=[0.023, 0.04],
             housing_gap=0.003, housing_thickness=0.002, top_distance=0.005, bottom_distance=0.01, rotation=0.0)
QUADS = (16.0, 24.0)          # the quads' own rotation on top of R
BASE_PHI = 90.0               # beam rotation of the unrotated system

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--name", required=True)
p.add_argument("--comsol", default=None, help="COMSOL text export in the Baseline frame (phase 1)")
p.add_argument("--bfield", default=None, help="an existing field pickle instead of phase 1")
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--gap-azimuth", type=float, default=34.0)
p.add_argument("--half-gap", type=float, default=0.005, help="[m] along the exit face normal")
p.add_argument("--no-q2-exit-plate", action="store_true", help="quad 2 terminated by the inflector entrance plate")
p.add_argument("--quad-len2", type=float, default=None, help="override quad 2 length [m] (e.g. 0.0705 with --no-q2-exit-plate)")
p.add_argument("--fix-truncations", type=float, nargs=2, default=[0.34, 0.77])
p.add_argument("--res", type=float, default=0.0025, help="basis-field resolution [m]")
p.add_argument("--geo-res", type=float, default=0.005, help="resolution of the design-orbit optimizer [m]")
p.add_argument("--maxiter", type=int, default=8)
p.add_argument("--q1", type=float, nargs=3, default=[4075, 8075, 5], metavar=("MIN", "MAX", "N"))
p.add_argument("--q2", type=float, nargs=3, default=[-9450, -5450, 5], metavar=("MIN", "MAX", "N"))
p.add_argument("--fine-step", type=float, default=500.0, help="fine grid: 3x3 at this step around the coarse best")
p.add_argument("--n-scan", type=int, default=1500)
p.add_argument("--n-final", type=int, default=10 ** 7, help="particles of the final runs (>= file size = all)")
p.add_argument("--current-ma", type=float, default=8.0)
p.add_argument("--sc-min-particles", type=int, default=200)
p.add_argument("--handoff-distance", type=float, default=0.0, help="[m] of design path past the electrode exit")
p.add_argument("--r-tol", type=float, default=0.2, help="re-optimize until the rotation moves by less than this [deg]")
p.add_argument("--crop-xy", type=float, default=0.15)
p.add_argument("--quick", action="store_true", help="rehearsal settings: coarse grids, few particles, 5 mm fields")
p.add_argument("--phases", default="1-11", help="e.g. 6-11")
args = p.parse_args()

if args.quick:
    args.res, args.geo_res, args.maxiter, args.n_scan = 0.005, 0.005, 4, 500
    args.q1, args.q2, args.fine_step = [5075, 7075, 3], [-8450, -6450, 3], 0.0
    args.n_final = 3000
lo, hi = (args.phases.split("-") + [args.phases])[:2]
PHASES = set(range(int(lo), int(hi) + 1))

OUT = os.path.join(DECK, "Results", "final", args.name)
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log_run.txt"), "a", encoding="utf-8")
T_START = time.time()


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def done(n):
    return os.path.exists(os.path.join(OUT, "phase_{}.done".format(n)))


def mark(n, info=None):
    with open(os.path.join(OUT, "phase_{}.done".format(n)), "w") as fh:
        json.dump({"time": time.strftime("%Y-%m-%d %H:%M:%S"), "elapsed_s": time.time() - T_START, "info": info}, fh, indent=2)


def run(cmd, log_name, cwd=SCRIPTS):
    stamp("  > {}".format(" ".join(os.path.basename(c) if os.sep in str(c) and str(c).endswith(".py") else str(c) for c in cmd)))
    with open(os.path.join(OUT, log_name), "w", encoding="utf-8") as fh:
        r = subprocess.run([PY, "-u"] + [str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT, cwd=cwd)
    if r.returncode != 0:
        raise SystemExit("{} failed (exit {}), see {}".format(cmd[0], r.returncode, os.path.join(OUT, log_name)))


def jload(fn):
    with open(fn) as fh:
        return json.load(fh)


def state_of(build):
    return jload(os.path.join(build, "summary.json"))


knobs = dict(KNOBS)
if args.no_q2_exit_plate:
    knobs["q2_exit_plate"] = False
if args.quad_len2 is not None:
    knobs["quad_len2"] = args.quad_len2
BFIELD = args.bfield or os.path.join(DECK, "Fields", "{}_baseline_1mm.pickle".format(args.name))
stamp("FINAL PUSH '{}' -> {}".format(args.name, OUT))
stamp("  phases {}, field {}, particles {}{}".format(sorted(PHASES), BFIELD, os.path.basename(args.particles),
                                                    " [QUICK]" if args.quick else ""))

from spyral_inflector.tracking.deck import particle_rows, load_bfield  # noqa: E402
from spyral_inflector.tracking.frames import to_deck_frame, rotate_bfield_z, bfield_frame  # noqa: E402
from spyral_inflector.tracking.geometry_point import build_geometry  # noqa: E402
from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402
from spyral_inflector.tracking.export import steps_to_baseline_mm, field_to_machine_frame  # noqa: E402
from PyPATools.field import Field  # noqa: E402

ENERGY = float(particle_rows(args.particles, core=True)[:, 8].mean())
stamp("  design energy {:.6f} MeV (core mean)".format(ENERGY))

# ---------------------------------------------------------------- 1. field
if 1 in PHASES and not done(1):
    if args.comsol is None:
        if not os.path.exists(BFIELD):
            raise SystemExit("phase 1 needs --comsol (or --bfield to skip it)")
        stamp("phase 1: using the existing pickle {}".format(BFIELD))
    else:
        stamp("phase 1: converting {} ({:.1f} GB)".format(os.path.basename(args.comsol), os.path.getsize(args.comsol) / 1e9))
        run([os.path.join(SCRIPTS, "convert_bfield.py"), args.comsol, BFIELD, "--mode", "none", "--fix-seams",
             "--crop-xy", args.crop_xy], "log_convert.txt")
    f = Field.from_file(BFIELD)
    z = np.asarray(f.grid["z"])
    frame = bfield_frame(f)
    bz0 = float(f(np.array([[0.0, 0.0, 0.0]]))[0, 2])
    stamp("  map z {:+.3f}..{:+.3f} m, {} x {} x {}, frame {}, Bz(0,0,0) = {:+.4f} T".format(
        z.min(), z.max(), len(f.grid["x"]), len(f.grid["y"]), len(z), frame, bz0))
    if bz0 >= 0:
        raise SystemExit("Bz(0,0,0) is not negative: positive ions would not circulate counter-clockwise; check the coil sense")
    fd = to_deck_frame(f, log=stamp)
    prof = {}
    for zz in (-0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.0):
        prof["{:+.2f}".format(zz)] = float(fd(np.array([[0.0, 0.0, zz]]))[0, 2])
    stamp("  deck-frame on-axis Bz: " + ", ".join("{} m: {:+.4f} T".format(k, v) for k, v in prof.items()))
    # cell-to-cell jumps: iron inside the map shows as tesla-size steps (the HCHC-60 map has
    # 2.5 T steps at r ~ 7 cm, |z| ~ 8 cm); the bore the beam uses has to be mT-level
    B = np.stack([np.asarray(fd.grid_values[k], float) for k in "xyz"], -1)
    Xg, Yg, Zg = np.meshgrid(*(np.asarray(fd.grid[k], float) for k in "xyz"), indexing="ij")
    Rg = np.hypot(Xg, Yg)
    del Xg, Yg
    jumps = {}
    for name, m in (("bore r<5cm, z<2cm", (Rg < 0.05) & (Zg < 0.02)), ("r<8cm, z<2cm", (Rg < 0.08) & (Zg < 0.02)),
                    ("exit r<12cm, z -1..+2cm", (Rg < 0.12) & (Zg >= -0.01) & (Zg < 0.02))):
        worst = 0.0
        for ax in range(3):
            d = np.linalg.norm(np.diff(B, axis=ax), axis=-1)
            sl = [slice(None)] * 3
            sl[ax] = slice(None, -1)
            dd = d[m[tuple(sl)]]
            if dd.size:
                worst = max(worst, float(dd.max()))
        jumps[name] = worst
    del B, Rg, Zg
    stamp("  largest cell-to-cell |dB|: " + ", ".join("{}: {:.4f} T".format(k, v) for k, v in jumps.items()))
    if jumps["bore r<5cm, z<2cm"] > 0.05:
        stamp("  WARNING: tesla-size field steps inside the bore -- inspect the map before trusting the tracking")
    old = os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle")
    cmp_ = None
    if os.path.exists(old) and os.path.abspath(old) != os.path.abspath(BFIELD):
        fo = Field.from_file(old)
        cmp_ = {k: float(fo(np.array([[0.0, 0.0, float(k)]]))[0, 2]) for k in prof}
        stamp("  previous map on axis : " + ", ".join("{} m: {:+.4f} T".format(k, v) for k, v in cmp_.items()))
    mark(1, {"pickle": BFIELD, "frame": frame, "bz0": bz0, "on_axis_deck": prof, "previous_map": cmp_, "cell_jumps_T": jumps})

# ---------------------------------------------------------------- 2/3/4/5. geometry, rotation, field rotation, iterate
GEO0, STEPS0 = os.path.join(OUT, "geometry0"), os.path.join(OUT, "steps0")
if 2 in PHASES and not done(2):
    stamp("phase 2: geometry + design orbit at R = 0 ({} knobs)".format(len(knobs)))
    build_geometry(GEO0, STEPS0, BFIELD, ENERGY, knobs=knobs, fix_truncations=tuple(args.fix_truncations),
                   maxiter=args.maxiter, res=args.geo_res, log=stamp)
    mark(2, state_of(GEO0)["optimizer"])
if 3 in PHASES and not done(3):
    stamp("phase 3: exit plate and rotation")
    spec0 = exit_plane_spec(STEPS0, os.path.join(GEO0, "state.pickle"), args.gap_azimuth, args.half_gap,
                            out_json=os.path.join(OUT, "exit_plane_R0.json"), log=stamp)
    mark(3, {"R0": spec0["rotation_deg"]})
R = jload(os.path.join(OUT, "exit_plane_R0.json"))["rotation_deg"] if os.path.exists(os.path.join(OUT, "exit_plane_R0.json")) else None
GEO, STEPS = os.path.join(OUT, "geometry"), os.path.join(OUT, "steps")
ROT_BF = os.path.join(OUT, "bfield_rotated.pickle")
if 4 in PHASES and not done(4):
    stamp("phase 4: the field a system rotated by {:+.3f} deg sees".format(R))
    fd = to_deck_frame(Field.from_file(BFIELD), log=stamp)
    rotate_bfield_z(fd, R, log=stamp).save(ROT_BF)
    mark(4, {"R": R, "file": ROT_BF})
if 5 in PHASES and not done(5):
    history = [R]
    for it in range(3):
        stamp("phase 5.{}: re-optimize the geometry in the field rotated by {:+.3f} deg".format(it + 1, R))
        if os.path.isdir(GEO):
            shutil.rmtree(GEO)
        build_geometry(GEO, STEPS, ROT_BF, ENERGY, knobs=knobs, fix_truncations=tuple(args.fix_truncations),
                       maxiter=args.maxiter, res=args.geo_res, log=stamp)
        spec = exit_plane_spec(STEPS, os.path.join(GEO, "state.pickle"), args.gap_azimuth, args.half_gap,
                               out_json=os.path.join(OUT, "exit_plane.json"), log=stamp)
        R_new = spec["rotation_deg"]
        stamp("  rotation {:+.3f} -> {:+.3f} deg (change {:+.3f})".format(R, R_new, R_new - R))
        history.append(R_new)
        if abs(R_new - R) <= args.r_tol:
            R = R_new
            break
        R = R_new
        fd = to_deck_frame(Field.from_file(BFIELD), log=stamp)
        rotate_bfield_z(fd, R, log=stamp).save(ROT_BF)
    mark(5, {"R_history": history, "R": R, "optimizer": state_of(GEO)["optimizer"]})
if os.path.exists(os.path.join(OUT, "exit_plane.json")):
    R = jload(os.path.join(OUT, "exit_plane.json"))["rotation_deg"]
PHI = BASE_PHI + R
stamp("  system rotation R = {:+.3f} deg, beam angle {:+.3f} deg".format(R, PHI))

# ---------------------------------------------------------------- 6. basis fields of the rotated system
BASIS = os.path.join(OUT, "basis")
os.makedirs(BASIS, exist_ok=True)
common = ["--step-dir", STEPS, "--voltages", os.path.join(GEO, "voltages.csv"), "--state", os.path.join(GEO, "state.pickle"),
          "--out-dir", BASIS, "--res", args.res, "--bfield", BFIELD, "--rotate-all", R, "--rotate-quads", QUADS[0], QUADS[1],
          "--energy-mev", ENERGY]
if 6 in PHASES and not done(6):
    stamp("phase 6: basis solves at {:.1f} mm, R = {:+.3f} deg".format(1e3 * args.res, R))
    for tag, extra in (("spiral", ["--quad-voltages", 0, 0]),
                       ("q1", ["--spiral-voltage", 0, "--quad-voltages", 3500, 0, "--no-test"]),
                       ("q2", ["--spiral-voltage", 0, "--quad-voltages", 0, 3500, "--no-test"])):
        if not os.path.exists(os.path.join(BASIS, "ef_itp_{}.pickle".format(tag))):
            run(["-m", "spyral_inflector.tracking.bem_reload", "--tag", tag] + extra + common, "log_basis_{}.txt".format(tag))
    rl = jload(os.path.join(BASIS, "reload_spiral.json"))
    mark(6, {"test_particle": rl.get("test_particle"), "quad_field_check": rl.get("quad_field_check")})

# ---------------------------------------------------------------- 7. quad retune
scan_common = ["--reload-dir", BASIS, "--out-dir", BASIS, "--step-dir", STEPS, "--bfield", BFIELD,
               "--basis", "spiral", "q1", "q1", "q2", "q2", "--phi", PHI, "--particles", args.particles, "--n", args.n_scan]
if 7 in PHASES and not done(7):
    stamp("phase 7: quad retune, coarse {}x{}".format(int(args.q1[2]), int(args.q2[2])))
    if not os.path.exists(os.path.join(BASIS, "scan2_coarse.json")):
        run([os.path.join(SCRIPTS, "BunchScan2.py")] + scan_common +
            ["--q1", args.q1[0], args.q1[1], int(args.q1[2]), "--q2", args.q2[0], args.q2[1], int(args.q2[2]), "--tag", "coarse"],
            "log_scan_coarse.txt")
    best = jload(os.path.join(BASIS, "scan2_coarse.json"))["best"]
    stamp("  coarse best: q1 {:+.0f}, q2 {:+.0f}, {:.1f} %".format(best["q1"], best["q2"], 100 * best["transmission"]))
    if args.fine_step > 0:
        if not os.path.exists(os.path.join(BASIS, "scan2_fine.json")):
            s = args.fine_step
            run([os.path.join(SCRIPTS, "BunchScan2.py")] + scan_common +
                ["--q1", best["q1"] - s, best["q1"] + s, 3, "--q2", best["q2"] - s, best["q2"] + s, 3, "--tag", "fine"],
                "log_scan_fine.txt")
        best = jload(os.path.join(BASIS, "scan2_fine.json"))["best"]
        stamp("  fine best  : q1 {:+.0f}, q2 {:+.0f}, {:.1f} %".format(best["q1"], best["q2"], 100 * best["transmission"]))
    else:
        for fn in ("scan2_{}.json", "ef_itp_{}_best.pickle", "si_state_{}_best.pickle"):
            shutil.copyfile(os.path.join(BASIS, fn.format("coarse")), os.path.join(BASIS, fn.format("fine")))
    mark(7, best)
BEST = jload(os.path.join(BASIS, "scan2_fine.json"))["best"] if os.path.exists(os.path.join(BASIS, "scan2_fine.json")) else None

# ---------------------------------------------------------------- 8/9. final runs, every particle
bunch_common = ["-m", "spyral_inflector.tracking.bunch", "--tag", "fine_best", "--reload-dir", BASIS, "--step-dir", STEPS,
                "--particles", args.particles, "--bfield", BFIELD, "--phi", PHI, "--n", args.n_final, "--stragglers",
                "--handoff-frame", "machine", "--handoff-distance", args.handoff_distance, "--current-ma", args.current_ma]
if 8 in PHASES and not done(8):
    stamp("phase 8: vacuum, every particle of {} (core bunch + injected tail)".format(os.path.basename(args.particles)))
    run(bunch_common + ["--out-tag", "vac_all", "--save-openpmd", os.path.join(BASIS, "handoff_vac_all.h5")], "log_bunch_vac_all.txt")
    mark(8, jload(os.path.join(BASIS, "bunch_vac_all.json")).get("tail"))
if 9 in PHASES and not done(9):
    stamp("phase 9: space charge {:g} mA, every particle".format(args.current_ma))
    run(bunch_common + ["--sc", "--sc-min-particles", args.sc_min_particles, "--reference", os.path.join(BASIS, "bunch_vac_all.json"),
                        "--out-tag", "sc_all", "--save-openpmd", os.path.join(BASIS, "handoff_sc_all.h5")], "log_bunch_sc_all.txt")
    mark(9, jload(os.path.join(BASIS, "bunch_sc_all.json")).get("tail"))

# ---------------------------------------------------------------- 10. exports
EXP = os.path.join(OUT, "export_baseline")
if 10 in PHASES and not done(10):
    stamp("phase 10: exports in the Baseline frame (mm)")
    os.makedirs(EXP, exist_ok=True)
    files = steps_to_baseline_mm(STEPS, os.path.join(EXP, "steps_baseline_mm"), rotation_deg=R, quad_rotation=QUADS,
                                 combined="HCHC60_inflector_{}_baseline_mm.step".format(args.name), log=stamp)
    field_to_machine_frame(os.path.join(BASIS, "ef_itp_fine_best.pickle"), os.path.join(EXP, "efield_fine_best_baseline.pickle"))
    for fn in ("handoff_vac_all.h5", "handoff_sc_all.h5"):
        if os.path.exists(os.path.join(BASIS, fn)):
            shutil.copyfile(os.path.join(BASIS, fn), os.path.join(EXP, fn))
    shutil.copyfile(os.path.join(OUT, "exit_plane.json"), os.path.join(EXP, "exit_plane_spec.json"))
    with open(os.path.join(EXP, "README.txt"), "w", encoding="utf-8") as fh:
        fh.write("Baseline (machine) frame: right-handed, +z up, the beam enters from +z; azimuth counter-clockwise from above.\n"
                 "System rotation R = {:+.4f} deg about z; quads {:+.1f} / {:+.1f} deg on top.\n"
                 "steps_baseline_mm/       one STEP per electrode + the combined assembly, MILLIMETRES (header and coordinates)\n"
                 "efield_fine_best_baseline.pickle   E-field of the optimized point (spiral + quads), V/m on a metre grid\n"
                 "handoff_*.h5             openPMD hand-off at {:.0f} mm past the electrode exit (vacuum / space charge)\n"
                 "exit_plane_spec.json     the exit plate's outer face at R\n"
                 "The B-field is the map given to this run ({}), already in the Baseline frame.\n".format(
                     R, QUADS[0], QUADS[1], 1e3 * args.handoff_distance, os.path.basename(BFIELD)))
    mark(10, {"files": files})

# ---------------------------------------------------------------- 11. report
if 11 in PHASES and not done(11):
    stamp("phase 11: report")
    spec = jload(os.path.join(OUT, "exit_plane.json"))
    geo = state_of(GEO)
    vac = jload(os.path.join(BASIS, "bunch_vac_all.json"))
    sc = jload(os.path.join(BASIS, "bunch_sc_all.json")) if os.path.exists(os.path.join(BASIS, "bunch_sc_all.json")) else None
    p5 = jload(os.path.join(OUT, "phase_5.done"))["info"]
    L = ["# Final push: {}".format(args.name), "",
         "Field `{}` ({}). Geometry optimized in the field rotated by R; the rigidly rotated system at R = **{:+.3f} deg** "
         "for everything downstream (history {}).".format(os.path.basename(BFIELD), "quick rehearsal" if args.quick else "full settings", R,
                                                          ", ".join("{:+.3f}".format(x) for x in p5["R_history"])), "",
         "## Design orbit", "",
         "voltage {:.1f} V, dz {:+.2f} mm, residuals: angle {:+.3f} deg, centering {:+.2f} mm, z {:+.2f} mm, width {:+.2f} mm ({}, {} evaluations)".format(
             geo["optimizer"]["voltage_V"], geo["optimizer"]["dz_mm"], geo["optimizer"]["residual_final"]["angle_deg"],
             geo["optimizer"]["residual_final"]["centering_mm"], geo["optimizer"]["residual_final"]["z_offset_mm"],
             geo["optimizer"]["residual_final"]["width_mm"], geo["optimizer"]["status"], geo["optimizer"]["n_evaluations"]), "",
         "## Exit plane (Baseline frame)", "",
         "outer face crossing at r {:.2f} mm, azimuth {:+.3f} deg, z {:+.2f} mm; normal azimuth {:+.2f} deg (vertical plate); "
         "the point {:.1f} mm along the normal sits at azimuth {:+.3f} deg (gap at {:g} deg).".format(
             spec["baseline"]["r_mm"], spec["baseline"]["azimuth_deg"], spec["baseline"]["point_mm"][2],
             spec["baseline"]["normal_azimuth_deg"], spec["half_gap_mm"], spec["pushed_point_azimuth_deg"], spec["gap_azimuth_deg"]), "",
         "## Quads", "", "q1 {:+.0f} V, q2 {:+.0f} V ({:.1f} % on the scan bunch)".format(BEST["q1"], BEST["q2"], 100 * BEST["transmission"]), "",
         "## Every particle of the RFQ file", "",
         "| run | particles | transmitted | core | tail | intercepted |", "|---|---|---|---|---|---|"]
    for name, b in (("vacuum", vac), ("space charge {:g} mA".format(args.current_ma), sc)):
        if b is None:
            continue
        t = b.get("tail") or {}
        L.append("| {} | {:,d} | {:.2f} % | {} | {} | {:,d} |".format(
            name, b["n_particles"], 100 * b["transmission"],
            "{:.2f} % of {:,d}".format(100 * t["core_transmission"], t["n_core"]) if t else "-",
            "{:.2f} % of {:,d}".format(100 * t["tail_transmission"], t["n_tail"]) if t else "-", b["n_lost"]))
    for name, b in (("vacuum", vac), ("space charge", sc)):
        if b is None or not b.get("tail"):
            continue
        L += ["", "### Where the tail terminates ({})".format(name), "", "| electrode | tail particles | % of tail | mean z of the hit [mm, deck] |", "|---|---|---|---|"]
        for el, cnt in sorted(b["tail"]["tail_losses_by_electrode"].items(), key=lambda kv: -kv[1]):
            L.append("| {} | {:,d} | {:.1f} | {:+.1f} |".format(el, cnt, 100.0 * cnt / b["tail"]["n_tail"],
                                                                 1e3 * b["tail"]["tail_loss_z_mean_by_electrode"].get(el, float("nan"))))
    L += ["", "## Losses by electrode (all particles)", "", "| electrode | vacuum | space charge |", "|---|---|---|"]
    names = sorted(set(vac["losses_by_electrode"]) | set((sc or {}).get("losses_by_electrode", {})),
                   key=lambda k: -vac["losses_by_electrode"].get(k, 0))
    for k in names:
        L.append("| {} | {:.2f} % | {} |".format(k, 100.0 * vac["losses_by_electrode"].get(k, 0) / vac["n_particles"],
                                                  "{:.2f} %".format(100.0 * sc["losses_by_electrode"].get(k, 0) / sc["n_particles"]) if sc else "-"))
    L += ["", "Exports: `{}`".format(EXP), "", "Total wall time {:.0f} min.".format((time.time() - T_START) / 60)]
    with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    mark(11)
    stamp("report -> {}".format(os.path.join(OUT, "report.md")))
stamp("done ({:.0f} min)".format((time.time() - T_START) / 60))
