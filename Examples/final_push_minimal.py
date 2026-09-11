"""The final-push pipeline end to end on synthetic inputs: a rehearsal that runs anywhere.

This is Scripts/RunFinalPush.py of the HCHC-60 deck without the deck: the B-field map and
the RFQ particle file are generated here (a cyclotron-like field in the Baseline frame, a
Gaussian bunch with a slow, late tail), then the same eleven phases run at rehearsal
settings, writing the same folders, files and report as the real run.

    python final_push_minimal.py                       # -> ./final_push_minimal_out/
    python final_push_minimal.py --out /tmp/fp --phases 6-11    # resume from the basis solves

Phases (each writes phase_N.done and is skipped on a re-run):
  1  field intake: frame detection (a Baseline map, bore at +z, is mirrored in memory), checks
  2  geometry + design orbit at R = 0 (build_geometry: knobs -> optimizer -> STEP + state)
  3  the housing's exit plate: rotation R0 that puts it half a gap before the 34 deg gap
  4  the field a system rotated by R0 sees (rotate_bfield_z)
  5  re-optimize in that field, re-solve R (fixed point on R)
  6  three basis solves of the rigidly rotated system (bem_reload --rotate-all)
  7  quad retune on a small voltage grid by superposing the basis fields (BunchScan2)
  8  vacuum run: core bunch + the tail injected at its arrival time (bunch --stragglers)
  9  the same with space charge
 10  Baseline-frame exports in mm: STEP per electrode + assembly, E-field, hand-offs, exit plane
 11  plots (assembly rotated to match the tracks) and report.md

About an hour on a GPU machine at the default settings (two design-orbit solves of 7 min,
basis solves 12 min, quad grid 10 min, the two bunch runs 12 min each). Needs spyral_inflector, PyPATools,
py_electrodes, bempp and dfols; the quad retune calls Examples/study/BunchScan2.py.

Frames, in one sentence: the deck integrates the mirror image of the machine (beam from
-z), the Baseline (machine) frame has +z up with the beam coming down; field maps in either
frame are accepted and everything exported goes out in the Baseline frame, in millimetres.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.join(HERE, "study")
PY = sys.executable

# the HCHC-60 final-push knobs (Results/final/final1); aperture radius 47 mm is the library default
KNOBS = dict(volt=12000.0, gap=0.019, tilt=31.0, dx=0.01, sigma=0.0004, aspect=2.4, gamma=11.0, angling=7.0,
             quad_bore=0.018, quad_z1=-0.264, quad_z2=-0.194, quad_len=0.055, quad_len2=0.06, shared_plates=True,
             aper_hole=0.0175, entrance_hole=None, plate_gap=0.005, plate_thickness=0.005,
             slot_width=0.019, slot_length=0.04, exit_opening=[0.023, 0.04],
             housing_gap=0.003, housing_thickness=0.002, top_distance=0.005, bottom_distance=0.01, rotation=0.0,
             q2_exit_plate=False)
QUADS = (16.0, 24.0)          # the quads' own rotation on top of R
BASE_PHI = 90.0               # beam rotation of the unrotated system
RF_MHZ = 32.8
E0_MEV = 0.0684               # H2+ from the RFQ, 34.2 keV/u

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--out", default=os.path.join(HERE, "final_push_minimal_out"))
p.add_argument("--phases", default="1-11", help="e.g. 6-11")
p.add_argument("--n-core", type=int, default=3000, help="particles in the synthetic bunch core")
p.add_argument("--n-tail", type=int, default=300, help="slow, late particles (the RFQ's unaccelerated tail)")
p.add_argument("--maxiter", type=int, default=1, help="design-orbit optimizer iterations per build (each ~100 s)")
p.add_argument("--res", type=float, default=0.005, help="field grid resolution [m] (optimizer and basis fields)")
p.add_argument("--h", type=float, default=0.005, help="surface mesh size [m]")
p.add_argument("--fix-truncations", type=float, nargs=2, default=[0.34, 0.77], help="entrance/exit truncation [deg]")
p.add_argument("--gap-azimuth", type=float, default=34.0)
p.add_argument("--half-gap", type=float, default=0.005, help="[m] along the exit face normal")
p.add_argument("--r-tol", type=float, default=1.0, help="stop re-optimizing when R moves by less than this [deg]")
p.add_argument("--current-ma", type=float, default=8.0)
p.add_argument("--flutter", type=float, default=0.08, help="4-fold azimuthal field variation at r = 70 mm (fraction)")
p.add_argument("--seed", type=int, default=1)
args = p.parse_args()
lo, hi = (args.phases.split("-") + [args.phases])[:2]
PHASES = set(range(int(lo), int(hi) + 1))

OUT = os.path.abspath(args.out)
INP = os.path.join(OUT, "inputs")
os.makedirs(INP, exist_ok=True)
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


def run(cmd, log_name, cwd=STUDY):
    stamp("  > {}".format(" ".join(os.path.basename(str(c)) if str(c).endswith(".py") else str(c) for c in cmd)))
    with open(os.path.join(OUT, log_name), "w", encoding="utf-8") as fh:
        r = subprocess.run([PY, "-u"] + [str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT, cwd=cwd)
    if r.returncode != 0:
        raise SystemExit("{} failed (exit {}), see {}".format(cmd[0], r.returncode, os.path.join(OUT, log_name)))


def jload(fn):
    with open(fn) as fh:
        return json.load(fh)


from PyPATools.field import Field  # noqa: E402
from spyral_inflector.tracking.deck import particle_rows  # noqa: E402
from spyral_inflector.tracking.frames import bfield_frame, to_deck_frame, rotate_bfield_z  # noqa: E402
from spyral_inflector.tracking.geometry_point import build_geometry  # noqa: E402
from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402
from spyral_inflector.tracking.export import steps_to_baseline_mm, field_to_machine_frame  # noqa: E402
from spyral_inflector.tracking.plots import main as plot_main  # noqa: E402

# ---------------------------------------------------------------- synthetic inputs
BFIELD = os.path.join(INP, "bfield_synthetic_baseline.pickle")
PARTICLES = os.path.join(INP, "rfq_synthetic.txt")


def make_bfield(path, flutter, step=0.005):
    """A cyclotron-like B-field on a grid in the BASELINE frame: bore at positive z, Bz < 0
    in the median plane (positive ions circulate counter-clockwise seen from +z).

    Axisymmetric part: Bz = B0 g(z) with g = 1 / (1 + (z / z0)^2), which reproduces the
    HCHC-60 axis (-0.91 T at z = 0, -0.71 T at 50 mm); Br = -(r/2) B0 g'(z) and the r^2 term
    of Bz make it Maxwell-consistent to second order in r. On top a 3-fold flutter
    B = grad psi with psi = eps B0 (r^3 cos 3phi / r0^3) z g(z)^2 -- eps of the main field
    at r0 in the median plane, decaying with the main field away from it (the z g^2 factor
    keeps it bounded, at the price of a small curl). It is there so that rotating the
    inflector in the field, phases 3 to 5, changes what it sees, as in the real magnet;
    3-fold rather than 4-fold because the rotation comes out near -90 deg, a period of a
    4-fold pattern, which would make the re-optimization a no-op."""
    B0, z0, r0 = -0.91, 0.094, 0.07
    x = np.arange(-0.12, 0.12 + 1e-9, step)
    y = np.arange(-0.12, 0.12 + 1e-9, step)
    z = np.arange(-0.06, 0.40 + 1e-9, step)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    R2 = X ** 2 + Y ** 2
    u = Z / z0
    g = 1.0 / (1.0 + u ** 2)
    dg = -2.0 * u / z0 * g ** 2
    d2g = (6.0 * u ** 2 - 2.0) / z0 ** 2 * g ** 3
    bz = B0 * (g - 0.25 * R2 * d2g)
    br_over_r = -0.5 * B0 * dg                      # Br / r
    bx, by = br_over_r * X, br_over_r * Y
    # flutter: psi = k q(x, y) w(z), q = r^3 cos 3phi = x^3 - 3 x y^2, w = z g^2
    k = flutter * B0 / r0 ** 3
    q = X ** 3 - 3.0 * X * Y ** 2
    w = Z * g ** 2
    dw = g ** 2 + 2.0 * Z * g * dg
    bx += k * w * (3.0 * X ** 2 - 3.0 * Y ** 2)
    by += k * w * (-6.0 * X * Y)
    bz += k * q * dw
    field = Field.from_arrays(grid={"x": x, "y": y, "z": z}, values={"x": bx, "y": by, "z": bz},
                              label="synthetic cyclotron field, Baseline frame", dim=3, units="m",
                              interpolator_backend="scipy")
    field.save(path)
    return field


def make_particles(path, n_core, n_tail, seed):
    """A synthetic RFQ output in the 10-column text format the package reads
    (x mm, x' mrad, y mm, y' mrad, z mm, z' mrad, phase deg, time s, energy MeV, loss):
    a Gaussian core at E0 and a tail of slow particles arriving 200..720 deg late. The
    package's core selection (|phase| <= 180 deg, energy >= half the median) separates them."""
    rng = np.random.default_rng(seed)
    ion_mass_mev = 1877.268
    gamma = 1.0 + E0_MEV / ion_mass_mev
    beta = np.sqrt(1.0 - 1.0 / gamma ** 2)
    beta_lambda = beta * 299792458.0 / (RF_MHZ * 1e6)
    phase_core = rng.normal(0.0, 20.0, n_core)
    core = np.column_stack([rng.normal(0.0, 1.5, n_core), rng.normal(0.0, 10.0, n_core),
                            rng.normal(0.0, 1.5, n_core), rng.normal(0.0, 10.0, n_core),
                            -phase_core / 360.0 * beta_lambda * 1e3, np.zeros(n_core),
                            phase_core, np.zeros(n_core), rng.normal(E0_MEV, 0.004 * E0_MEV, n_core), np.zeros(n_core)])
    tail = np.column_stack([rng.normal(0.0, 2.5, n_tail), rng.normal(0.0, 20.0, n_tail),
                            rng.normal(0.0, 2.5, n_tail), rng.normal(0.0, 20.0, n_tail),
                            np.zeros(n_tail), np.zeros(n_tail),
                            rng.uniform(200.0, 720.0, n_tail), np.zeros(n_tail),
                            rng.uniform(0.10, 0.45, n_tail) * E0_MEV, np.zeros(n_tail)])
    rows = np.vstack([core, tail])
    np.savetxt(path, rows, header="x(mm) x'(mrad) y(mm) y'(mrad) z(mm) z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss",
               fmt="%.6e", comments="")
    return rows


stamp("FINAL PUSH (minimal, synthetic inputs) -> {}".format(OUT))
if not os.path.exists(BFIELD):
    stamp("making the synthetic B-field ({:.0%} flutter at 70 mm)".format(args.flutter))
    make_bfield(BFIELD, args.flutter)
if not os.path.exists(PARTICLES):
    stamp("making the synthetic RFQ file ({} core + {} tail particles)".format(args.n_core, args.n_tail))
    make_particles(PARTICLES, args.n_core, args.n_tail, args.seed)
ENERGY = float(particle_rows(PARTICLES, core=True)[:, 8].mean())
stamp("  phases {}, design energy {:.5f} MeV (core mean)".format(sorted(PHASES), ENERGY))

# ---------------------------------------------------------------- 1. field intake
if 1 in PHASES and not done(1):
    f = Field.from_file(BFIELD)
    z = np.asarray(f.grid["z"])
    frame = bfield_frame(f)
    bz0 = float(f(np.array([[0.0, 0.0, 0.0]]))[0, 2])
    stamp("phase 1: map z {:+.3f}..{:+.3f} m, frame '{}', Bz(0,0,0) = {:+.4f} T".format(z.min(), z.max(), frame, bz0))
    if bz0 >= 0:
        raise SystemExit("Bz(0,0,0) must be negative for counter-clockwise circulation of positive ions")
    fd = to_deck_frame(f, log=stamp)          # Baseline map -> mirrored in memory; a deck map is left alone
    prof = {"{:+.2f}".format(zz): float(fd(np.array([[0.0, 0.0, zz]]))[0, 2]) for zz in (-0.30, -0.20, -0.10, -0.05, 0.0)}
    stamp("  deck-frame on-axis Bz: " + ", ".join("{} m: {:+.4f} T".format(k, v) for k, v in prof.items()))
    mark(1, {"pickle": BFIELD, "frame": frame, "bz0": bz0, "on_axis_deck": prof})

# ---------------------------------------------------------------- 2-5. geometry, exit plate, rotation, iterate
GEO0, STEPS0 = os.path.join(OUT, "geometry0"), os.path.join(OUT, "steps0")
if 2 in PHASES and not done(2):
    stamp("phase 2: geometry + design orbit at R = 0")
    build_geometry(GEO0, STEPS0, BFIELD, ENERGY, knobs=KNOBS, fix_truncations=tuple(args.fix_truncations),
                   maxiter=args.maxiter, res=args.res, h=args.h, log=stamp)
    mark(2, jload(os.path.join(GEO0, "summary.json"))["optimizer"])
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
    rotate_bfield_z(to_deck_frame(Field.from_file(BFIELD), log=stamp), R, log=stamp).save(ROT_BF)
    mark(4, {"R": R, "file": ROT_BF})
if 5 in PHASES and not done(5):
    history = [R]
    for it in range(2):
        stamp("phase 5.{}: re-optimize in the field rotated by {:+.3f} deg".format(it + 1, R))
        if os.path.isdir(GEO):
            shutil.rmtree(GEO)
        build_geometry(GEO, STEPS, ROT_BF, ENERGY, knobs=KNOBS, fix_truncations=tuple(args.fix_truncations),
                       maxiter=args.maxiter, res=args.res, h=args.h, log=stamp)
        spec = exit_plane_spec(STEPS, os.path.join(GEO, "state.pickle"), args.gap_azimuth, args.half_gap,
                               out_json=os.path.join(OUT, "exit_plane.json"), log=stamp)
        R_new = spec["rotation_deg"]
        stamp("  rotation {:+.3f} -> {:+.3f} deg (change {:+.3f})".format(R, R_new, R_new - R))
        history.append(R_new)
        if abs(R_new - R) <= args.r_tol:
            R = R_new
            break
        R = R_new
        rotate_bfield_z(to_deck_frame(Field.from_file(BFIELD), log=stamp), R, log=stamp).save(ROT_BF)
    mark(5, {"R_history": history, "R": R, "optimizer": jload(os.path.join(GEO, "summary.json"))["optimizer"]})
if os.path.exists(os.path.join(OUT, "exit_plane.json")):
    R = jload(os.path.join(OUT, "exit_plane.json"))["rotation_deg"]
PHI = BASE_PHI + R
stamp("  system rotation R = {:+.3f} deg, beam angle {:+.3f} deg".format(R, PHI))

# ---------------------------------------------------------------- 6. basis fields of the rotated system
BASIS = os.path.join(OUT, "basis")
os.makedirs(BASIS, exist_ok=True)
common = ["--step-dir", STEPS, "--voltages", os.path.join(GEO, "voltages.csv"), "--state", os.path.join(GEO, "state.pickle"),
          "--out-dir", BASIS, "--res", args.res, "--h", args.h, "--bfield", BFIELD, "--rotate-all", R,
          "--rotate-quads", QUADS[0], QUADS[1], "--energy-mev", ENERGY]
if 6 in PHASES and not done(6):
    stamp("phase 6: basis solves at {:.1f} mm, R = {:+.3f} deg".format(1e3 * args.res, R))
    for tag, extra in (("spiral", ["--quad-voltages", 0, 0]),
                       ("q1", ["--spiral-voltage", 0, "--quad-voltages", 3500, 0, "--no-test"]),
                       ("q2", ["--spiral-voltage", 0, "--quad-voltages", 0, 3500, "--no-test"])):
        if not os.path.exists(os.path.join(BASIS, "ef_itp_{}.pickle".format(tag))):
            run(["-m", "spyral_inflector.tracking.bem_reload", "--tag", tag] + extra + common, "log_basis_{}.txt".format(tag))
    mark(6, jload(os.path.join(BASIS, "reload_spiral.json")).get("test_particle"))

# ---------------------------------------------------------------- 7. quad retune (3 x 3, superposition)
if 7 in PHASES and not done(7):
    stamp("phase 7: quad retune, 3x3 grid by superposition")
    if not os.path.exists(os.path.join(BASIS, "scan2_coarse.json")):
        run([os.path.join(STUDY, "BunchScan2.py"), "--reload-dir", BASIS, "--out-dir", BASIS, "--step-dir", STEPS,
             "--bfield", BFIELD, "--basis", "spiral", "q1", "q1", "q2", "q2", "--phi", PHI, "--particles", PARTICLES,
             "--n", 300, "--q1", 5075, 7075, 3, "--q2", -8450, -6450, 3, "--tag", "coarse"], "log_scan_coarse.txt")
    best = jload(os.path.join(BASIS, "scan2_coarse.json"))["best"]
    stamp("  best: q1 {:+.0f}, q2 {:+.0f}, {:.1f} %".format(best["q1"], best["q2"], 100 * best["transmission"]))
    for fn in ("scan2_{}.json", "ef_itp_{}_best.pickle", "si_state_{}_best.pickle"):
        shutil.copyfile(os.path.join(BASIS, fn.format("coarse")), os.path.join(BASIS, fn.format("fine")))
    mark(7, best)

# ---------------------------------------------------------------- 8/9. every particle: core bunch + injected tail
bunch_common = ["-m", "spyral_inflector.tracking.bunch", "--tag", "fine_best", "--reload-dir", BASIS, "--step-dir", STEPS,
                "--particles", PARTICLES, "--bfield", BFIELD, "--phi", PHI, "--n", 10 ** 7, "--stragglers", "--rf-mhz", RF_MHZ,
                "--handoff-frame", "machine", "--handoff-distance", 0.0, "--current-ma", args.current_ma]
if 8 in PHASES and not done(8):
    stamp("phase 8: vacuum, every particle (core bunch + tail injected by phase)")
    run(bunch_common + ["--out-tag", "vac_all", "--save-openpmd", os.path.join(BASIS, "handoff_vac_all.h5")], "log_bunch_vac_all.txt")
    mark(8, jload(os.path.join(BASIS, "bunch_vac_all.json")).get("tail"))
if 9 in PHASES and not done(9):
    stamp("phase 9: space charge {:g} mA, every particle".format(args.current_ma))
    run(bunch_common + ["--sc", "--h", 0.004, "--resolve-every", 16, "--sc-min-particles", 100,
                        "--reference", os.path.join(BASIS, "bunch_vac_all.json"),
                        "--out-tag", "sc_all", "--save-openpmd", os.path.join(BASIS, "handoff_sc_all.h5")], "log_bunch_sc_all.txt")
    mark(9, jload(os.path.join(BASIS, "bunch_sc_all.json")).get("tail"))

# ---------------------------------------------------------------- 10. exports, Baseline frame, mm
EXP = os.path.join(OUT, "export_baseline")
if 10 in PHASES and not done(10):
    stamp("phase 10: exports in the Baseline frame (mm)")
    os.makedirs(EXP, exist_ok=True)
    files = steps_to_baseline_mm(STEPS, os.path.join(EXP, "steps_baseline_mm"), rotation_deg=R, quad_rotation=QUADS,
                                 combined="inflector_minimal_baseline_mm.step", log=stamp)
    field_to_machine_frame(os.path.join(BASIS, "ef_itp_fine_best.pickle"), os.path.join(EXP, "efield_fine_best_baseline.pickle"))
    for fn in ("handoff_vac_all.h5", "handoff_sc_all.h5"):
        if os.path.exists(os.path.join(BASIS, fn)):
            shutil.copyfile(os.path.join(BASIS, fn), os.path.join(EXP, fn))
    shutil.copyfile(os.path.join(OUT, "exit_plane.json"), os.path.join(EXP, "exit_plane_spec.json"))
    with open(os.path.join(EXP, "README.txt"), "w", encoding="utf-8") as fh:
        fh.write("Baseline (machine) frame: right-handed, +z up, the beam enters from +z; azimuth counter-clockwise from above.\n"
                 "System rotation R = {:+.4f} deg about z; quads {:+.1f} / {:+.1f} deg on top.\n"
                 "steps_baseline_mm/  one STEP per electrode + the combined assembly, MILLIMETRES\n"
                 "efield_fine_best_baseline.pickle  E-field of the tuned point (spiral + quads), V/m on a metre grid\n"
                 "handoff_*.h5  openPMD hand-off at the electrode exit (vacuum / space charge)\n"
                 "exit_plane_spec.json  the exit plate's outer face at R\n".format(R, QUADS[0], QUADS[1]))
    mark(10, {"files": files})

# ---------------------------------------------------------------- 11. plots and report
if 11 in PHASES and not done(11):
    stamp("phase 11: plots and report")
    for tag in ("vac", "sc"):
        npz = os.path.join(BASIS, "bunch_{}_all.npz".format(tag))
        if os.path.exists(npz):
            plot_main([npz, "--step-dir", STEPS, "--out", os.path.join(OUT, "geometry_trajectories_{}_all.png".format(tag)),
                       "--side", os.path.join(OUT, "geometry_side_{}_all.png".format(tag)),
                       "--title", "minimal final push, {}: every particle (core + tail)".format("vacuum" if tag == "vac" else "space charge")])
    spec = jload(os.path.join(OUT, "exit_plane.json"))
    geo = jload(os.path.join(GEO, "summary.json"))["optimizer"]
    best = jload(os.path.join(BASIS, "scan2_fine.json"))["best"]
    vac = jload(os.path.join(BASIS, "bunch_vac_all.json"))
    sc = jload(os.path.join(BASIS, "bunch_sc_all.json")) if os.path.exists(os.path.join(BASIS, "bunch_sc_all.json")) else None
    L = ["# Minimal final push (synthetic field and bunch)", "",
         "System rotation R = **{:+.3f} deg** (history {}).".format(R, ", ".join("{:+.3f}".format(x) for x in jload(os.path.join(OUT, "phase_5.done"))["info"]["R_history"])), "",
         "Design orbit: voltage {:.1f} V, dz {:+.2f} mm; residuals angle {:+.3f} deg, z {:+.2f} mm ({})".format(
             geo["voltage_V"], geo["dz_mm"], geo["residual_final"]["angle_deg"], geo["residual_final"]["z_offset_mm"], geo["status"]), "",
         "Exit plate (Baseline frame): r {:.2f} mm, azimuth {:+.3f} deg, z {:+.2f} mm; the point {:.1f} mm along the normal at azimuth {:+.3f} deg (gap at {:g} deg).".format(
             spec["baseline"]["r_mm"], spec["baseline"]["azimuth_deg"], spec["baseline"]["point_mm"][2],
             spec["half_gap_mm"], spec["pushed_point_azimuth_deg"], spec["gap_azimuth_deg"]), "",
         "Quads: q1 {:+.0f} V, q2 {:+.0f} V ({:.1f} % on the scan bunch)".format(best["q1"], best["q2"], 100 * best["transmission"]), "",
         "| run | particles | transmitted | core | tail |", "|---|---|---|---|---|"]
    for name, b in (("vacuum", vac), ("space charge {:g} mA".format(args.current_ma), sc)):
        if b is None:
            continue
        t = b.get("tail") or {}
        L.append("| {} | {:,d} | {:.2f} % | {} | {} |".format(
            name, b["n_particles"], 100 * b["transmission"],
            "{:.2f} % of {:,d}".format(100 * t["core_transmission"], t["n_core"]) if t else "-",
            "{:.2f} % of {:,d}".format(100 * t["tail_transmission"], t["n_tail"]) if t else "-"))
    L += ["", "![iso](geometry_trajectories_sc_all.png)", "", "![side](geometry_side_sc_all.png)", "",
          "Exports: `export_baseline/`. Total wall time {:.0f} min.".format((time.time() - T_START) / 60)]
    with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    mark(11)
    stamp("report -> {}".format(os.path.join(OUT, "report.md")))
stamp("done ({:.0f} min)".format((time.time() - T_START) / 60))
