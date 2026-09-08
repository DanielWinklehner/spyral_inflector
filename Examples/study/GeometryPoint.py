"""One geometry point of the wiggle study.

Builds the deck inflector with knob overrides, lets the trajectory optimizer set the
spiral voltage and dz for the design particle (truncations held at the deck values),
exports the STEP files, solves three basis fields from them (spiral only, quad 1 at
3500 V, quad 2 at 3500 V, with the quads rotated as requested), re-tunes the two quad
voltages on a small grid by superposition with the rotated Bevatech beam, tracks a
larger bunch at the best point and records transmission, losses by category and the
exit-beam metrics. Everything goes to Results/wiggle/<name>/.

    python GeometryPoint.py --name base --phi 90 --rotate-quads 0 8 --q1 3750 6750 3 --q2 -6750 -3750 3
    python GeometryPoint.py --name sigma_p --sigma 0.0032 ...
"""
import argparse
import csv
import datetime
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exit_metrics import exit_metrics, fmt as fmt_metrics  # noqa: E402
from spyral_inflector import *  # noqa: F401,F403,E402

DECK = os.environ.get("SI_DECK") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SCRIPTS = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

parser = argparse.ArgumentParser()
parser.add_argument("--name", required=True)
parser.add_argument("--out-root", default=os.path.join(DECK, "Results", "wiggle"))
parser.add_argument("--steps-root", default=os.path.join(DECK, "Geometry", "wiggle"), help="folder for <name>_steps")
parser.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit_core_as_txt.txt"))
parser.add_argument("--bfield", default=os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle"))
parser.add_argument("--phi", type=float, default=0.0, help="beam rotation [deg]")
parser.add_argument("--swap-xy", action="store_true", help="exchange x and y of the input beam (before the rotation)")
parser.add_argument("--rotate-quads", type=float, nargs=2, default=[0.0, 0.0], metavar=("A1", "A2"))
parser.add_argument("--q1", type=float, nargs=3, default=[3750, 6750, 3], metavar=("MIN", "MAX", "N"))
parser.add_argument("--q2", type=float, nargs=3, default=[-6750, -3750, 3], metavar=("MIN", "MAX", "N"))
parser.add_argument("--n-scan", type=int, default=1500)
parser.add_argument("--n-bunch", type=int, default=5000)
parser.add_argument("--res", type=float, default=0.005, help="potential grid resolution [m], optimizer and basis fields")
parser.add_argument("--maxiter", type=int, default=8)
parser.add_argument("--fix-truncations", type=float, nargs=2, default=[0.34, 0.77], metavar=("ENT", "EXIT"))
parser.add_argument("--fix-dz", type=float, default=None, help="hold the axial shift at this value [m]; only the voltage is optimized")
parser.add_argument("--energy-mev", type=float, default=None,
                    help="design-particle kinetic energy [MeV] (default: mean of the particle file)")
# geometry knobs (deck values as defaults)
parser.add_argument("--volt", type=float, default=12000.0)
parser.add_argument("--gap", type=float, default=0.019)
parser.add_argument("--tilt", type=float, default=31.0, help="k' tilt [deg] (kept fixed in the study)")
parser.add_argument("--dx", type=float, default=0.01)
parser.add_argument("--sigma", type=float, default=0.0022, help="V-depth [m]")
parser.add_argument("--aspect", type=float, default=2.4, help="electrode width / gap")
parser.add_argument("--gamma", type=float, default=5.0, help="gammaAng, exit tilt [deg]")
parser.add_argument("--angling", type=float, default=11.0, help="anglingAng, inner face angling [deg]")
parser.add_argument("--quad-bore", type=float, default=0.013, help="quadrupole hyperbola vertex radius a = b [m]")
parser.add_argument("--quad-z1", type=float, default=-0.27, help="quad 1 start z [m]")
parser.add_argument("--quad-z2", type=float, default=-0.19, help="quad 2 start z [m]")
parser.add_argument("--quad-len", type=float, default=0.045, help="quad 1 length [m] (and quad 2 unless --quad-len2)")
parser.add_argument("--quad-len2", type=float, default=None, help="quad 2 length [m]")
parser.add_argument("--shared-plates", action="store_true",
                    help="one grounded plate between the quads; quad 2 then starts at z1 + len1 + 2 gap + plate thickness (overrides --quad-z2)")
parser.add_argument("--aper-hole", type=float, default=0.0125, help="quad aperture hole radius [m] (deck: aper_rad 0.025 = 12.5 mm hole)")
parser.add_argument("--entrance-hole", type=float, default=None,
                    help="hole radius of the first plate (facing the RFQ exit) [m]; default = --aper-hole")
parser.add_argument("--plate-gap", type=float, default=0.001, help="axial gap between a quad's grounded aperture plates and its pole ends [m] (deck: 1 mm)")
parser.add_argument("--plate-thickness", type=float, default=0.005, help="thickness of the quad aperture plates [m]")
parser.add_argument("--slot-width", type=float, default=0.015, help="inflector entrance aperture slot width (gap direction) [m]")
parser.add_argument("--slot-length", type=float, default=0.040, help="inflector entrance aperture slot length [m]")
parser.add_argument("--exit-opening", type=float, nargs=2, default=None, metavar=("ACROSS", "ALONG"),
                    help="housing exit opening [m]: across the tilted gap direction (default = slot width) and along it (default = slot length)")
args = parser.parse_args()
if args.quad_len2 is None:
    args.quad_len2 = args.quad_len
if args.shared_plates:
    args.quad_z2 = args.quad_z1 + args.quad_len + 2.0 * args.plate_gap + args.plate_thickness

out_dir = os.path.join(args.out_root, args.name)
os.makedirs(out_dir, exist_ok=True)
steps_dir = os.path.join(args.steps_root, "{}_steps".format(args.name))
os.makedirs(args.steps_root, exist_ok=True)
t_start = time.time()
knobs = {k: getattr(args, k) for k in ("volt", "gap", "tilt", "dx", "sigma", "aspect", "gamma", "angling",
                                        "quad_bore", "quad_z1", "quad_z2", "quad_len", "aper_hole", "slot_width", "slot_length", "exit_opening",
                                        "plate_gap", "plate_thickness", "quad_len2", "shared_plates", "entrance_hole")}
print("=" * 78)
print("GEOMETRY POINT '{}': {}".format(args.name, knobs))
print("  beam rotation {:.1f} deg, quad rotations {} deg".format(args.phi, args.rotate_quads))
print("=" * 78, flush=True)
summary = {"name": args.name, "knobs": knobs, "phi": args.phi, "swap_xy": bool(args.swap_xy), "rotate_quads": args.rotate_quads,
           "fixed_truncations_deg": args.fix_truncations, "res_m": args.res}


def log(step, msg):
    print("[{:>5.0f} s] {}: {}".format(time.time() - t_start, step, msg), flush=True)


def run(cmd, logname):
    """Run a script as a subprocess, log to a file, fail loudly."""
    with open(os.path.join(out_dir, logname), "w", encoding="utf-8") as fh:
        p = subprocess.run([PY, "-u"] + cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=out_dir)
    if p.returncode != 0:
        raise RuntimeError("{} failed (exit {}), see {}".format(cmd[0], p.returncode, logname))


# ---------------------------------------------------------------- build + optimize
H2P = ParticleDistribution(species=IonSpecies("H2_1+"))
if args.energy_mev is None:
    _e = np.loadtxt(args.particles, skiprows=1, usecols=8)
    args.energy_mev = float(np.mean(_e))
H2P.set_mean_energy_z_mev(args.energy_mev)
summary["energy_mev"] = args.energy_mev
print("  design energy {:.5f} MeV".format(args.energy_mev), flush=True)
si = SpiralInflector(ion=H2P, method="numerical", solver="bempp",
                     volt=args.volt, gap=args.gap, tilt=args.tilt, dx=args.dx, sigma=args.sigma,
                     vee_shape="parabolic", ns=100, aspect_ratio=args.aspect, rotation=0.0,
                     debug=False, gammaAng=args.gamma, anglingAng=args.angling)
si.load_bfield(bfield=args.bfield)
si.initialize()
si.set_parameter(key="h", value=0.005)
si.set_parameter(key="make_aperture", value=True)
si.set_parameter(key="aperture_params", value={"thickness": 4e-3, "radius": 50e-3, "length": args.slot_length, "width": args.slot_width,
                                               "top_distance": 5e-3, "bottom_distance": 10e-3,
                                               "hole_type": "rectangle", "voltage": 0.0,
                                               **({"exit_width": args.exit_opening[0], "exit_length": args.exit_opening[1]}
                                                  if args.exit_opening else {})})
si.set_parameter(key="make_housing", value=True)
si.set_parameter(key="housing_params", value={"zmin": -0.12, "zmax": 0.03, "span": True, "gap": 6e-3,
                                              "thickness": 4e-3, "voltage": 0.0, "experimental": True})
si.set_parameter(key="make_quadrupoles", value=True)
si.set_parameter(key="quadrupole_params", value={"a": args.quad_bore, "b": args.quad_bore, "radius": 0.04,
                                                 "z_starts": [args.quad_z1, args.quad_z2],
                                                 "lengths": [args.quad_len, args.quad_len2],
                                                 "voltages": [3500, 3500], "aper_rad": 2.0 * args.aper_hole,
                                                 "plate_gap": args.plate_gap, "plate_thickness": args.plate_thickness,
                                                 "shared_plates": bool(args.shared_plates),
                                                 "entrance_aper_rad": 2.0 * (args.entrance_hole if args.entrance_hole else args.aper_hole)})
si.generate_geometry()
log("geometry", "generated")

result = si.optimize_trajectory(maxiter=args.maxiter, solver="dfols", res=args.res,
                                initial_guess=[args.fix_truncations[0], args.fix_truncations[1],
                                               1.85e-3 if args.fix_dz is None else args.fix_dz, 0.97],
                                fixed=({0: args.fix_truncations[0], 1: args.fix_truncations[1]} if args.fix_dz is None else
                                       {0: args.fix_truncations[0], 1: args.fix_truncations[1], 2: args.fix_dz}),
                                bounds=((0.0, 15.0), (0.0, 15.0), (-15.0e-3, 15.0e-3), (0.85, 1.3)),
                                exclude_quadrupoles=True)
m = result["measurements"]
summary["optimizer"] = {
    "status": result["status"], "converged": result["converged"], "n_evaluations": result["n_evaluations"],
    "dz_mm": 1e3 * result["dz"], "volt_scale": result["volt_scale"], "voltage_V": result["voltage"],
    "residual_final": {"angle_deg": float(result["residual_final"][0]), "centering_mm": 1e3 * float(result["residual_final"][1]),
                       "z_offset_mm": 1e3 * float(result["residual_final"][2]), "width_mm": 1e3 * float(result["residual_final"][3])},
    "clearance_exit_mm": {"anode": 1e3 * m["clearance_anode_exit"], "cathode": 1e3 * m["clearance_cathode_exit"]},
    "min_clearance_mm": 1e3 * m["min_clearance"], "exit_point_mm": (1e3 * m["exit_point"]).tolist()}
log("optimizer", "voltage {:.1f} V, dz {:+.2f} mm, final residuals angle {:+.3f} deg z {:+.2f} mm ({})".format(
    result["voltage"], 1e3 * result["dz"], result["residual_final"][0], 1e3 * result["residual_final"][2], result["status"]))
with open(os.path.join(out_dir, "summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)

# ---------------------------------------------------------------- export (as OptimizeTrajectory.py)
def sanitize(name):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(name)).strip("_")


assembly = si.numerical_variables["objects"]
names, used = {}, set()
for i, electrode in enumerate(assembly.electrodes.values()):
    name = sanitize(electrode.name) or "electrode_{:03d}".format(i)
    unique, k = name, 2
    while unique.lower() in used:
        unique = "{}_{}".format(name, k)
        k += 1
    used.add(unique.lower())
    names[i] = unique
electrode_voltages = {names[i]: float(e.voltage) for i, e in enumerate(assembly.electrodes.values())}
mesh = si.numerical_variables["full mesh"]
domain_names = {int(e.bempp_domain): names[i] for i, e in enumerate(assembly.electrodes.values())}
ctx = result["context"]
state = {"trj_design": si.analytic_variables["trj_design"] + result["shift"],
         "v_design": si.analytic_variables["v_design"],
         "trj_design_full": si.analytic_variables["trj_design_full"] + result["shift"],
         "v_design_full": si.analytic_variables["v_design_full"],
         "b_lim_deg": result["b_lim_deg"], "shift": result["shift"], "shift_lab": result["shift_lab"],
         "volt_scale": result["volt_scale"], "voltage": result["voltage"], "electrode_voltages": electrode_voltages,
         "z_start": ctx["z_start"], "v0": ctx["v0"], "dt": ctx["dt"], "nsteps": ctx["nsteps"],
         "track_r": m["r"], "track_v": m["v"],
         "mesh": {"verts": np.asarray(mesh["verts"]), "elems": np.asarray(mesh["elems"]),
                  "domns": np.asarray(mesh["domns"]), "names": domain_names},
         "residual_final": np.asarray(result["residual_final"]), "res_m": args.res, "knobs": knobs}
state_fn = os.path.join(out_dir, "state.pickle")
with open(state_fn, "wb") as fh:
    pickle.dump(state, fh)
volt_fn = os.path.join(out_dir, "voltages.csv")
with open(volt_fn, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Electrode", "Voltage (V)"])
    for name, volt in electrode_voltages.items():
        w.writerow([name, "{:.3f}".format(volt)])
if os.path.isdir(steps_dir):
    shutil.rmtree(steps_dir)
os.makedirs(steps_dir)
for i, electrode in enumerate(assembly.electrodes.values()):
    if electrode.export(os.path.join(steps_dir, "{:03d}_{}.step".format(i, names[i]))) != 0:
        raise RuntimeError("STEP export failed for {}".format(electrode.name))
log("export", "{} STEP files -> {}".format(len(names), steps_dir))
del si  # free the BEM solution / GPU memory before the subprocesses start

# ---------------------------------------------------------------- basis fields from the STEP files
common = ["--step-dir", steps_dir, "--voltages", volt_fn, "--state", state_fn, "--out-dir", out_dir,
          "--res", str(args.res), "--bfield", args.bfield, "--rotate-quads", str(args.rotate_quads[0]), str(args.rotate_quads[1])]
run([os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "spiral", "--quad-voltages", "0", "0"] + common, "log_reload_spiral.txt")
rl = json.load(open(os.path.join(out_dir, "reload_spiral.json")))
summary["reload_test_particle"] = rl.get("test_particle")
log("basis", "spiral field: test particle exit point differs by {:.3f} mm from the optimizer track".format(
    rl["test_particle"]["exit_point_diff_mm"]))
run([os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "q1", "--spiral-voltage", "0", "--quad-voltages", "3500", "0", "--no-test"] + common, "log_reload_q1.txt")
run([os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "q2", "--spiral-voltage", "0", "--quad-voltages", "0", "3500", "--no-test"] + common, "log_reload_q2.txt")
log("basis", "three basis fields solved")
with open(os.path.join(out_dir, "summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)

# ---------------------------------------------------------------- inner quad scan by superposition
run([os.path.join(SCRIPTS, "BunchScan2.py"), "--reload-dir", out_dir, "--out-dir", out_dir, "--step-dir", steps_dir,
     "--bfield", args.bfield, "--basis", "spiral", "q1", "q1", "q2", "q2", "--phi", str(args.phi)] + (["--swap-xy"] if args.swap_xy else []) + [
     "--q1", str(args.q1[0]), str(args.q1[1]), str(int(args.q1[2])), "--q2", str(args.q2[0]), str(args.q2[1]), str(int(args.q2[2])),
     "--particles", args.particles, "--n", str(args.n_scan), "--tag", "inner"], "log_scan_inner.txt")
sc = json.load(open(os.path.join(out_dir, "scan2_inner.json")))
summary["inner_scan"] = {"best": sc["best"], "points": [{k: r[k] for k in ("q1", "q2", "transmission", "lost_spiral", "lost_quads", "lost_apertures")} for r in sc["results"]]}
log("inner scan", "best q1 {:.0f} q2 {:.0f}: {:.1f} % ({} particles)".format(sc["best"]["q1"], sc["best"]["q2"], 100 * sc["best"]["transmission"], args.n_scan))

# ---------------------------------------------------------------- bunch at the best point
run([os.path.join(SCRIPTS, "BunchTrack.py"), "--tag", "inner_best", "--reload-dir", out_dir, "--step-dir", steps_dir,
     "--particles", args.particles, "--bfield", args.bfield, "--phi", str(args.phi), "--n", str(args.n_bunch)] + (["--swap-xy"] if args.swap_xy else []) + [
     "--out-tag", "bunch"], "log_bunch.txt")
metrics = exit_metrics(os.path.join(out_dir, "bunch_bunch.npz"), os.path.join(out_dir, "si_state_inner_best.pickle"))
summary["bunch"] = metrics
summary["wall_s"] = time.time() - t_start
with open(os.path.join(out_dir, "summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)
log("done", fmt_metrics(metrics))
print("wrote {}".format(os.path.join(out_dir, "summary.json")), flush=True)
