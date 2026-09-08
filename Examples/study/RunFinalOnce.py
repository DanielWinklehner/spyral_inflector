"""Run a finished inflector geometry once: direct BEM solve of the STEP files at the final
voltages, then the full core beam through it without space charge (and, with --sc, with
PyAMG space charge), exit metrics, 3-D geometry/trajectory figures and the openPMD hand-off
files for the central region.

The pg5L geometry (55/60 mm quads with a shared plate, 5 mm plate gap, 23 x 40 mm housing
exit opening, 19 mm gap, dz +3.50 mm already baked into the STEP files):

    python RunFinalOnce.py --step-dir "..\\Geometry\\wiggle\\pg5Lb_steps" --run-dir "..\\Results\\wiggle\\pg5Lb" ^
        --out "..\\Results\\jarrett_pg5L" --q1 6575 --q2 -8450 --rotate-quads 16 24 --spiral-voltage 11594.558 --sc

Inputs it needs from the deck: the STEP folder, the run folder's voltages.csv and state.pickle
(design orbit of that geometry), the B-field map (Fields/HCHC-60_CentralBField_z-40to5cm_1mm.pickle)
and the Bevatech core beam (Particles/MIT RFQ Beamdynamics/ext3_exit_core_as_txt.txt); the
scripts TrackFromStep.py, BunchTrack.py, BunchTrackSC.py, track_inflector.py, exit_metrics.py
and plot_geometry_trajectories.py next to this one; packages py_electrodes, PyPATools
(openPMD writer), bempp-cl, warp. The beam angle is 90 deg (RFQ mounting angle); the design
energy is the mean of the core file (68.45 keV).
"""
import argparse
import json
import os
import subprocess
import sys
import time

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)

parser = argparse.ArgumentParser()
parser.add_argument("--step-dir", required=True, help="folder of the exported STEP files (one per electrode)")
parser.add_argument("--run-dir", required=True, help="folder with voltages.csv and state.pickle of this geometry")
parser.add_argument("--out", required=True, help="output folder")
parser.add_argument("--tag", default="final")
parser.add_argument("--q1", type=float, required=True, help="quad 1 voltage [V] (D0,D1 = +q1, D2,D3 = -q1)")
parser.add_argument("--q2", type=float, required=True, help="quad 2 voltage [V] (D4,D5 = +q2, D6,D7 = -q2)")
parser.add_argument("--rotate-quads", type=float, nargs=2, default=[0.0, 0.0], metavar=("A1", "A2"), help="quad rotations about z [deg]")
parser.add_argument("--spiral-voltage", type=float, required=True, help="spiral electrode voltage [V] (anode +V, cathode -V)")
parser.add_argument("--phi", type=float, default=90.0, help="rotation of the input beam about z [deg]")
parser.add_argument("--n", type=int, default=43969, help="particles (the core file has 43,969)")
parser.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit_core_as_txt.txt"))
parser.add_argument("--bfield", default=os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle"))
parser.add_argument("--res", type=float, default=0.0025, help="BEM potential grid resolution [m]")
parser.add_argument("--sc", action="store_true", help="also run with PyAMG space charge")
parser.add_argument("--h", type=float, default=0.0015, help="space-charge cell size [m] (2 mm reads +1.9, 1.5 mm +1.3, 1 mm +0.8 points above converged)")
parser.add_argument("--current-ma", type=float, default=8.0)
parser.add_argument("--no-handoff", action="store_true", help="skip the openPMD hand-off files")
parser.add_argument("--n-traj", type=int, default=300, help="trajectories drawn in the figures")
args = parser.parse_args()

PY = sys.executable
os.makedirs(args.out, exist_ok=True)


def stamp(msg):
    print("[{}] {}".format(time.strftime("%H:%M:%S"), msg), flush=True)


def run(name, cmd, wait=True):
    log = os.path.join(args.out, "log_{}.txt".format(name))
    stamp("{} -> {}".format(name, os.path.basename(log)))
    fh = open(log, "w", encoding="utf-8")
    p = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS)
    if not wait:
        return p, fh
    p.wait()
    fh.close()
    if p.returncode != 0:
        sys.exit("{} failed (exit {}), see {}".format(name, p.returncode, log))
    return None


# 1. direct BEM solve of the STEP files at the final voltages (reused if present)
if not os.path.exists(os.path.join(args.out, "ef_itp_{}.pickle".format(args.tag))):
    run("direct_solve", [PY, "-u", os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", args.tag, "--step-dir", args.step_dir,
                         "--voltages", os.path.join(args.run_dir, "voltages.csv"), "--state", os.path.join(args.run_dir, "state.pickle"),
                         "--out-dir", args.out, "--res", str(args.res), "--bfield", args.bfield,
                         "--quad-voltages", str(args.q1), str(args.q2), "--rotate-quads", str(args.rotate_quads[0]), str(args.rotate_quads[1]),
                         "--spiral-voltage", str(args.spiral_voltage)])

# 2. the core beam, without and (optionally) with space charge, in parallel
jobs = []
if not os.path.exists(os.path.join(args.out, "bunch_e_nosc.json")):
    cmd = [PY, "-u", os.path.join(SCRIPTS, "BunchTrack.py"), "--tag", args.tag, "--reload-dir", args.out, "--step-dir", args.step_dir,
           "--out-tag", "e_nosc", "--n", str(args.n), "--phi", str(args.phi), "--particles", args.particles, "--bfield", args.bfield,
           "--current-ma", str(args.current_ma), "--record", str(max(args.n_traj, 1000))]
    if not args.no_handoff:
        cmd += ["--save-openpmd", os.path.join(args.out, "handoff_nosc.h5"), "--save-mode", "both"]
    jobs.append(("e_nosc",) + run("e_nosc", cmd, wait=False))
if args.sc and not os.path.exists(os.path.join(args.out, "bunch_e_sc.json")):
    cmd = [PY, "-u", os.path.join(SCRIPTS, "BunchTrackSC.py"), "--tag", args.tag, "--reload-dir", args.out, "--step-dir", args.step_dir,
           "--out-dir", args.out, "--out-tag", "e_sc", "--n", str(args.n), "--h", str(args.h), "--phi", str(args.phi),
           "--particles", args.particles, "--bfield", args.bfield, "--current-ma", str(args.current_ma), "--record", str(max(args.n_traj, 1000))]
    if not args.no_handoff:
        cmd += ["--save-openpmd", os.path.join(args.out, "handoff_sc.h5"), "--save-mode", "both"]
    jobs.append(("e_sc",) + run("e_sc", cmd, wait=False))
for name, p, fh in jobs:
    p.wait()
    fh.close()
    if p.returncode != 0:
        sys.exit("{} failed (exit {}), see log_{}.txt".format(name, p.returncode, name))
    stamp("{} finished".format(name))

# 3. exit metrics and figures
results = {}
state = os.path.join(args.out, "si_state_{}.pickle".format(args.tag))
for key in ("e_nosc", "e_sc"):
    npz = os.path.join(args.out, "bunch_{}.npz".format(key))
    if not os.path.exists(npz):
        continue
    js = os.path.join(args.out, "exit_metrics_{}.json".format(key))
    with open(os.path.join(args.out, "exit_metrics_{}.txt".format(key)), "w", encoding="utf-8") as fh:
        subprocess.call([PY, os.path.join(SCRIPTS, "exit_metrics.py"), npz, state, js], cwd=SCRIPTS, stdout=fh, stderr=subprocess.STDOUT)
    if os.path.exists(js):
        with open(js) as fh:
            results[key] = json.load(fh)
    with open(os.path.join(args.out, "log_plot_{}.txt".format(key)), "w", encoding="utf-8") as fh:
        subprocess.call([PY, os.path.join(SCRIPTS, "plot_geometry_trajectories.py"), npz, "--step-dir", args.step_dir, "--n-traj", str(args.n_traj),
                         "--out", os.path.join(args.out, "geometry_trajectories_{}.png".format(key)),
                         "--side", os.path.join(args.out, "geometry_side_{}.png".format(key))], cwd=SCRIPTS, stdout=fh, stderr=subprocess.STDOUT)

lines = ["# {} : {} particles, beam angle {:.0f} deg".format(os.path.basename(os.path.normpath(args.out)), args.n, args.phi), "",
         "STEP files `{}`; quads {:+.0f}/{:+.0f} V rotated {:g}/{:g} deg; spiral {:.1f} V; fields at {:.2f} mm.".format(
             args.step_dir, args.q1, args.q2, args.rotate_quads[0], args.rotate_quads[1], args.spiral_voltage, 1e3 * args.res), "",
         "| run | through housing [%] | exit plane [%] | housing exit [%] | spiral [%] | asymptotic z [mm] | asymptotic vert. angle [deg] |",
         "|---|---|---|---|---|---|---|"]
for key, label in (("e_nosc", "no space charge"), ("e_sc", "{} mA PyAMG, {:.1f} mm cells".format(args.current_ma, 1e3 * args.h))):
    m = results.get(key)
    if m:
        lines.append("| {} | **{:.1f}** | {:.1f} | {:.2f} | {:.1f} | {:+.2f} +- {:.2f} | {:+.2f} +- {:.2f} |".format(
            label, 100 * m["transmission_through_housing"], 100 * m["transmission_exit_plane"], 100 * m["lost_housing_after_exit"],
            100 * m.get("lost_spiral", float("nan")), m["z_asym_mean_mm"], m["z_asym_rms_mm"], m["vert_angle_asym_mean_deg"], m["vert_angle_asym_rms_deg"]))
lines += ["", "Figures: geometry_trajectories_<run>.png (isometric, -z up), geometry_side_<run>.png (side view, median plane drawn).",
          "Hand-off files (openPMD, see Docs/openpmd_plane_crossing_handoff.md): handoff_nosc.h5 / handoff_sc.h5 (+ _lab6d.h5)."]
with open(os.path.join(args.out, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
stamp("done; outputs in {}".format(args.out))
