"""System rotation scan by rigid rotation of an EXPORTED geometry (the generator's own
`rotation` knob produces geometries gmsh cannot mesh).

For every angle R: the STEP assembly (electrodes, housing, apertures, quads) and the design
orbit are rotated about the injection axis by R (TrackFromStep --rotate-all), the quads keep
their 16/24 deg on top, the beam angle becomes 90 + R; three basis fields are solved at
--res, the quad voltages retuned by superposition on a small bunch (BunchScan2) and a
larger bunch tracked at the best point (BunchTrack + exit metrics). Finished angles are
reused; the table goes to Results/rotation/<name>/report.md.

    python RotationScan2.py --name rot2 --angles -110 -105 -100 -95 -90 -85 -80 -60 -30
"""
import argparse
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exit_metrics import exit_metrics  # noqa: E402

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
PY = sys.executable

parser = argparse.ArgumentParser()
parser.add_argument("--name", default="rot2")
parser.add_argument("--angles", type=float, nargs="+", default=[-110, -105, -100, -95, -90, -85, -80, -60, -30])
parser.add_argument("--steps", default=os.path.join(DECK, "Geometry", "housing", "housing3_steps"), help="STEP folder of the base geometry")
parser.add_argument("--run-dir", default=os.path.join(DECK, "Results", "housing", "housing3", "geometry"), help="its voltages.csv and state.pickle")
parser.add_argument("--base-quads", type=float, nargs=2, default=[16.0, 24.0])
parser.add_argument("--base-phi", type=float, default=90.0)
parser.add_argument("--q1", type=float, nargs=3, default=[5075, 8075, 4], metavar=("MIN", "MAX", "N"))
parser.add_argument("--q2", type=float, nargs=3, default=[-9950, -5950, 5], metavar=("MIN", "MAX", "N"))
parser.add_argument("--res", type=float, default=0.005)
parser.add_argument("--n-scan", type=int, default=1500)
parser.add_argument("--n-bunch", type=int, default=5000)
parser.add_argument("--parallel", type=int, default=3)
parser.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit_core_as_txt.txt"))
parser.add_argument("--bfield", default=os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle"))
args = parser.parse_args()

out = os.path.join(DECK, "Results", "rotation", args.name)
os.makedirs(out, exist_ok=True)


def stamp(msg):
    print("[{}] {}".format(time.strftime("%H:%M:%S"), msg), flush=True)


def name_of(a):
    return "rot{:+g}".format(a).replace("+", "p").replace("-", "m")


def point_script(a):
    """A small python driver per angle (run as a subprocess so angles can go in parallel)."""
    d = os.path.join(out, name_of(a))
    os.makedirs(d, exist_ok=True)
    q1, q2 = args.base_quads[0], args.base_quads[1]
    code = r'''
import json, os, subprocess, sys, time
sys.path.insert(0, {scripts!r})
from exit_metrics import exit_metrics
PY, S, d, steps, run_dir, R = sys.executable, {scripts!r}, {d!r}, {steps!r}, {run_dir!r}, {angle!r}
def run(cmd, log):
    with open(os.path.join(d, log), "w", encoding="utf-8") as fh:
        p = subprocess.run([PY, "-u"] + cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=S)
    if p.returncode != 0:
        raise SystemExit("{{}} failed, see {{}}".format(cmd[0], log))
common = ["--step-dir", steps, "--voltages", os.path.join(run_dir, "voltages.csv"), "--state", os.path.join(run_dir, "state.pickle"),
          "--out-dir", d, "--res", {res!r}, "--bfield", {bfield!r}, "--rotate-all", str(R), "--rotate-quads", {q1!r}, {q2!r}]
t0 = time.time()
for tag, extra in (("spiral", ["--quad-voltages", "0", "0"]),
                   ("q1", ["--spiral-voltage", "0", "--quad-voltages", "3500", "0", "--no-test"]),
                   ("q2", ["--spiral-voltage", "0", "--quad-voltages", "0", "3500", "--no-test"])):
    if not os.path.exists(os.path.join(d, "ef_itp_" + tag + ".pickle")):
        run([os.path.join(S, "TrackFromStep.py"), "--tag", tag] + extra + common, "log_reload_" + tag + ".txt")
if not os.path.exists(os.path.join(d, "scan2_inner.json")):
    run([os.path.join(S, "BunchScan2.py"), "--reload-dir", d, "--out-dir", d, "--step-dir", steps, "--bfield", {bfield!r},
         "--basis", "spiral", "q1", "q1", "q2", "q2", "--phi", str({phi!r} + R), "--q1", {q1a!r}, {q1b!r}, {q1n!r}, "--q2", {q2a!r}, {q2b!r}, {q2n!r},
         "--particles", {particles!r}, "--n", {nscan!r}, "--tag", "inner"], "log_scan_inner.txt")
sc = json.load(open(os.path.join(d, "scan2_inner.json")))
best = sc["best"]
if not os.path.exists(os.path.join(d, "bunch_bunch.npz")):
    run([os.path.join(S, "BunchTrack.py"), "--tag", "inner_best", "--reload-dir", d, "--step-dir", steps, "--particles", {particles!r},
         "--bfield", {bfield!r}, "--phi", str({phi!r} + R), "--n", {nbunch!r}, "--out-tag", "bunch"], "log_bunch.txt")
m = exit_metrics(os.path.join(d, "bunch_bunch.npz"), os.path.join(d, "si_state_inner_best.pickle"))
rl = json.load(open(os.path.join(d, "reload_spiral.json")))
json.dump({{"rotation_deg": R, "inner_best": best, "bunch": m, "reload_test_particle": rl.get("test_particle"),
           "surface_field": rl.get("surface_field"), "wall_s": time.time() - t0}}, open(os.path.join(d, "summary.json"), "w"), indent=2)
print("done", R, m["transmission"])
'''.format(scripts=SCRIPTS, d=d, steps=args.steps, run_dir=args.run_dir, angle=float(a), res=str(args.res), bfield=args.bfield,
           q1=str(q1), q2=str(q2), phi=float(args.base_phi), q1a=str(args.q1[0]), q1b=str(args.q1[1]), q1n=str(int(args.q1[2])),
           q2a=str(args.q2[0]), q2b=str(args.q2[1]), q2n=str(int(args.q2[2])), particles=args.particles, nscan=str(args.n_scan), nbunch=str(args.n_bunch))
    fn = os.path.join(d, "point.py")
    with open(fn, "w", encoding="utf-8") as fh:
        fh.write(code)
    return fn


queue = [a for a in args.angles if not os.path.exists(os.path.join(out, name_of(a), "summary.json"))]
stamp("rotation scan '{}' (rigid rotation of {}): angles {} ({} to run)".format(args.name, os.path.basename(args.steps), args.angles, len(queue)))
running = []
while queue or running:
    while queue and len(running) < args.parallel:
        a = queue.pop(0)
        fh = open(os.path.join(out, "log_{}.txt".format(name_of(a))), "w", encoding="utf-8")
        running.append((a, subprocess.Popen([PY, "-u", point_script(a)], stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS), fh))
        stamp("start {:+g} deg".format(a))
    for item in list(running):
        a, p, fh = item
        if p.poll() is not None:
            fh.close()
            running.remove(item)
            stamp("{}: {:+g} deg (exit {})".format("done" if p.returncode == 0 else "FAILED", a, p.returncode))
    if running:
        time.sleep(20)

rows = []
for a in args.angles:
    fn = os.path.join(out, name_of(a), "summary.json")
    if not os.path.exists(fn):
        continue
    s = json.load(open(fn))
    b, tp = s["bunch"], s.get("reload_test_particle") or {}
    rows.append({"rotation_deg": a, "transmission": b.get("transmission_through_housing", b["transmission"]), "lost_spiral": b["lost_spiral"],
                 "lost_apertures": b["lost_apertures"], "q1": s["inner_best"]["q1"], "q2": s["inner_best"]["q2"],
                 "exit_azimuth_deg": b.get("design_azimuth_deg"), "z_asym_mm": b.get("z_asym_mean_mm"), "angle_asym_deg": b.get("vert_angle_asym_mean_deg"),
                 "design_outside_angle_deg": tp.get("outside_angle_new_deg"), "design_outside_z_mm": tp.get("outside_z_new_mm")})
lines = ["# System rotation scan (rigid): {}".format(args.name), "",
         "Base geometry `{}` rotated about the injection axis; quads keep {}/{} deg on top, beam angle {} + R. {} particles at the best quads, fields at {:.0f} mm (ranking only).".format(
             args.steps, args.base_quads[0], args.base_quads[1], args.base_phi, args.n_bunch, 1e3 * args.res), "",
         "| rotation [deg] | exit azimuth [deg] | transmission [%] | spiral | apertures | q1/q2 [V] | bunch asym z [mm] | bunch asym angle [deg] | design particle outside angle [deg] / z [mm] |",
         "|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    f = lambda v, fmt: (fmt.format(v) if v is not None else "-")  # noqa: E731
    lines.append("| {:+g} | {} | {:.1f} | {:.1f} | {:.1f} | {:+.0f}/{:+.0f} | {} | {} | {} / {} |".format(
        r["rotation_deg"], f(r["exit_azimuth_deg"], "{:.1f}"), 100 * r["transmission"], 100 * r["lost_spiral"], 100 * r["lost_apertures"],
        r["q1"], r["q2"], f(r["z_asym_mm"], "{:+.2f}"), f(r["angle_asym_deg"], "{:+.2f}"), f(r["design_outside_angle_deg"], "{:+.2f}"), f(r["design_outside_z_mm"], "{:+.2f}")))
with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
with open(os.path.join(out, "summary.json"), "w") as fh:
    json.dump({"args": vars(args), "rows": rows}, fh, indent=2)
stamp("done: " + ", ".join("{:+g} deg {:.1f} %".format(r["rotation_deg"], 100 * r["transmission"]) for r in rows))
