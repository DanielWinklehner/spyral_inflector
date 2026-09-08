"""Bunch-based fringe optimization of the spiral inflector (Daniel, 2026-09-07).

The single-particle optimizer sets the exit truncation, the axial shift dz and the spiral
voltage so that the DESIGN particle leaves level, centred and on the median plane. The bunch
wants something else: its losses sit on the anode in the last 15 mm of the electrodes, and a
higher voltage moves the bunch away from the anode but over-turns it, which a slightly
shorter exit undoes. This script optimizes those three knobs on the bunch:

  * objective   : transmission through the housing exit opening; within 1 point of the best,
                  the lowest spiral-electrode loss (anode + cathode) wins
  * constraints : asymptotic centroid (31 mm past the crossing, outside the fringe) level and
                  on the median plane; the vertical angle is zeroed by the spiral voltage at
                  each exit truncation, the height by dz afterwards (dz shifts it 1:1)
  * bound       : spiral electrode voltage <= --vmax (12 kV sparking limit)

For each exit truncation: build the geometry at the reference dz and reference voltage
(single-particle optimizer with every knob held, i.e. one evaluation), export STEP, solve the
three basis fields (spiral, quad 1, quad 2 with the given rotations) at --res, and scan the
spiral voltage scale by superposition on --n-scan particles. The voltage that zeroes the
asymptotic angle is interpolated, the losses there are read off, and the truncation with the
best transmission (spiral loss as tie-break) under the voltage cap wins. The winner is rebuilt with dz corrected by
minus its asymptotic height, solved directly with all voltages, and checked on --n-bunch
particles.

    python FringeOptBunch.py --name fr_pg5 --steps-from pg5b --q1 7150 --q2 -10225 --rotate-quads 8 16 ...
"""
import argparse
import csv
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from PyPATools.beam import ParticleDistribution  # noqa: E402
from PyPATools.species import IonSpecies  # noqa: E402
from spyral_inflector import SpiralInflector  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exit_metrics import exit_metrics, fmt as fmt_metrics  # noqa: E402

DECK = os.environ.get("SI_DECK") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SCRIPTS = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

parser = argparse.ArgumentParser()
parser.add_argument("--name", required=True)
parser.add_argument("--out-root", default=os.path.join(DECK, "Results", "fringe"))
parser.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit_core_as_txt.txt"))
parser.add_argument("--bfield", default=os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle"))
parser.add_argument("--phi", type=float, default=90.0)
parser.add_argument("--rotate-quads", type=float, nargs=2, default=[0.0, 0.0], metavar=("A1", "A2"))
parser.add_argument("--q1", type=float, required=True, help="quad 1 voltage [V] (fixed)")
parser.add_argument("--q2", type=float, required=True, help="quad 2 voltage [V] (fixed)")
parser.add_argument("--ent-trunc", type=float, default=0.34, help="entrance truncation [deg], frozen")
parser.add_argument("--exit-truncs", type=float, nargs="+", default=[0.77, 1.5, 2.5, 3.5, 5.0], help="exit truncations to try [deg]")
parser.add_argument("--dz-ref", type=float, required=True, help="reference axial shift of the whole system [m]")
parser.add_argument("--v-ref", type=float, required=True, help="reference spiral voltage [V] (the chain's), vscale = 1 there")
parser.add_argument("--vscales", type=float, nargs="+", default=[0.97, 0.985, 1.0, 1.01, 1.02, 1.03, 1.045, 1.06])
parser.add_argument("--vmax", type=float, default=12000.0, help="spiral electrode voltage limit [V]")
parser.add_argument("--level-iters", type=int, default=3,
                    help="re-levelling rounds on the shifted geometry: the dz shift changes the exit angle (~0.3 deg/mm), so "
                         "voltage and dz are iterated until the asymptotic angle and height are within tolerance (0 = single shift)")
parser.add_argument("--tol-angle", type=float, default=0.15, help="re-levelling converged when |asymptotic angle| < this [deg] ...")
parser.add_argument("--tol-z", type=float, default=0.3, help="... and |asymptotic z| < this [mm]")
parser.add_argument("--iter-vscales", type=float, nargs="+", default=[0.97, 0.985, 1.0, 1.015, 1.03],
                    help="voltage scales relative to the round's voltage scanned in each re-levelling round")
parser.add_argument("--n-scan", type=int, default=1500)
parser.add_argument("--n-bunch", type=int, default=10000)
parser.add_argument("--res", type=float, default=0.0025)
parser.add_argument("--volt", type=float, default=12000.0, help="design voltage of the analytic geometry [V]")
parser.add_argument("--gap", type=float, default=0.019)
parser.add_argument("--tilt", type=float, default=31.0)
parser.add_argument("--dx", type=float, default=0.01)
parser.add_argument("--sigma", type=float, default=0.0022)
parser.add_argument("--aspect", type=float, default=2.4)
parser.add_argument("--gamma", type=float, default=5.0)
parser.add_argument("--angling", type=float, default=11.0)
parser.add_argument("--quad-bore", type=float, default=0.013)
parser.add_argument("--quad-z1", type=float, default=-0.27)
parser.add_argument("--quad-z2", type=float, default=-0.19)
parser.add_argument("--quad-len", type=float, default=0.045)
parser.add_argument("--quad-len2", type=float, default=None)
parser.add_argument("--shared-plates", action="store_true")
parser.add_argument("--aper-hole", type=float, default=0.0125)
parser.add_argument("--entrance-hole", type=float, default=None)
parser.add_argument("--plate-gap", type=float, default=0.001)
parser.add_argument("--plate-thickness", type=float, default=0.005)
parser.add_argument("--slot-width", type=float, default=0.015)
parser.add_argument("--slot-length", type=float, default=0.040)
parser.add_argument("--exit-opening", type=float, nargs=2, default=None, metavar=("ACROSS", "ALONG"))
parser.add_argument("--energy-mev", type=float, default=None)
args = parser.parse_args()
if args.quad_len2 is None:
    args.quad_len2 = args.quad_len
if args.shared_plates:
    args.quad_z2 = args.quad_z1 + args.quad_len + 2.0 * args.plate_gap + args.plate_thickness
if args.energy_mev is None:
    args.energy_mev = float(np.mean(np.loadtxt(args.particles, skiprows=1, usecols=8)))

out_dir = os.path.join(args.out_root, args.name)
os.makedirs(out_dir, exist_ok=True)
t_start = time.time()
vs_ref = args.v_ref / args.volt          # scale of the reference voltage relative to the design voltage
summary = {"name": args.name, "args": vars(args), "energy_mev": args.energy_mev, "vs_ref": vs_ref, "cases": []}


def log(step, msg):
    print("[{:>5.0f} s] {}: {}".format(time.time() - t_start, step, msg), flush=True)


def run(cmd, logname, cwd):
    with open(os.path.join(cwd, logname), "w", encoding="utf-8") as fh:
        p = subprocess.run([PY, "-u"] + cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=cwd)
    if p.returncode != 0:
        raise RuntimeError("{} failed (exit {}), see {}".format(os.path.basename(cmd[0]), p.returncode, logname))


def sanitize(name):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(name)).strip("_")


def build_and_export(tag, exit_trunc, dz, vscale):
    """Build the inflector at fixed (entrance, exit truncation, dz, voltage scale), export STEP.
    Returns (case_dir, steps_dir, state_fn, volt_fn, spiral_V)."""
    case_dir = os.path.join(out_dir, tag)
    os.makedirs(case_dir, exist_ok=True)
    steps_dir = os.path.join(DECK, "Geometry", "wiggle", "{}_{}_steps".format(args.name, tag))
    ion = ParticleDistribution(species=IonSpecies("H2_1+"))
    ion.set_mean_energy_z_mev(args.energy_mev)
    si = SpiralInflector(ion=ion, method="numerical", solver="bempp", volt=args.volt, gap=args.gap, tilt=args.tilt,
                         dx=args.dx, sigma=args.sigma, vee_shape="parabolic", ns=100, aspect_ratio=args.aspect,
                         rotation=0.0, debug=False, gammaAng=args.gamma, anglingAng=args.angling)
    si.load_bfield(bfield=args.bfield)
    si.initialize()
    si.set_parameter(key="h", value=0.005)
    si.set_parameter(key="make_aperture", value=True)
    si.set_parameter(key="aperture_params", value={"thickness": 4e-3, "radius": 50e-3, "length": args.slot_length, "width": args.slot_width,
                                                   "top_distance": 5e-3, "bottom_distance": 10e-3, "hole_type": "rectangle", "voltage": 0.0,
                                                   **({"exit_width": args.exit_opening[0], "exit_length": args.exit_opening[1]} if args.exit_opening else {})})
    si.set_parameter(key="make_housing", value=True)
    si.set_parameter(key="housing_params", value={"zmin": -0.12, "zmax": 0.03, "span": True, "gap": 6e-3, "thickness": 4e-3, "voltage": 0.0,
                                                  "experimental": True})
    si.set_parameter(key="make_quadrupoles", value=True)
    si.set_parameter(key="quadrupole_params", value={"a": args.quad_bore, "b": args.quad_bore, "radius": 0.04,
                                                     "z_starts": [args.quad_z1, args.quad_z2], "lengths": [args.quad_len, args.quad_len2],
                                                     "voltages": [3500, 3500], "aper_rad": 2.0 * args.aper_hole,
                                                     "plate_gap": args.plate_gap, "plate_thickness": args.plate_thickness,
                                                     "shared_plates": bool(args.shared_plates),
                                                     "entrance_aper_rad": 2.0 * (args.entrance_hole if args.entrance_hole else args.aper_hole)})
    si.generate_geometry()
    # one evaluation at the held knobs (no optimization): geometry, BEM solve, design-particle check
    result = si.optimize_trajectory(maxiter=0, solver="dfols", res=0.005, vary_voltage=False,
                                    initial_guess=[args.ent_trunc, exit_trunc, dz, vscale],
                                    fixed={0: args.ent_trunc, 1: exit_trunc, 2: dz},
                                    bounds=((0.0, 15.0), (0.0, 15.0), (-15.0e-3, 15.0e-3), (0.85, 1.3)),
                                    exclude_quadrupoles=True)
    m = result["measurements"]
    assembly = si.numerical_variables["objects"]
    names, used = {}, set()
    for i, e in enumerate(assembly.electrodes.values()):
        nm = sanitize(e.name) or "electrode_{:03d}".format(i)
        u, k = nm, 2
        while u.lower() in used:
            u = "{}_{}".format(nm, k)
            k += 1
        used.add(u.lower())
        names[i] = u
    electrode_voltages = {names[i]: float(e.voltage) for i, e in enumerate(assembly.electrodes.values())}
    mesh = si.numerical_variables["full mesh"]
    ctx = result["context"]
    state = {"trj_design": si.analytic_variables["trj_design"] + result["shift"], "v_design": si.analytic_variables["v_design"],
             "trj_design_full": si.analytic_variables["trj_design_full"] + result["shift"], "v_design_full": si.analytic_variables["v_design_full"],
             "b_lim_deg": result["b_lim_deg"], "shift": result["shift"], "shift_lab": result["shift_lab"],
             "volt_scale": result["volt_scale"], "voltage": result["voltage"], "electrode_voltages": electrode_voltages,
             "z_start": ctx["z_start"], "v0": ctx["v0"], "dt": ctx["dt"], "nsteps": ctx["nsteps"], "track_r": m["r"], "track_v": m["v"],
             "mesh": {"verts": np.asarray(mesh["verts"]), "elems": np.asarray(mesh["elems"]), "domns": np.asarray(mesh["domns"]),
                      "names": {int(e.bempp_domain): names[i] for i, e in enumerate(assembly.electrodes.values())}},
             "residual_final": np.asarray(result["residual_final"]), "res_m": 0.005,
             "knobs": {"ent_trunc": args.ent_trunc, "exit_trunc": exit_trunc, "dz": dz, "vscale": vscale}}
    state_fn = os.path.join(case_dir, "state.pickle")
    with open(state_fn, "wb") as fh:
        pickle.dump(state, fh)
    volt_fn = os.path.join(case_dir, "voltages.csv")
    with open(volt_fn, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Electrode", "Voltage (V)"])
        for nm, v in electrode_voltages.items():
            w.writerow([nm, "{:.3f}".format(v)])
    if os.path.isdir(steps_dir):
        shutil.rmtree(steps_dir)
    os.makedirs(steps_dir)
    for i, e in enumerate(assembly.electrodes.values()):
        if e.export(os.path.join(steps_dir, "{:03d}_{}.step".format(i, names[i]))) != 0:
            raise RuntimeError("STEP export failed for {}".format(e.name))
    design = {"exit_angle_deg": float(result["residual_final"][0]), "centering_mm": 1e3 * float(result["residual_final"][1]),
              "z_offset_mm": 1e3 * float(result["residual_final"][2]), "width_mm": 1e3 * float(result["residual_final"][3]),
              "spiral_V": float(result["voltage"]), "clearance_anode_mm": 1e3 * m["clearance_anode_exit"],
              "clearance_cathode_mm": 1e3 * m["clearance_cathode_exit"]}
    del si
    return case_dir, steps_dir, state_fn, volt_fn, design


def basis_fields(case_dir, steps_dir, state_fn, volt_fn):
    common = ["--step-dir", steps_dir, "--voltages", volt_fn, "--state", state_fn, "--out-dir", case_dir, "--res", str(args.res),
              "--bfield", args.bfield, "--rotate-quads", str(args.rotate_quads[0]), str(args.rotate_quads[1])]
    run([os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "spiral", "--quad-voltages", "0", "0"] + common, "log_reload_spiral.txt", case_dir)
    run([os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "q1", "--spiral-voltage", "0", "--quad-voltages", "3500", "0", "--no-test"] + common, "log_reload_q1.txt", case_dir)
    run([os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "q2", "--spiral-voltage", "0", "--quad-voltages", "0", "3500", "--no-test"] + common, "log_reload_q2.txt", case_dir)


def interp_zero(x, y):
    """x where the piecewise-linear y(x) crosses zero (first crossing), or None."""
    for i in range(len(x) - 1):
        if y[i] == 0:
            return x[i]
        if y[i] * y[i + 1] < 0:
            return x[i] + (x[i + 1] - x[i]) * (-y[i]) / (y[i + 1] - y[i])
    return None


# ------------------------------------------------------------------ exit-truncation cases
for et in args.exit_truncs:
    tag = "t{:.2f}".format(et).replace(".", "p")
    case_dir = os.path.join(out_dir, tag)
    steps_dir = os.path.join(DECK, "Geometry", "wiggle", "{}_{}_steps".format(args.name, tag))
    scan_fn = os.path.join(case_dir, "scan2_vscan.json")
    prev = None
    if os.path.exists(scan_fn) and os.path.isdir(steps_dir) and os.path.exists(os.path.join(case_dir, "ef_itp_spiral.pickle")):
        prev = json.load(open(scan_fn))
        have = sorted(r["vscale"] for r in prev["results"])
        missing = [v for v in args.vscales if not any(abs(v - h) < 1e-6 for h in have)]
        design = json.load(open(os.path.join(case_dir, "design.json"))) if os.path.exists(os.path.join(case_dir, "design.json")) else {}
        log(tag, "existing case reused ({} scanned voltages{})".format(len(have), ", scanning {} more".format(len(missing)) if missing else ""))
        if missing:
            run([os.path.join(SCRIPTS, "BunchScan2.py"), "--reload-dir", case_dir, "--out-dir", case_dir, "--step-dir", steps_dir, "--bfield", args.bfield,
                 "--basis", "spiral", "q1", "q1", "q2", "q2", "--phi", str(args.phi), "--points", str(args.q1), str(args.q2),
                 "--vscale"] + [str(v) for v in missing] + ["--particles", args.particles, "--n", str(args.n_scan), "--tag", "vscan_more"], "log_vscan_more.txt", case_dir)
            more = json.load(open(os.path.join(case_dir, "scan2_vscan_more.json")))
            prev["results"] += more["results"]
            json.dump(prev, open(scan_fn, "w"), indent=2)
        sc = prev
    else:
        log(tag, "exit truncation {:.2f} deg: build at dz {:+.2f} mm, spiral {:.0f} V".format(et, 1e3 * args.dz_ref, args.v_ref))
        case_dir, steps_dir, state_fn, volt_fn, design = build_and_export(tag, et, args.dz_ref, vs_ref)
        json.dump(design, open(os.path.join(case_dir, "design.json"), "w"), indent=2)
        log(tag, "design particle: exit angle {:+.3f} deg, centring {:+.2f} mm, z {:+.2f} mm, clearance A/C {:.2f}/{:.2f} mm".format(
            design["exit_angle_deg"], design["centering_mm"], design["z_offset_mm"], design["clearance_anode_mm"], design["clearance_cathode_mm"]))
        basis_fields(case_dir, steps_dir, state_fn, volt_fn)
        run([os.path.join(SCRIPTS, "BunchScan2.py"), "--reload-dir", case_dir, "--out-dir", case_dir, "--step-dir", steps_dir, "--bfield", args.bfield,
             "--basis", "spiral", "q1", "q1", "q2", "q2", "--phi", str(args.phi), "--points", str(args.q1), str(args.q2),
             "--vscale"] + [str(v) for v in args.vscales] + ["--particles", args.particles, "--n", str(args.n_scan), "--tag", "vscan"], "log_vscan.txt", case_dir)
        sc = json.load(open(scan_fn))
    pts = sorted(sc["results"], key=lambda r: r["vscale"])
    vs = np.array([r["vscale"] for r in pts])
    ang = np.array([r["vert_angle_asym_mean_deg"] if r["vert_angle_asym_mean_deg"] is not None else np.nan for r in pts])
    zas = np.array([r["z_asym_mean_mm"] if r["z_asym_mean_mm"] is not None else np.nan for r in pts])
    tr = np.array([r["transmission"] for r in pts])
    sp = np.array([r["lost_spiral"] for r in pts])
    vs_zero = interp_zero(vs, ang)
    case = {"exit_trunc_deg": et, "design": design, "vscales": vs.tolist(), "transmission": tr.tolist(), "lost_spiral": sp.tolist(),
            "angle_asym_deg": ang.tolist(), "z_asym_mm": zas.tolist(), "steps_dir": steps_dir, "case_dir": case_dir}
    fallback = vs_zero is None
    if fallback and np.isfinite(ang).any():
        # no zero crossing in the scanned range: take the scanned voltage closest to level (within the cap)
        ok = np.isfinite(ang) & (args.v_ref * vs <= args.vmax)
        if ok.any():
            vs_zero = float(vs[ok][np.argmin(np.abs(ang[ok]))])
            log(tag, "WARNING: no zero crossing of the asymptotic angle ({:+.2f} .. {:+.2f} deg); nearest level point taken".format(
                np.nanmin(ang), np.nanmax(ang)))
    if vs_zero is not None:
        V = args.v_ref * vs_zero
        capped = V > args.vmax
        if capped:
            vs_zero = args.vmax / args.v_ref
            V = args.vmax
        case.update({"vs_level": float(vs_zero), "spiral_V_level": float(V), "capped": bool(capped), "fallback": bool(fallback),
                     "transmission_level": float(np.interp(vs_zero, vs, tr)), "lost_spiral_level": float(np.interp(vs_zero, vs, sp)),
                     "angle_level_deg": float(np.interp(vs_zero, vs, ang)), "z_level_mm": float(np.interp(vs_zero, vs, zas))})
        log(tag, "level beam at vscale {:.4f} ({:.0f} V{}): transmission {:.1f} %, spiral loss {:.2f} %, asymptotic z {:+.2f} mm".format(
            vs_zero, V, ", CAPPED at the limit" if capped else "", 100 * case["transmission_level"], 100 * case["lost_spiral_level"], case["z_level_mm"]))
    else:
        log(tag, "no zero crossing of the asymptotic angle in the scanned voltage range ({:+.2f} .. {:+.2f} deg)".format(ang[0], ang[-1]))
    summary["cases"].append(case)
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

# ------------------------------------------------------------------ choose and confirm
valid = [c for c in summary["cases"] if "vs_level" in c and not c["capped"] and not c["fallback"]]
if not valid:
    valid = [c for c in summary["cases"] if "vs_level" in c and not c["capped"]]
if not valid:
    valid = [c for c in summary["cases"] if "vs_level" in c]
if not valid:
    raise SystemExit("no case produced a usable level-beam point; see summary.json")
# overall losses first; within TIE of the best transmission, the lowest spiral-electrode loss wins
TIE = 0.01
top = max(c["transmission_level"] for c in valid)
best = min([c for c in valid if c["transmission_level"] >= top - TIE], key=lambda c: (c["lost_spiral_level"], -c["transmission_level"]))
dz_final = args.dz_ref - 1e-3 * best["z_level_mm"]
vs_final = best["vs_level"]                      # relative to v_ref
log("choice", "exit truncation {:.2f} deg, vscale {:.4f} ({:.0f} V), dz {:+.2f} mm (was {:+.2f}): spiral loss {:.2f} %, transmission {:.1f} %".format(
    best["exit_trunc_deg"], vs_final, best["spiral_V_level"], 1e3 * dz_final, 1e3 * args.dz_ref, 100 * best["lost_spiral_level"], 100 * best["transmission_level"]))

# --- re-level on the shifted geometry. Shifting the whole system by dz changes the exit angle
#     (about 0.3 deg per mm: the electrodes move relative to the magnetic field), so the voltage
#     that levelled the beam at dz_ref no longer does at dz_final. Each round builds the shifted
#     geometry, scans the voltage around the current one by superposition, and updates voltage
#     (angle -> 0) and dz (height -> 0) until both are within tolerance.
history = []
for it in range(1, args.level_iters + 1):
    tag = "final_it{}".format(it)
    c_dir, s_dir, st_fn, v_fn, dsg = build_and_export(tag, best["exit_trunc_deg"], dz_final, vs_final * vs_ref)
    basis_fields(c_dir, s_dir, st_fn, v_fn)
    run([os.path.join(SCRIPTS, "BunchScan2.py"), "--reload-dir", c_dir, "--out-dir", c_dir, "--step-dir", s_dir, "--bfield", args.bfield,
         "--basis", "spiral", "q1", "q1", "q2", "q2", "--phi", str(args.phi), "--points", str(args.q1), str(args.q2),
         "--vscale"] + [str(v) for v in args.iter_vscales] + ["--particles", args.particles, "--n", str(args.n_scan), "--tag", "relevel"],
        "log_relevel.txt", c_dir)
    rl = json.load(open(os.path.join(c_dir, "scan2_relevel.json")))
    pts = sorted(rl["results"], key=lambda r: r["vscale"])
    s = np.array([r["vscale"] for r in pts])                       # relative to this round's voltage
    ang = np.array([r["vert_angle_asym_mean_deg"] if r["vert_angle_asym_mean_deg"] is not None else np.nan for r in pts])
    zas = np.array([r["z_asym_mean_mm"] if r["z_asym_mean_mm"] is not None else np.nan for r in pts])
    tr = np.array([r["transmission"] for r in pts])
    sp = np.array([r["lost_spiral"] for r in pts])
    ang_now, z_now = float(np.interp(1.0, s, ang)), float(np.interp(1.0, s, zas))
    rec = {"iteration": it, "dz_m": dz_final, "vscale": vs_final, "spiral_V": args.v_ref * vs_final, "angle_deg": ang_now, "z_mm": z_now,
           "transmission": float(np.interp(1.0, s, tr)), "lost_spiral": float(np.interp(1.0, s, sp)), "case_dir": c_dir,
           "scan": {"vscales_rel": s.tolist(), "angle_deg": ang.tolist(), "z_mm": zas.tolist(), "transmission": tr.tolist(), "lost_spiral": sp.tolist()}}
    history.append(rec)
    summary["relevel"] = history
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    log(tag, "at dz {:+.2f} mm, {:.0f} V: asymptotic angle {:+.2f} deg, z {:+.2f} mm; transmission {:.1f} %, spiral loss {:.2f} %".format(
        1e3 * dz_final, args.v_ref * vs_final, ang_now, z_now, 100 * rec["transmission"], 100 * rec["lost_spiral"]))
    if abs(ang_now) < args.tol_angle and abs(z_now) < args.tol_z:
        log(tag, "converged (|angle| < {:.2f} deg, |z| < {:.2f} mm)".format(args.tol_angle, args.tol_z))
        break
    s_zero = interp_zero(s, ang)
    if s_zero is None:
        s_zero = float(s[np.nanargmin(np.abs(ang))])
        log(tag, "WARNING: no zero crossing of the angle in the re-level scan ({:+.2f} .. {:+.2f} deg); nearest point taken".format(np.nanmin(ang), np.nanmax(ang)))
    z_at = float(np.interp(s_zero, s, zas))
    vs_new = vs_final * s_zero
    if args.v_ref * vs_new > args.vmax:
        vs_new = args.vmax / args.v_ref
        log(tag, "voltage capped at the {:.0f} V limit".format(args.vmax))
    log(tag, "re-level: voltage {:.0f} -> {:.0f} V, dz {:+.2f} -> {:+.2f} mm".format(
        args.v_ref * vs_final, args.v_ref * vs_new, 1e3 * dz_final, 1e3 * (dz_final - 1e-3 * z_at)))
    vs_final = vs_new
    dz_final = dz_final - 1e-3 * z_at
best_V = args.v_ref * vs_final
case_dir, steps_dir, state_fn, volt_fn, design = build_and_export("final", best["exit_trunc_deg"], dz_final, vs_final * vs_ref)
common = ["--step-dir", steps_dir, "--voltages", volt_fn, "--state", state_fn, "--out-dir", case_dir, "--res", str(args.res), "--bfield", args.bfield]
run([os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "direct", "--quad-voltages", str(args.q1), str(args.q2),
     "--rotate-quads", str(args.rotate_quads[0]), str(args.rotate_quads[1]), "--spiral-voltage", str(best_V)] + common, "log_direct.txt", case_dir)
run([os.path.join(SCRIPTS, "BunchTrack.py"), "--tag", "direct", "--reload-dir", case_dir, "--step-dir", steps_dir, "--particles", args.particles,
     "--bfield", args.bfield, "--phi", str(args.phi), "--n", str(args.n_bunch), "--out-tag", "bunch"], "log_bunch.txt", case_dir)
metrics = exit_metrics(os.path.join(case_dir, "bunch_bunch.npz"), os.path.join(case_dir, "si_state_direct.pickle"))
with open(os.path.join(out_dir, "exit_metrics_final.json"), "w") as fh:
    json.dump(metrics, fh, indent=2)
summary["final"] = {"exit_trunc_deg": best["exit_trunc_deg"], "vscale": vs_final, "spiral_V": best_V, "dz_m": dz_final,
                    "scan_choice": {"vscale": best["vs_level"], "spiral_V": best["spiral_V_level"], "dz_m": args.dz_ref - 1e-3 * best["z_level_mm"]},
                    "relevel_rounds": len(history), "design": design, "steps_dir": steps_dir, "metrics": metrics, "wall_s": time.time() - t_start}
with open(os.path.join(out_dir, "summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)
log("final", fmt_metrics(metrics))

# ------------------------------------------------------------------ figure
fig, axs = plt.subplots(1, 3, figsize=(16, 5))
for c in summary["cases"]:
    lab = "exit trunc {:.2f} deg".format(c["exit_trunc_deg"])
    V = args.v_ref * np.array(c["vscales"])
    axs[0].plot(V, 100 * np.array(c["transmission"]), "o-", label=lab)
    axs[1].plot(V, 100 * np.array(c["lost_spiral"]), "o-", label=lab)
    axs[2].plot(V, c["angle_asym_deg"], "o-", label=lab)
for ax in axs:
    ax.axvline(args.vmax, color="tab:red", ls="--", lw=1)
    ax.set_xlabel("spiral electrode voltage [V]")
    ax.grid(alpha=0.3)
axs[0].set_ylabel("transmission through the housing [%]")
axs[1].set_ylabel("spiral-electrode losses [%]")
axs[2].set_ylabel("asymptotic vertical angle of the centroid [deg]")
axs[2].axhline(0, color="k", lw=0.8)
axs[0].legend(fontsize=8)
fig.suptitle("bunch-based fringe optimization '{}': chosen exit truncation {:.2f} deg at {:.0f} V (scan) -> {:.0f} V, dz {:+.2f} mm after {} re-levelling round(s)".format(
    args.name, best["exit_trunc_deg"], best["spiral_V_level"], best_V, 1e3 * dz_final, len(history)))
fig.tight_layout()
fig.savefig(os.path.join(out_dir, "fringe_opt.png"), dpi=140)
print("wrote {}".format(os.path.join(out_dir, "summary.json")), flush=True)
