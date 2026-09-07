"""Run a finished geometry once: direct BEM solve of the STEP files at the final voltages,
the full RFQ core beam through it without space charge (and, with sc=True, with PyAMG
space charge), exit metrics, 3-D geometry/trajectory figures, the openPMD hand-off files
and a short report.md.

The solve and the two bunch runs are separate processes (python -m
spyral_inflector.tracking.bem_reload / .bunch), so bempp, warp and the Poisson solver
never share a process and the two bunch runs go in parallel.

    python -m spyral_inflector.tracking.run_final --step-dir STEPS --voltages voltages.csv --state state.pickle \\
        --out OUT --q1 6575 --q2 -8450 --rotate-quads 16 24 --particles core.txt --bfield B.pickle --sc
"""
import argparse
import json
import os
import subprocess
import sys
import time

from .deck import read_voltages


def _stamp(msg):
    print("[{}] {}".format(time.strftime("%H:%M:%S"), msg), flush=True)


def _launch(name, cmd, log_path):
    fh = open(log_path, "w", encoding="utf-8")
    return name, subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT), fh


def _wait(jobs):
    for name, p, fh in jobs:
        p.wait()
        fh.close()
        if p.returncode != 0:
            raise RuntimeError("{} failed (exit {}), see {}".format(name, p.returncode, fh.name))
        _stamp("{} finished".format(name))


def run_final(step_dir, voltages, state, out, q1, q2, rotate_quads=(0.0, 0.0), spiral_voltage=None, particles=None, bfield=None,
              phi=90.0, n=43969, res=0.0025, sc=False, h=1.5e-3, current_ma=8.0, rf_mhz=32.8, handoff=True, n_traj=300,
              tag="final", python=None):
    """Solve + track + report for one geometry. voltages/state: the geometry build's
    voltages.csv and state.pickle (the spiral voltage defaults to the csv's anode voltage).
    Returns {"metrics": {run: exit metrics}, "report": path, ...}."""
    from .metrics import exit_metrics
    from .plots import electrode_meshes, plot_geometry_trajectories, plot_side_view

    python = python or sys.executable
    if particles is None or bfield is None:
        raise ValueError("particles and bfield are required")
    if spiral_voltage is None:
        spiral_voltage = abs(read_voltages(voltages)["SI_Anode"])
    os.makedirs(out, exist_ok=True)
    common_env = dict(os.environ, PYTHONUNBUFFERED="1")
    _stamp("run_final -> {}".format(out))

    # 1. direct BEM solve (reused if present)
    if not os.path.exists(os.path.join(out, "ef_itp_{}.pickle".format(tag))):
        cmd = [python, "-m", "spyral_inflector.tracking.bem_reload", "--step-dir", step_dir, "--voltages", voltages, "--state", state,
               "--out-dir", out, "--tag", tag, "--bfield", bfield, "--res", str(res), "--quad-voltages", str(q1), str(q2),
               "--rotate-quads", str(rotate_quads[0]), str(rotate_quads[1]), "--spiral-voltage", str(spiral_voltage)]
        _stamp("BEM solve at quads {:+.0f}/{:+.0f} V rotated {:g}/{:g} deg, spiral {:.1f} V".format(q1, q2, rotate_quads[0], rotate_quads[1], spiral_voltage))
        fh = open(os.path.join(out, "log_solve.txt"), "w", encoding="utf-8")
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=common_env)
        fh.close()
        if p.returncode != 0:
            raise RuntimeError("BEM solve failed (exit {}), see log_solve.txt".format(p.returncode))

    # 2. the core beam without and (optionally) with space charge, in parallel
    base = [python, "-m", "spyral_inflector.tracking.bunch", "--tag", tag, "--reload-dir", out, "--step-dir", step_dir, "--out-dir", out,
            "--particles", particles, "--bfield", bfield, "--phi", str(phi), "--n", str(n), "--current-ma", str(current_ma),
            "--rf-mhz", str(rf_mhz), "--record", str(max(n_traj, 1000))]
    jobs = []
    if not os.path.exists(os.path.join(out, "bunch_e_nosc.json")):
        cmd = base + ["--out-tag", "e_nosc"] + (["--save-openpmd", os.path.join(out, "handoff_nosc.h5"), "--save-mode", "both"] if handoff else [])
        jobs.append(_launch("e_nosc", cmd, os.path.join(out, "log_e_nosc.txt")))
    if sc and not os.path.exists(os.path.join(out, "bunch_e_sc.json")):
        cmd = base + ["--out-tag", "e_sc", "--sc", "--h", str(h), "--reference", os.path.join(out, "bunch_e_nosc.json")] + (
            ["--save-openpmd", os.path.join(out, "handoff_sc.h5"), "--save-mode", "both"] if handoff else [])
        jobs.append(_launch("e_sc", cmd, os.path.join(out, "log_e_sc.txt")))
    if jobs:
        _stamp("tracking {:,d} particles: {}".format(n, ", ".join(j[0] for j in jobs)))
    _wait(jobs)

    # 3. exit metrics and figures
    results = {}
    state_fn = os.path.join(out, "si_state_{}.pickle".format(tag))
    meshes = None
    for key in ("e_nosc", "e_sc"):
        npz = os.path.join(out, "bunch_{}.npz".format(key))
        if not os.path.exists(npz):
            continue
        m = exit_metrics(npz, state_fn)
        results[key] = m
        with open(os.path.join(out, "exit_metrics_{}.json".format(key)), "w") as fh:
            json.dump(m, fh, indent=2)
        if meshes is None:
            meshes = electrode_meshes(step_dir)
        title = "{}: {:,d} particles, transmission {:.1f} %".format(key, m["n"], 100 * m.get("transmission_through_housing", m["transmission"]))
        plot_geometry_trajectories(npz, step_dir, os.path.join(out, "geometry_trajectories_{}.png".format(key)), title=title, n_traj=n_traj, meshes=meshes)
        plot_side_view(npz, step_dir, os.path.join(out, "geometry_side_{}.png".format(key)), title=title, n_traj=n_traj, meshes=meshes)

    lines = ["# {}: {:,d} particles, beam angle {:.0f} deg".format(os.path.basename(os.path.normpath(out)), n, phi), "",
             "STEP files `{}`; quads {:+.0f}/{:+.0f} V rotated {:g}/{:g} deg; spiral {:.1f} V; fields at {:.2f} mm.".format(
                 step_dir, q1, q2, rotate_quads[0], rotate_quads[1], spiral_voltage, 1e3 * res), "",
             "| run | through housing [%] | exit plane [%] | housing exit [%] | spiral [%] | asymptotic z [mm] | asymptotic vert. angle [deg] |",
             "|---|---|---|---|---|---|---|"]
    for key, label in (("e_nosc", "no space charge"), ("e_sc", "{} mA PyAMG, {:.1f} mm cells".format(current_ma, 1e3 * h))):
        m = results.get(key)
        if m and "z_asym_mean_mm" in m:
            lines.append("| {} | **{:.1f}** | {:.1f} | {:.2f} | {:.1f} | {:+.2f} +- {:.2f} | {:+.2f} +- {:.2f} |".format(
                label, 100 * m["transmission_through_housing"], 100 * m["transmission_exit_plane"], 100 * m["lost_housing_after_exit"],
                100 * m["lost_spiral"], m["z_asym_mean_mm"], m["z_asym_rms_mm"], m["vert_angle_asym_mean_deg"], m["vert_angle_asym_rms_deg"]))
        elif m:
            lines.append("| {} | **{:.1f}** | | | {:.1f} | | |".format(label, 100 * m["transmission"], 100 * m["lost_spiral"]))
    lines += ["", "Figures: geometry_trajectories_<run>.png (isometric, -z up), geometry_side_<run>.png (side view, median plane drawn).",
              "Hand-off files (openPMD, PyPATools/documents/openpmd_plane_crossing_handoff.md): handoff_nosc.h5 / handoff_sc.h5 (+ _lab6d.h5)."]
    report = os.path.join(out, "report.md")
    with open(report, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)
    _stamp("done; outputs in {}".format(out))
    return {"metrics": results, "report": report, "out": out, "spiral_voltage": spiral_voltage}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--step-dir", required=True)
    p.add_argument("--voltages", required=True, help="voltages.csv of the geometry build")
    p.add_argument("--state", required=True, help="state.pickle of the geometry build")
    p.add_argument("--out", required=True)
    p.add_argument("--q1", type=float, required=True)
    p.add_argument("--q2", type=float, required=True)
    p.add_argument("--rotate-quads", type=float, nargs=2, default=[0.0, 0.0], metavar=("A1", "A2"))
    p.add_argument("--spiral-voltage", type=float, default=None, help="default: the anode voltage of voltages.csv")
    p.add_argument("--particles", required=True)
    p.add_argument("--bfield", required=True)
    p.add_argument("--phi", type=float, default=90.0)
    p.add_argument("--n", type=int, default=43969)
    p.add_argument("--res", type=float, default=0.0025)
    p.add_argument("--sc", action="store_true")
    p.add_argument("--h", type=float, default=1.5e-3)
    p.add_argument("--current-ma", type=float, default=8.0)
    p.add_argument("--rf-mhz", type=float, default=32.8)
    p.add_argument("--no-handoff", action="store_true")
    p.add_argument("--n-traj", type=int, default=300)
    p.add_argument("--tag", default="final")
    a = p.parse_args(argv)
    return run_final(a.step_dir, a.voltages, a.state, a.out, a.q1, a.q2, rotate_quads=tuple(a.rotate_quads), spiral_voltage=a.spiral_voltage,
                     particles=a.particles, bfield=a.bfield, phi=a.phi, n=a.n, res=a.res, sc=a.sc, h=a.h, current_ma=a.current_ma,
                     rf_mhz=a.rf_mhz, handoff=not a.no_handoff, n_traj=a.n_traj, tag=a.tag)


if __name__ == "__main__":
    main()
