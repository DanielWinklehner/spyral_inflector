"""Space-charge-aware quad retune on a finished chain geometry.

Grid of quad rotations x quad voltages, each point a BunchTrackSC run (PyAMG space charge,
default 8 mA, 2 mm cells, 10k particles) on the vacuum field superposed from the basis solves
of the run folder (spiral, q1, q1skew, q2, q2skew -- the skew ones are solved first if
missing). Points run in parallel; finished points are reused on a relaunch. The best point
(transmission through the housing, spiral loss as the tie-break within 1 point) is confirmed
with a direct BEM solve and the full core beam with and without space charge (the SC run at
--confirm-h cells), and both confirmation runs write the openPMD hand-off files.

    python SCRetune.py --name pg5L --run-dir ..\\Results\\wiggle\\pg5Lb --step-dir ..\\Geometry\\wiggle\\pg5Lb_steps \\
        --alphas 8,16 16,24 24,32 --q1 5825 6575 7325 --q2 -9200 -8450 -7700
"""
import argparse
import itertools
import json
import os
import subprocess
import sys
import time

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)

parser = argparse.ArgumentParser()
parser.add_argument("--name", required=True)
parser.add_argument("--run-dir", required=True, help="chain run folder with the basis solves and the d_direct state")
parser.add_argument("--step-dir", required=True)
parser.add_argument("--state-tag", default="d_direct", help="si_state_<tag>.pickle of the design orbit")
parser.add_argument("--out-root", default=os.path.join(DECK, "Results", "scretune"))
parser.add_argument("--alphas", nargs="+", default=["8,16", "16,24", "24,32"], help="quad rotations a1,a2 [deg]")
parser.add_argument("--q1", type=float, nargs="+", default=[5825, 6575, 7325])
parser.add_argument("--q2", type=float, nargs="+", default=[-9200, -8450, -7700])
parser.add_argument("--vscale", type=float, default=1.0)
parser.add_argument("--n", type=int, default=10000)
parser.add_argument("--h", type=float, default=0.002)
parser.add_argument("--parallel", type=int, default=3)
parser.add_argument("--current-ma", type=float, default=8.0)
parser.add_argument("--phi", type=float, default=90.0)
parser.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit_core_as_txt.txt"))
parser.add_argument("--bfield", default=os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle"))
parser.add_argument("--res", type=float, default=0.0025, help="BEM potential grid of the basis / confirmation solves [m]")
parser.add_argument("--confirm-n", type=int, default=43969)
parser.add_argument("--confirm-h", type=float, default=0.0015)
parser.add_argument("--skip-confirm", action="store_true")
parser.add_argument("--reference", default=None, help="exit_metrics json of the chain's own SC run, for the report")
args = parser.parse_args()

PY = sys.executable
out = os.path.join(args.out_root, args.name)
points_dir = os.path.join(out, "points")
os.makedirs(points_dir, exist_ok=True)
T0 = time.time()


def stamp(msg):
    print("[{}] {}".format(time.strftime("%H:%M:%S"), msg), flush=True)


def wait_all(procs, poll=20):
    """procs: list of (name, Popen, logfile handle); returns the failed names."""
    failed = []
    while procs:
        for item in list(procs):
            name, p, fh = item
            if p.poll() is not None:
                fh.close()
                procs.remove(item)
                if p.returncode != 0:
                    failed.append(name)
                    stamp("FAILED: {} (exit {})".format(name, p.returncode))
                else:
                    stamp("done: {}".format(name))
        if procs:
            time.sleep(poll)
    return failed


def launch(name, cmd, log_path, cwd=SCRIPTS):
    fh = open(log_path, "w", encoding="utf-8")
    p = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=cwd)
    return (name, p, fh)


# ------------------------------------------------------------------ 1. basis solves
# All five basis fields are solved into out/basis with UNROTATED quads (the run folder's own
# spiral/q1/q2 solves carry the chain's final rotation baked in, which the cos 2a / sin 2a
# superposition with the 45-degree skew solves cannot use).
basis_dir = os.path.join(out, "basis")
os.makedirs(basis_dir, exist_ok=True)
stamp("retune '{}': run folder {}, {} rotations x {} q1 x {} q2 = {} points, {} in parallel".format(
    args.name, args.run_dir, len(args.alphas), len(args.q1), len(args.q2), len(args.alphas) * len(args.q1) * len(args.q2), args.parallel))
common = ["--step-dir", args.step_dir, "--voltages", os.path.join(args.run_dir, "voltages.csv"),
          "--state", os.path.join(args.run_dir, "state.pickle"), "--out-dir", basis_dir, "--res", str(args.res), "--bfield", args.bfield]
BASIS = (("spiral", ["--quad-voltages", "0", "0", "--rotate-quads", "0", "0"]),
         ("q1", ["--spiral-voltage", "0", "--quad-voltages", "3500", "0", "--rotate-quads", "0", "0", "--no-test"]),
         ("q2", ["--spiral-voltage", "0", "--quad-voltages", "0", "3500", "--rotate-quads", "0", "0", "--no-test"]),
         ("q1skew", ["--spiral-voltage", "0", "--quad-voltages", "3500", "0", "--rotate-quads", "45", "0", "--no-test"]),
         ("q2skew", ["--spiral-voltage", "0", "--quad-voltages", "0", "3500", "--rotate-quads", "0", "45", "--no-test"]))
procs = []
for tag, extra in BASIS:
    if not os.path.exists(os.path.join(basis_dir, "ef_itp_{}.pickle".format(tag))):
        stamp("solving the {} basis field (unrotated quads)".format(tag))
        procs.append(launch(tag, [PY, "-u", os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", tag] + extra + common,
                            os.path.join(out, "log_basis_{}.txt".format(tag))))
if wait_all(procs):
    sys.exit("basis solve failed")
with open(os.path.join(basis_dir, "reload_spiral.json"), encoding="utf-8-sig") as fh:
    v_ref = float(json.load(fh)["voltages"]["SI_Anode"])
stamp("basis fields in {}; spiral basis {:.1f} V".format(basis_dir, v_ref))

# consistency check of the superposition: the grid centre without space charge, 2000 particles
alphas = [tuple(float(x) for x in a.split(",")) for a in args.alphas]
centre = (alphas[len(alphas) // 2], args.q1[len(args.q1) // 2], args.q2[len(args.q2) // 2])
chk = os.path.join(out, "bunch_check_vacuum.json")
if not os.path.exists(chk):
    a, q1, q2 = centre
    cmd = [PY, "-u", os.path.join(SCRIPTS, "BunchTrackSC.py"), "--tag", args.state_tag, "--reload-dir", args.run_dir, "--step-dir", args.step_dir,
           "--out-dir", out, "--out-tag", "check_vacuum", "--no-sc", "--superpose", str(q1), str(a[0]), str(q2), str(a[1]),
           "--vscale", str(args.vscale), "--basis-dir", basis_dir, "--n", "2000", "--phi", str(args.phi), "--particles", args.particles,
           "--bfield", args.bfield, "--record", "10"]
    if wait_all([launch("check_vacuum", cmd, os.path.join(out, "log_check_vacuum.txt"))]):
        sys.exit("vacuum check failed")
with open(chk, encoding="utf-8") as fh:
    d = json.load(fh)
stamp("vacuum check at alpha {:g}/{:g}, q {:+.0f}/{:+.0f} by superposition: {:.1f} % of 2000 (expect the chain's ~90 %){}".format(
    centre[0][0], centre[0][1], centre[1], centre[2], 100 * d["transmission"], "" if d["transmission"] > 0.8 else "   <-- SUPERPOSITION SUSPECT"))
if d["transmission"] <= 0.8:
    sys.exit("superposed vacuum field does not reproduce the chain result; check the basis solves")

# ------------------------------------------------------------------ 2. the grid
alphas = [tuple(float(x) for x in a.split(",")) for a in args.alphas]
grid = list(itertools.product(alphas, args.q1, args.q2))


def point_tag(a, q1, q2):
    return "a{:g}_{:g}_q{:g}_{:g}".format(a[0], a[1], q1, q2)


todo = [(a, q1, q2) for a, q1, q2 in grid if not os.path.exists(os.path.join(points_dir, "bunch_{}.json".format(point_tag(a, q1, q2))))]
stamp("{} of {} points still to run".format(len(todo), len(grid)))
running = []
failed = []
queue = list(todo)
while queue or running:
    while queue and len(running) < args.parallel:
        a, q1, q2 = queue.pop(0)
        tag = point_tag(a, q1, q2)
        cmd = [PY, "-u", os.path.join(SCRIPTS, "BunchTrackSC.py"), "--tag", args.state_tag, "--reload-dir", args.run_dir,
               "--step-dir", args.step_dir, "--out-dir", points_dir, "--out-tag", tag,
               "--superpose", str(q1), str(a[0]), str(q2), str(a[1]), "--vscale", str(args.vscale), "--basis-dir", basis_dir,
               "--n", str(args.n), "--h", str(args.h), "--phi", str(args.phi), "--particles", args.particles,
               "--bfield", args.bfield, "--current-ma", str(args.current_ma), "--record", "200"]
        stamp("start {}".format(tag))
        running.append(launch(tag, cmd, os.path.join(points_dir, "log_{}.txt".format(tag))))
    for item in list(running):
        name, p, fh = item
        if p.poll() is not None:
            fh.close()
            running.remove(item)
            if p.returncode != 0:
                failed.append(name)
                stamp("FAILED: {} (exit {})".format(name, p.returncode))
            else:
                try:
                    with open(os.path.join(points_dir, "bunch_{}.json".format(name)), encoding="utf-8") as jf:
                        d = json.load(jf)
                    lb = d["losses_by_electrode"]
                    n = d["n_particles"]
                    stamp("done: {} -> {:.1f} % (spiral {:.1f}, housing {:.1f} %)".format(
                        name, 100 * d["transmission"], 100.0 * (lb.get("SI_Anode", 0) + lb.get("SI_Cathode", 0)) / n,
                        100.0 * lb.get("Housing_exit", 0) / n))
                except Exception as exc:  # noqa: BLE001
                    stamp("done: {} (no summary: {})".format(name, exc))
    if running:
        time.sleep(20)

# ------------------------------------------------------------------ 3. collect
rows = []
for a, q1, q2 in grid:
    tag = point_tag(a, q1, q2)
    fn = os.path.join(points_dir, "bunch_{}.json".format(tag))
    if not os.path.exists(fn):
        continue
    with open(fn, encoding="utf-8") as fh:
        d = json.load(fh)
    lb, n = d["losses_by_electrode"], d["n_particles"]
    quads = sum(v for k, v in lb.items() if k in ("D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7"))
    spiral = lb.get("SI_Anode", 0) + lb.get("SI_Cathode", 0)
    rows.append({"tag": tag, "alpha1": a[0], "alpha2": a[1], "q1": q1, "q2": q2, "n": n,
                 "transmission": d["transmission"], "lost_spiral": spiral / n, "lost_anode": lb.get("SI_Anode", 0) / n,
                 "lost_cathode": lb.get("SI_Cathode", 0) / n, "lost_quads": quads / n,
                 "lost_housing_exit": lb.get("Housing_exit", 0) / n,
                 "lost_apertures": (n - d["n_transmitted"] - spiral - quads) / n,
                 "z_exit_mean_mm": 1e3 * d["z_exit_mean_m"] if d["z_exit_mean_m"] is not None else None,
                 "wall_s": d["wall_time_s"], "sc": d.get("sc", {})})
with open(os.path.join(out, "scan_sc.json"), "w") as fh:
    json.dump({"args": vars(args), "v_ref": v_ref, "basis_dir": basis_dir, "rows": rows, "failed": failed}, fh, indent=2)
if not rows:
    sys.exit("no results")

best_T = max(r["transmission"] for r in rows)
cands = [r for r in rows if r["transmission"] >= best_T - 0.01]
best = min(cands, key=lambda r: r["lost_spiral"])
lines = ["# SC-aware quad retune: {}".format(args.name), "",
         "{} particles per point, PyAMG space charge {} mA at {:.1f} mm cells, vacuum field superposed from the basis solves of `{}` "
         "(spiral {:.1f} V x {:.4f}), beam angle {:.0f} deg. Transmission through the housing exit opening.".format(
             args.n, args.current_ma, 1e3 * args.h, args.run_dir, v_ref, args.vscale, args.phi), "",
         "| alpha1/alpha2 [deg] | q1 [V] | q2 [V] | transmission [%] | spiral | anode | quads | apertures | housing exit | exit z [mm] |",
         "|---|---|---|---|---|---|---|---|---|---|"]
for r in sorted(rows, key=lambda r: (r["alpha1"], r["alpha2"], r["q1"], r["q2"])):
    mark = "**" if r is best else ""
    lines.append("| {:g}/{:g} | {:+.0f} | {:+.0f} | {}{:.1f}{} | {:.1f} | {:.1f} | {:.1f} | {:.1f} | {:.1f} | {} |".format(
        r["alpha1"], r["alpha2"], r["q1"], r["q2"], mark, 100 * r["transmission"], mark, 100 * r["lost_spiral"], 100 * r["lost_anode"],
        100 * r["lost_quads"], 100 * r["lost_apertures"], 100 * r["lost_housing_exit"],
        "{:+.2f}".format(r["z_exit_mean_mm"]) if r["z_exit_mean_mm"] is not None else "-"))
lines += ["", "Best (transmission first, lowest spiral loss within 1 point): alpha {:g}/{:g}, q1 {:+.0f} V, q2 {:+.0f} V: {:.1f} %, spiral {:.1f} %.".format(
    best["alpha1"], best["alpha2"], best["q1"], best["q2"], 100 * best["transmission"], 100 * best["lost_spiral"])]
if failed:
    lines += ["", "Failed points: {}".format(", ".join(failed))]
stamp("best: alpha {:g}/{:g}, q1 {:+.0f}, q2 {:+.0f}: {:.1f} % (spiral {:.1f})".format(
    best["alpha1"], best["alpha2"], best["q1"], best["q2"], 100 * best["transmission"], 100 * best["lost_spiral"]))

# heat maps per rotation
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(alphas), figsize=(4.2 * len(alphas), 4), squeeze=False)
    vals = [100 * r["transmission"] for r in rows]
    for ax, a in zip(axes[0], alphas):
        M = np.full((len(args.q2), len(args.q1)), np.nan)
        for r in rows:
            if (r["alpha1"], r["alpha2"]) == a:
                M[args.q2.index(r["q2"]), args.q1.index(r["q1"])] = 100 * r["transmission"]
        im = ax.imshow(M, origin="lower", vmin=min(vals), vmax=max(vals), cmap="viridis", aspect="auto")
        for i in range(len(args.q2)):
            for j in range(len(args.q1)):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, "{:.1f}".format(M[i, j]), ha="center", va="center", color="w", fontsize=9)
        ax.set_xticks(range(len(args.q1)))
        ax.set_xticklabels(["{:+.0f}".format(q) for q in args.q1])
        ax.set_yticks(range(len(args.q2)))
        ax.set_yticklabels(["{:+.0f}".format(q) for q in args.q2])
        ax.set_xlabel("q1 [V]")
        ax.set_ylabel("q2 [V]")
        ax.set_title("quads rotated {:g}/{:g} deg".format(*a))
    fig.colorbar(im, ax=axes[0].tolist(), label="transmission [%], {} mA".format(args.current_ma))
    fig.suptitle("{}: SC-aware quad retune ({} particles, {:.0f} mm cells)".format(args.name, args.n, 1e3 * args.h))
    fig.savefig(os.path.join(out, "scan_sc.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
except Exception as exc:  # noqa: BLE001
    stamp("plot failed: {}".format(exc))

# ------------------------------------------------------------------ 4. confirmation
confirm = {}
if not args.skip_confirm:
    cdir = os.path.join(out, "confirm")
    os.makedirs(cdir, exist_ok=True)
    spiral_v = v_ref * args.vscale
    if not os.path.exists(os.path.join(cdir, "ef_itp_sc_best.pickle")):
        stamp("confirmation: direct solve at alpha {:g}/{:g}, q {:+.0f}/{:+.0f}, spiral {:.1f} V".format(
            best["alpha1"], best["alpha2"], best["q1"], best["q2"], spiral_v))
        cmd = [PY, "-u", os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "sc_best", "--quad-voltages", str(best["q1"]), str(best["q2"]),
               "--rotate-quads", str(best["alpha1"]), str(best["alpha2"]), "--spiral-voltage", str(spiral_v),
               "--step-dir", args.step_dir, "--voltages", os.path.join(args.run_dir, "voltages.csv"),
               "--state", os.path.join(args.run_dir, "state.pickle"), "--out-dir", cdir, "--res", str(args.res), "--bfield", args.bfield]
        if wait_all([launch("direct", cmd, os.path.join(cdir, "log_direct.txt"))]):
            sys.exit("direct solve failed")
    procs = []
    if not os.path.exists(os.path.join(cdir, "bunch_e_sc_best.json")):
        stamp("confirmation: {} particles, {} mA at {:.1f} mm cells (and the vacuum run in parallel)".format(args.confirm_n, args.current_ma, 1e3 * args.confirm_h))
        procs.append(launch("e_sc_best", [PY, "-u", os.path.join(SCRIPTS, "BunchTrackSC.py"), "--tag", "sc_best", "--reload-dir", cdir,
                                          "--step-dir", args.step_dir, "--out-dir", cdir, "--out-tag", "e_sc_best", "--n", str(args.confirm_n),
                                          "--h", str(args.confirm_h), "--phi", str(args.phi), "--particles", args.particles, "--bfield", args.bfield,
                                          "--current-ma", str(args.current_ma), "--save-openpmd", os.path.join(cdir, "handoff_sc_best.h5"),
                                          "--save-mode", "both"], os.path.join(cdir, "log_e_sc_best.txt")))
    if not os.path.exists(os.path.join(cdir, "bunch_e_nosc_best.json")):
        procs.append(launch("e_nosc_best", [PY, "-u", os.path.join(SCRIPTS, "BunchTrack.py"), "--tag", "sc_best", "--reload-dir", cdir,
                                            "--step-dir", args.step_dir, "--out-tag", "e_nosc_best", "--n", str(args.confirm_n), "--phi", str(args.phi),
                                            "--particles", args.particles, "--bfield", args.bfield, "--current-ma", str(args.current_ma),
                                            "--save-openpmd", os.path.join(cdir, "handoff_nosc_best.h5"), "--save-mode", "both"],
                            os.path.join(cdir, "log_e_nosc_best.txt")))
    wait_all(procs)
    for key in ("e_sc_best", "e_nosc_best"):
        npz = os.path.join(cdir, "bunch_{}.npz".format(key))
        if os.path.exists(npz):
            js = os.path.join(cdir, "exit_metrics_{}.json".format(key))
            subprocess.call([PY, os.path.join(SCRIPTS, "exit_metrics.py"), npz, os.path.join(cdir, "si_state_sc_best.pickle"), js],
                            cwd=SCRIPTS, stdout=open(os.path.join(cdir, "exit_metrics_{}.txt".format(key)), "w"), stderr=subprocess.STDOUT)
            if os.path.exists(js):
                with open(js) as fh:
                    confirm[key] = json.load(fh)
            subprocess.call([PY, os.path.join(SCRIPTS, "plot_geometry_trajectories.py"), npz, "--step-dir", args.step_dir,
                             "--out", os.path.join(cdir, "geometry_trajectories_{}.png".format(key)),
                             "--side", os.path.join(cdir, "geometry_side_{}.png".format(key))],
                            cwd=SCRIPTS, stdout=open(os.path.join(cdir, "log_plot_{}.txt".format(key)), "w"), stderr=subprocess.STDOUT)
    ref = None
    if args.reference and os.path.exists(args.reference):
        with open(args.reference) as fh:
            ref = json.load(fh)
    lines += ["", "## Confirmation at the best point ({} particles, direct BEM solve, spiral {:.1f} V)".format(args.confirm_n, spiral_v), "",
              "| run | transmission [%] | exit plane [%] | housing exit [%] | spiral [%] | asymptotic z [mm] | asymptotic vert. angle [deg] |",
              "|---|---|---|---|---|---|---|"]
    for key, label in (("e_nosc_best", "no space charge"), ("e_sc_best", "{} mA, {:.1f} mm cells".format(args.current_ma, 1e3 * args.confirm_h))):
        m = confirm.get(key)
        if m:
            lines.append("| {} | **{:.1f}** | {:.1f} | {:.2f} | {:.1f} | {:+.2f} +- {:.2f} | {:+.2f} +- {:.2f} |".format(
                label, 100 * m["transmission_through_housing"], 100 * m["transmission_exit_plane"], 100 * m["lost_housing_after_exit"],
                100 * m.get("lost_spiral", float("nan")), m["z_asym_mean_mm"], m["z_asym_rms_mm"], m["vert_angle_asym_mean_deg"], m["vert_angle_asym_rms_deg"]))
    if ref:
        lines.append("| chain's own SC run (reference, 2 mm cells) | {:.1f} | {:.1f} | {:.2f} | {:.1f} | {:+.2f} +- {:.2f} | {:+.2f} +- {:.2f} |".format(
            100 * ref["transmission_through_housing"], 100 * ref["transmission_exit_plane"], 100 * ref["lost_housing_after_exit"],
            100 * ref.get("lost_spiral", float("nan")), ref["z_asym_mean_mm"], ref["z_asym_rms_mm"], ref["vert_angle_asym_mean_deg"], ref["vert_angle_asym_rms_deg"]))
    lines += ["", "Resolution note (final2 ladder): 2 mm cells read +1.9, 1.5 mm +1.3 and 1 mm +0.8 points above the converged value.",
              "Hand-off files: `confirm/handoff_sc_best.h5` (+ `_lab6d.h5`) and `confirm/handoff_nosc_best.h5`.",
              "", "![geometry](confirm/geometry_trajectories_e_sc_best.png)", "", "![side](confirm/geometry_side_e_sc_best.png)"]

with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
with open(os.path.join(out, "summary.json"), "w") as fh:
    json.dump({"args": vars(args), "v_ref": v_ref, "best": best, "confirm": confirm, "failed": failed, "wall_s": time.time() - T0}, fh, indent=2)
stamp("done in {:.1f} min; report {}".format((time.time() - T0) / 60, os.path.join(out, "report.md")))
