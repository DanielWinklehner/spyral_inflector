"""Re-optimization of the spiral inflector's internal focusing on the Bevatech core beam.

Grid over the electrode-shape knobs sigma (V-depth), anglingAng (internal angling) and
gammaAng (gamma wedge) on top of a base geometry (pg5L quads/apertures by default). Every
point is a GeometryPoint run: build, design-particle optimization of spiral voltage and dz,
STEP export, basis fields with the rotated quads, quad-voltage retune by superposition, and a
bunch; the stages are

  1. scan    : all points at --res (5 mm, ranking only), --parallel at a time
  2. refine  : the --refine-top best points again at --refine-res with a bigger bunch
  3. confirm : the best refined point with a direct BEM solve and the full core beam with
               (--confirm-h cells) and without space charge, exit metrics, figures, hand-off files

Finished points are reused on a relaunch.

    python FocusScan.py --name focus1 --sigma 0.0012 0.0022 0.0032 --angling 7 11 15 --gamma 2 5 8
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
parser.add_argument("--out-root", default=os.path.join(DECK, "Results", "focus"))
parser.add_argument("--steps-root", default=os.path.join(DECK, "Geometry", "focus"))
parser.add_argument("--sigma", type=float, nargs="+", default=[0.0012, 0.0022, 0.0032], help="V-depth values [m]")
parser.add_argument("--angling", type=float, nargs="+", default=[7.0, 11.0, 15.0], help="anglingAng values [deg]")
parser.add_argument("--gamma", type=float, nargs="+", default=[2.0, 5.0, 8.0], help="gammaAng values [deg]")
parser.add_argument("--base", nargs="*", default=["--quad-bore", "0.018", "--aper-hole", "0.0175", "--slot-width", "0.019",
                                                   "--exit-opening", "0.023", "0.040", "--gap", "0.019", "--plate-gap", "0.005",
                                                   "--quad-z1", "-0.264", "--quad-len", "0.055", "--quad-len2", "0.060", "--shared-plates"],
                    help="base geometry knobs passed to GeometryPoint.py (default: pg5L)")
parser.add_argument("--phi", type=float, default=90.0)
parser.add_argument("--rotate-quads", type=float, nargs=2, default=[16.0, 24.0])
parser.add_argument("--q1", type=float, nargs=3, default=[5075, 8075, 4], metavar=("MIN", "MAX", "N"))
parser.add_argument("--q2", type=float, nargs=3, default=[-9950, -6950, 4], metavar=("MIN", "MAX", "N"))
parser.add_argument("--n-scan", type=int, default=1500)
parser.add_argument("--n-bunch", type=int, default=5000)
parser.add_argument("--res", type=float, default=0.005)
parser.add_argument("--parallel", type=int, default=3)
parser.add_argument("--refine-top", type=int, default=3)
parser.add_argument("--refine-res", type=float, default=0.0025)
parser.add_argument("--refine-n-bunch", type=int, default=10000)
parser.add_argument("--confirm-n", type=int, default=43969)
parser.add_argument("--confirm-h", type=float, default=0.0015)
parser.add_argument("--current-ma", type=float, default=8.0)
parser.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit_core_as_txt.txt"))
parser.add_argument("--bfield", default=os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle"))
parser.add_argument("--reference", default=None, help="exit_metrics json of the base geometry's SC run, for the report")
parser.add_argument("--skip-refine", action="store_true")
parser.add_argument("--skip-confirm", action="store_true")
args = parser.parse_args()

PY = sys.executable
out = os.path.join(args.out_root, args.name)
os.makedirs(out, exist_ok=True)
T0 = time.time()


def stamp(msg):
    print("[{}] {}".format(time.strftime("%H:%M:%S"), msg), flush=True)


def launch(name, cmd, log_path):
    fh = open(log_path, "w", encoding="utf-8")
    return name, subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS), fh


def wait_all(procs, poll=20):
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


def point_name(sigma, ang, gam):
    return "s{:g}_a{:g}_g{:g}".format(1e3 * sigma, ang, gam)


def point_cmd(name, sigma, ang, gam, out_root, steps_root, res, n_bunch):
    return [PY, "-u", os.path.join(SCRIPTS, "GeometryPoint.py"), "--name", name, "--out-root", out_root, "--steps-root", steps_root,
            "--particles", args.particles, "--bfield", args.bfield, "--phi", str(args.phi),
            "--rotate-quads", str(args.rotate_quads[0]), str(args.rotate_quads[1]),
            "--q1", str(args.q1[0]), str(args.q1[1]), str(int(args.q1[2])), "--q2", str(args.q2[0]), str(args.q2[1]), str(int(args.q2[2])),
            "--n-scan", str(args.n_scan), "--n-bunch", str(n_bunch), "--res", str(res),
            "--sigma", str(sigma), "--angling", str(ang), "--gamma", str(gam)] + args.base


def run_points(points, out_root, steps_root, res, n_bunch, label):
    """Run GeometryPoint for every (sigma, angling, gamma) without a finished summary; pool of --parallel."""
    os.makedirs(out_root, exist_ok=True)
    queue = [(point_name(*pt), pt) for pt in points]
    queue = [(nm, pt) for nm, pt in queue if not _finished(os.path.join(out_root, nm, "summary.json"))]
    stamp("{}: {} of {} points to run".format(label, len(queue), len(points)))
    running, failed = [], []
    while queue or running:
        while queue and len(running) < args.parallel:
            nm, pt = queue.pop(0)
            stamp("start {}".format(nm))
            running.append(launch(nm, point_cmd(nm, *pt, out_root=out_root, steps_root=steps_root, res=res, n_bunch=n_bunch),
                                  os.path.join(out_root, "log_{}.txt".format(nm))))
        for item in list(running):
            nm, p, fh = item
            if p.poll() is not None:
                fh.close()
                running.remove(item)
                if p.returncode != 0:
                    failed.append(nm)
                    stamp("FAILED: {} (exit {})".format(nm, p.returncode))
                else:
                    r = _row(out_root, nm)
                    stamp("done: {} -> {:.1f} % (spiral {:.1f}, q {:+.0f}/{:+.0f}, V {:.0f}, dz {:+.2f} mm)".format(
                        nm, 100 * r["transmission"], 100 * r["lost_spiral"], r["q1"], r["q2"], r["spiral_V"], r["dz_mm"]) if r else "done: {} (no summary)".format(nm))
        if running:
            time.sleep(20)
    return failed


def _finished(summary_fn):
    if not os.path.exists(summary_fn):
        return False
    with open(summary_fn, encoding="utf-8") as fh:
        return "bunch" in json.load(fh)


def _row(out_root, nm):
    fn = os.path.join(out_root, nm, "summary.json")
    if not _finished(fn):
        return None
    with open(fn, encoding="utf-8") as fh:
        s = json.load(fh)
    b, k, o = s["bunch"], s["knobs"], s["optimizer"]
    return {"name": nm, "sigma_mm": 1e3 * k["sigma"], "angling": k["angling"], "gamma": k["gamma"],
            "transmission": b.get("transmission_through_housing", b["transmission"]), "transmission_exit_plane": b.get("transmission_exit_plane", b["transmission"]),
            "lost_spiral": b["lost_spiral"], "lost_quads": b["lost_quads"], "lost_apertures": b["lost_apertures"],
            "lost_housing_exit": b.get("lost_housing_after_exit", 0.0), "z_asym_mm": b.get("z_asym_mean_mm"), "angle_asym_deg": b.get("vert_angle_asym_mean_deg"),
            "q1": s["inner_scan"]["best"]["q1"], "q2": s["inner_scan"]["best"]["q2"], "inner_best_T": s["inner_scan"]["best"]["transmission"],
            "spiral_V": o["voltage_V"], "dz_mm": o["dz_mm"], "clearance_mm": o["min_clearance_mm"], "n_bunch": b["n"], "wall_s": s.get("wall_s")}


def rank(rows, tie=0.01):
    """Transmission first; within `tie` of the best, the lowest spiral loss."""
    if not rows:
        return []
    top = max(r["transmission"] for r in rows)
    cands = sorted([r for r in rows if r["transmission"] >= top - tie], key=lambda r: (r["lost_spiral"], -r["transmission"]))
    rest = sorted([r for r in rows if r["transmission"] < top - tie], key=lambda r: -r["transmission"])
    return cands + rest


def table(rows, best_names=()):
    lines = ["| point | sigma [mm] | angling [deg] | gamma [deg] | transmission [%] | spiral | quads | apertures | housing exit | q1/q2 [V] | spiral V | dz [mm] | asym z [mm] | asym angle [deg] |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: (r["gamma"], r["angling"], r["sigma_mm"])):
        m = "**" if r["name"] in best_names else ""
        lines.append("| {} | {:.1f} | {:g} | {:g} | {}{:.1f}{} | {:.1f} | {:.1f} | {:.1f} | {:.1f} | {:+.0f}/{:+.0f} | {:.0f} | {:+.2f} | {} | {} |".format(
            r["name"], r["sigma_mm"], r["angling"], r["gamma"], m, 100 * r["transmission"], m, 100 * r["lost_spiral"], 100 * r["lost_quads"],
            100 * r["lost_apertures"], 100 * r["lost_housing_exit"], r["q1"], r["q2"], r["spiral_V"], r["dz_mm"],
            "{:+.2f}".format(r["z_asym_mm"]) if r["z_asym_mm"] is not None else "-",
            "{:+.2f}".format(r["angle_asym_deg"]) if r["angle_asym_deg"] is not None else "-"))
    return lines


def heatmaps(rows, png, title):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # noqa: BLE001
        return
    gams = sorted(set(r["gamma"] for r in rows))
    sigs = sorted(set(r["sigma_mm"] for r in rows))
    angs = sorted(set(r["angling"] for r in rows))
    vals = [100 * r["transmission"] for r in rows]
    fig, axes = plt.subplots(1, len(gams), figsize=(4.4 * len(gams), 4.2), squeeze=False)
    for ax, g in zip(axes[0], gams):
        M = np.full((len(angs), len(sigs)), np.nan)
        for r in rows:
            if r["gamma"] == g:
                M[angs.index(r["angling"]), sigs.index(r["sigma_mm"])] = 100 * r["transmission"]
        im = ax.imshow(M, origin="lower", vmin=min(vals), vmax=max(vals), cmap="viridis", aspect="auto")
        for i in range(len(angs)):
            for j in range(len(sigs)):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, "{:.1f}".format(M[i, j]), ha="center", va="center", color="w", fontsize=9)
        ax.set_xticks(range(len(sigs)))
        ax.set_xticklabels(["{:.1f}".format(s) for s in sigs])
        ax.set_yticks(range(len(angs)))
        ax.set_yticklabels(["{:g}".format(a) for a in angs])
        ax.set_xlabel("sigma (V-depth) [mm]")
        ax.set_ylabel("angling [deg]")
        ax.set_title("gamma wedge {:g} deg".format(g))
    fig.colorbar(im, ax=axes[0].tolist(), label="transmission through the housing [%]")
    fig.suptitle(title)
    fig.savefig(png, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ 1. scan
points = list(itertools.product(args.sigma, args.angling, args.gamma))
scan_root = os.path.join(out, "scan")
scan_steps = os.path.join(args.steps_root, args.name, "scan")
stamp("focus scan '{}': {} sigma x {} angling x {} gamma = {} points at {:.1f} mm, quads {}/{} deg, base {}".format(
    args.name, len(args.sigma), len(args.angling), len(args.gamma), len(points), 1e3 * args.res, args.rotate_quads[0], args.rotate_quads[1], " ".join(args.base)))
failed = run_points(points, scan_root, scan_steps, args.res, args.n_bunch, "scan")
rows = [r for r in (_row(scan_root, point_name(*pt)) for pt in points) if r]
ranked = rank(rows)
lines = ["# Internal-focusing re-optimization: {}".format(args.name), "",
         "Base geometry `{}`; beam angle {:g} deg, quads rotated {:g}/{:g} deg with a {}x{} voltage scan per point; "
         "{} particles per bunch; fields at {:.1f} mm (ranking only: 5 mm reads ~18 points below 2.5 mm).".format(
             " ".join(args.base), args.phi, args.rotate_quads[0], args.rotate_quads[1], int(args.q1[2]), int(args.q2[2]), args.n_bunch, 1e3 * args.res), "",
         "## 1. Scan ({} points)".format(len(rows)), ""] + table(rows, [r["name"] for r in ranked[:args.refine_top]]) + \
        ["", "Ranking: transmission through the housing first, lowest spiral loss within 1 point. Top {}: {}.".format(
            args.refine_top, ", ".join("{} ({:.1f} %)".format(r["name"], 100 * r["transmission"]) for r in ranked[:args.refine_top]))]
if failed:
    lines.append("Failed points: {}".format(", ".join(failed)))
with open(os.path.join(out, "scan.json"), "w") as fh:
    json.dump({"args": vars(args), "rows": rows, "ranked": [r["name"] for r in ranked], "failed": failed}, fh, indent=2)
heatmaps(rows, os.path.join(out, "scan.png"), "{}: transmission at {:.0f} mm fields, {} particles".format(args.name, 1e3 * args.res, args.n_bunch))
with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
stamp("scan done: best {} ({:.1f} %, spiral {:.1f} %)".format(ranked[0]["name"], 100 * ranked[0]["transmission"], 100 * ranked[0]["lost_spiral"]) if ranked else "scan done: no results")

# ------------------------------------------------------------------ 2. refine
best = None
if ranked and not args.skip_refine:
    top = ranked[:args.refine_top]
    refine_root = os.path.join(out, "refine")
    refine_steps = os.path.join(args.steps_root, args.name, "refine")
    pts = [(1e-3 * r["sigma_mm"], r["angling"], r["gamma"]) for r in top]
    failed_r = run_points(pts, refine_root, refine_steps, args.refine_res, args.refine_n_bunch, "refine")
    rrows = [r for r in (_row(refine_root, point_name(*pt)) for pt in pts) if r]
    rranked = rank(rrows)
    best = rranked[0] if rranked else None
    lines += ["", "## 2. Refinement at {:.2f} mm, {} particles".format(1e3 * args.refine_res, args.refine_n_bunch), ""] + table(rrows, [best["name"]] if best else [])
    if failed_r:
        lines.append("Failed points: {}".format(", ".join(failed_r)))
    with open(os.path.join(out, "refine.json"), "w") as fh:
        json.dump({"rows": rrows, "ranked": [r["name"] for r in rranked], "failed": failed_r}, fh, indent=2)
    with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    if best:
        stamp("refine done: best {} ({:.1f} %, spiral {:.1f} %) at q {:+.0f}/{:+.0f}, {:.0f} V, dz {:+.2f} mm".format(
            best["name"], 100 * best["transmission"], 100 * best["lost_spiral"], best["q1"], best["q2"], best["spiral_V"], best["dz_mm"]))
elif ranked:
    best = ranked[0]
    refine_root, refine_steps = scan_root, scan_steps

# ------------------------------------------------------------------ 3. confirm
confirm = {}
if best and not args.skip_confirm:
    cdir = os.path.join(out, "confirm")
    os.makedirs(cdir, exist_ok=True)
    pdir = os.path.join(refine_root, best["name"])
    sdir = os.path.join(refine_steps, "{}_steps".format(best["name"]))
    if not os.path.exists(os.path.join(cdir, "ef_itp_best.pickle")):
        stamp("confirmation: direct solve of {} at q {:+.0f}/{:+.0f}, rotated {:g}/{:g}".format(best["name"], best["q1"], best["q2"], *args.rotate_quads))
        cmd = [PY, "-u", os.path.join(SCRIPTS, "TrackFromStep.py"), "--tag", "best", "--quad-voltages", str(best["q1"]), str(best["q2"]),
               "--rotate-quads", str(args.rotate_quads[0]), str(args.rotate_quads[1]), "--step-dir", sdir,
               "--voltages", os.path.join(pdir, "voltages.csv"), "--state", os.path.join(pdir, "state.pickle"), "--out-dir", cdir,
               "--res", str(args.refine_res), "--bfield", args.bfield]
        if wait_all([launch("direct", cmd, os.path.join(cdir, "log_direct.txt"))]):
            sys.exit("direct solve failed")
    procs = []
    if not os.path.exists(os.path.join(cdir, "bunch_e_sc.json")):
        stamp("confirmation: {} particles, {} mA at {:.1f} mm cells, and the vacuum run in parallel".format(args.confirm_n, args.current_ma, 1e3 * args.confirm_h))
        procs.append(launch("e_sc", [PY, "-u", os.path.join(SCRIPTS, "BunchTrackSC.py"), "--tag", "best", "--reload-dir", cdir, "--step-dir", sdir,
                                     "--out-dir", cdir, "--out-tag", "e_sc", "--n", str(args.confirm_n), "--h", str(args.confirm_h), "--phi", str(args.phi),
                                     "--particles", args.particles, "--bfield", args.bfield, "--current-ma", str(args.current_ma),
                                     "--save-openpmd", os.path.join(cdir, "handoff_sc.h5"), "--save-mode", "both"], os.path.join(cdir, "log_e_sc.txt")))
    if not os.path.exists(os.path.join(cdir, "bunch_e_nosc.json")):
        procs.append(launch("e_nosc", [PY, "-u", os.path.join(SCRIPTS, "BunchTrack.py"), "--tag", "best", "--reload-dir", cdir, "--step-dir", sdir,
                                       "--out-tag", "e_nosc", "--n", str(args.confirm_n), "--phi", str(args.phi), "--particles", args.particles,
                                       "--bfield", args.bfield, "--current-ma", str(args.current_ma),
                                       "--save-openpmd", os.path.join(cdir, "handoff_nosc.h5"), "--save-mode", "both"], os.path.join(cdir, "log_e_nosc.txt")))
    wait_all(procs)
    for key in ("e_nosc", "e_sc"):
        npz = os.path.join(cdir, "bunch_{}.npz".format(key))
        if os.path.exists(npz):
            js = os.path.join(cdir, "exit_metrics_{}.json".format(key))
            with open(os.path.join(cdir, "exit_metrics_{}.txt".format(key)), "w") as fh:
                subprocess.call([PY, os.path.join(SCRIPTS, "exit_metrics.py"), npz, os.path.join(cdir, "si_state_best.pickle"), js], cwd=SCRIPTS, stdout=fh, stderr=subprocess.STDOUT)
            if os.path.exists(js):
                with open(js) as fh:
                    confirm[key] = json.load(fh)
            with open(os.path.join(cdir, "log_plot_{}.txt".format(key)), "w") as fh:
                subprocess.call([PY, os.path.join(SCRIPTS, "plot_geometry_trajectories.py"), npz, "--step-dir", sdir,
                                 "--out", os.path.join(cdir, "geometry_trajectories_{}.png".format(key)),
                                 "--side", os.path.join(cdir, "geometry_side_{}.png".format(key))], cwd=SCRIPTS, stdout=fh, stderr=subprocess.STDOUT)
    ref = None
    if args.reference and os.path.exists(args.reference):
        with open(args.reference) as fh:
            ref = json.load(fh)
    lines += ["", "## 3. Confirmation of {} ({} particles, direct BEM solve, q {:+.0f}/{:+.0f} V, spiral {:.0f} V, dz {:+.2f} mm)".format(
        best["name"], args.confirm_n, best["q1"], best["q2"], best["spiral_V"], best["dz_mm"]), "",
        "| run | transmission [%] | exit plane [%] | housing exit [%] | spiral [%] | asymptotic z [mm] | asymptotic vert. angle [deg] |",
        "|---|---|---|---|---|---|---|"]
    for key, label in (("e_nosc", "no space charge"), ("e_sc", "{} mA, {:.1f} mm cells".format(args.current_ma, 1e3 * args.confirm_h))):
        m = confirm.get(key)
        if m and "z_asym_mean_mm" in m:
            lines.append("| {} | **{:.1f}** | {:.1f} | {:.2f} | {:.1f} | {:+.2f} +- {:.2f} | {:+.2f} +- {:.2f} |".format(
                label, 100 * m["transmission_through_housing"], 100 * m["transmission_exit_plane"], 100 * m["lost_housing_after_exit"],
                100 * m["lost_spiral"], m["z_asym_mean_mm"], m["z_asym_rms_mm"], m["vert_angle_asym_mean_deg"], m["vert_angle_asym_rms_deg"]))
    if ref:
        lines.append("| base geometry's own SC run (reference) | {:.1f} | {:.1f} | {:.2f} | {:.1f} | {:+.2f} +- {:.2f} | {:+.2f} +- {:.2f} |".format(
            100 * ref["transmission_through_housing"], 100 * ref["transmission_exit_plane"], 100 * ref["lost_housing_after_exit"],
            100 * ref.get("lost_spiral", float("nan")), ref["z_asym_mean_mm"], ref["z_asym_rms_mm"], ref["vert_angle_asym_mean_deg"], ref["vert_angle_asym_rms_deg"]))
    lines += ["", "Resolution note (final2 ladder): 2 mm cells read +1.9, 1.5 mm +1.3 and 1 mm +0.8 points above the converged value.",
              "The dz of each point is the design particle's; the bunch centroid was not re-levelled (run the fringe optimizer on the winner for that).",
              "Hand-off files: `confirm/handoff_sc.h5` (+ `_lab6d.h5`), `confirm/handoff_nosc.h5`.", "",
              "![geometry](confirm/geometry_trajectories_e_sc.png)", "", "![side](confirm/geometry_side_e_sc.png)"]

with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
with open(os.path.join(out, "summary.json"), "w") as fh:
    json.dump({"args": vars(args), "scan_best": ranked[0] if ranked else None, "best": best, "confirm": confirm, "wall_s": time.time() - T0}, fh, indent=2)
stamp("done in {:.1f} min; report {}".format((time.time() - T0) / 60, os.path.join(out, "report.md")))
