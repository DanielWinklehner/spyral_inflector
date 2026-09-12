"""Electrode-shape scan for VERTICAL focusing, in the Baseline frame on top of a finished
final-push run.

Every point of the grid sigma (V-depth) x anglingAng x gammaAng is built in the base run's
rotated field (the field a system rotated by R sees), its design orbit optimized with the
exit truncation and dz as knobs, its own rotation R solved from the exit-plate rule, three
basis fields of the rigidly rotated system solved, and the quads retuned by superposition
on a bunch -- ranked by the rms vertical angle vz / v_longitudinal of the transmitted bunch
at the asymptotic state (55 mm past the exit), among the quad settings within --rank-tol of
the best transmission. Points are then ranked the same way. The size envelope over
0..150 mm past the exit is recorded alongside.

    python VFocusScan.py --name vfocus1 --sigma 0.0004 0.0012 0.0022 --angling 7 11 15 --gamma 5 8 11
    python VFocusScan.py --name vfocus1 --final       # ... and run RunFinalPush for the winner

About 30 min per point (design-orbit solve 12, basis solves 8, 4x4 quad grid 10). Finished
points are reused on a relaunch. Outputs: Results/vfocus/<name>/points/<point>/,
report.md, scan.png, summary.json, winner_knobs.json.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
PY = sys.executable
sys.path.insert(0, SCRIPTS)

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--name", required=True)
p.add_argument("--base-run", default=os.path.join(DECK, "Results", "final", "final1"),
               help="finished RunFinalPush run: its knobs, rotated field and R are the starting point")
p.add_argument("--sigma", type=float, nargs="+", default=[0.0004, 0.0012, 0.0022], help="V-depth [m]")
p.add_argument("--angling", type=float, nargs="+", default=[7.0, 11.0, 15.0], help="anglingAng [deg]")
p.add_argument("--gamma", type=float, nargs="+", default=[5.0, 8.0, 11.0], help="gammaAng [deg]")
p.add_argument("--fix-entrance", type=float, default=0.34, help="entrance truncation [deg]")
p.add_argument("--fix-exit", type=float, default=None, help="exit truncation [deg]; default: free knob")
p.add_argument("--maxiter", type=int, default=4)
p.add_argument("--geo-res", type=float, default=0.005)
p.add_argument("--res", type=float, default=0.005, help="basis-field resolution [m]")
p.add_argument("--h", type=float, default=0.005)
p.add_argument("--q1", type=float, nargs=3, default=[4575, 7575, 4], metavar=("MIN", "MAX", "N"))
p.add_argument("--q2", type=float, nargs=3, default=[-8450, -5450, 4], metavar=("MIN", "MAX", "N"))
p.add_argument("--n-scan", type=int, default=1000)
p.add_argument("--rank-tol", type=float, default=0.03, help="transmission tolerance of the vertical-angle ranking (fraction)")
p.add_argument("--quads", type=float, nargs=2, default=[16.0, 24.0], help="quad rotations on top of R [deg]")
p.add_argument("--gap-azimuth", type=float, default=34.0)
p.add_argument("--half-gap", type=float, default=0.005)
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle"))
p.add_argument("--final", action="store_true", help="run RunFinalPush for the winner afterwards")
p.add_argument("--final-name", default=None)
args = p.parse_args()

OUT = os.path.join(DECK, "Results", "vfocus", args.name)
PTS = os.path.join(OUT, "points")
os.makedirs(PTS, exist_ok=True)
LOG = open(os.path.join(OUT, "log_scan.txt"), "a", encoding="utf-8")
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
        raise RuntimeError("{} failed (exit {}), see {}".format(os.path.basename(str(cmd[0])), r.returncode, log_path))


def jload(fn):
    with open(fn) as fh:
        return json.load(fh)


from spyral_inflector.tracking.deck import particle_rows  # noqa: E402
from spyral_inflector.tracking.geometry_point import build_geometry  # noqa: E402
from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402

base_knobs = jload(os.path.join(args.base_run, "geometry", "summary.json"))["knobs"]
ROT_BF = os.path.join(args.base_run, "bfield_rotated.pickle")
R_BASE = jload(os.path.join(args.base_run, "exit_plane.json"))["rotation_deg"]
ENERGY = float(particle_rows(args.particles, core=True)[:, 8].mean())
TRUNC = (args.fix_entrance, args.fix_exit)
stamp("VFOCUS SCAN '{}': {} sigma x {} angling x {} gamma = {} points; base {} (R = {:+.3f} deg), truncations {}".format(
    args.name, len(args.sigma), len(args.angling), len(args.gamma), len(args.sigma) * len(args.angling) * len(args.gamma),
    os.path.basename(args.base_run), R_BASE, TRUNC))
stamp("  design energy {:.5f} MeV, quad grid {}x{} at {} particles, ranked by the rms vertical angle within {:.0%} of the best transmission".format(
    ENERGY, int(args.q1[2]), int(args.q2[2]), args.n_scan, args.rank_tol))


def point_name(sig, ang, gam):
    return "s{:.1f}_a{:g}_g{:g}".format(1e3 * sig, ang, gam)


def run_point(sig, ang, gam):
    name = point_name(sig, ang, gam)
    pdir = os.path.join(PTS, name)
    done_fn = os.path.join(pdir, "point.json")
    if os.path.exists(done_fn):
        return jload(done_fn)
    os.makedirs(pdir, exist_ok=True)
    t0 = time.time()
    knobs = dict(base_knobs)
    knobs.update(sigma=sig, angling=ang, gamma=gam)
    GEO, STEPS, BASIS = (os.path.join(pdir, d) for d in ("geometry", "steps", "basis"))
    stamp("point {}: build + design orbit (exit truncation {})".format(name, "free" if TRUNC[1] is None else "fixed"))
    if not os.path.exists(os.path.join(GEO, "summary.json")):
        if os.path.isdir(GEO):
            shutil.rmtree(GEO)
        build_geometry(GEO, STEPS, ROT_BF, ENERGY, knobs=knobs, fix_truncations=TRUNC, maxiter=args.maxiter,
                       res=args.geo_res, h=args.h, log=stamp)
    opt = jload(os.path.join(GEO, "summary.json"))["optimizer"]
    spec = exit_plane_spec(STEPS, os.path.join(GEO, "state.pickle"), args.gap_azimuth, args.half_gap,
                           out_json=os.path.join(pdir, "exit_plane.json"), log=stamp)
    R = spec["rotation_deg"]
    stamp("  R = {:+.3f} deg (base {:+.3f}); V {:.0f} V, dz {:+.2f} mm, truncations {}, residuals angle {:+.3f} deg, z {:+.2f} mm".format(
        R, R_BASE, opt["voltage_V"], opt["dz_mm"], opt.get("truncations_deg"), opt["residual_final"]["angle_deg"],
        opt["residual_final"]["z_offset_mm"]))
    os.makedirs(BASIS, exist_ok=True)
    common = ["--step-dir", STEPS, "--voltages", os.path.join(GEO, "voltages.csv"), "--state", os.path.join(GEO, "state.pickle"),
              "--out-dir", BASIS, "--res", args.res, "--h", args.h, "--bfield", args.bfield, "--rotate-all", R,
              "--rotate-quads", args.quads[0], args.quads[1], "--energy-mev", ENERGY]
    for tag, extra in (("spiral", ["--quad-voltages", 0, 0]),
                       ("q1", ["--spiral-voltage", 0, "--quad-voltages", 3500, 0, "--no-test"]),
                       ("q2", ["--spiral-voltage", 0, "--quad-voltages", 0, 3500, "--no-test"])):
        if not os.path.exists(os.path.join(BASIS, "ef_itp_{}.pickle".format(tag))):
            stamp("  basis {}".format(tag))
            run(["-m", "spyral_inflector.tracking.bem_reload", "--tag", tag] + extra + common, os.path.join(pdir, "log_basis_{}.txt".format(tag)))
    if not os.path.exists(os.path.join(BASIS, "scan2_scan.json")):
        stamp("  quad grid {}x{}, ranked by the rms vertical angle".format(int(args.q1[2]), int(args.q2[2])))
        run([os.path.join(SCRIPTS, "BunchScan2.py"), "--reload-dir", BASIS, "--out-dir", BASIS, "--step-dir", STEPS,
             "--bfield", args.bfield, "--basis", "spiral", "q1", "q1", "q2", "q2", "--phi", 90.0 + R,
             "--particles", args.particles, "--n", args.n_scan, "--q1", args.q1[0], args.q1[1], int(args.q1[2]),
             "--q2", args.q2[0], args.q2[1], int(args.q2[2]), "--tag", "scan", "--rank", "vert", "--rank-tol", args.rank_tol],
            os.path.join(pdir, "log_scan.txt"))
    sc = jload(os.path.join(BASIS, "scan2_scan.json"))
    rows = [r for r in sc["results"] if r.get("vfom") is not None]
    t_best = max(r["transmission"] for r in rows)
    # the point's quad setting: the smallest rms vertical angle within --rank-tol of the best
    # transmission on its grid (re-picked from the stored grid, so the tolerance can change
    # without re-tracking)
    cands = [r for r in rows if r["transmission"] >= t_best - args.rank_tol]
    best = min(cands, key=lambda r: (r["vfom"], -r["transmission"]))
    t_choice = max(rows, key=lambda r: r["transmission"])
    point = {"name": name, "sigma_mm": 1e3 * sig, "angling": ang, "gamma": gam, "R_deg": R,
             "voltage_V": opt["voltage_V"], "dz_mm": opt["dz_mm"], "truncations_deg": opt.get("truncations_deg"),
             "residual_angle_deg": opt["residual_final"]["angle_deg"], "residual_z_mm": opt["residual_final"]["z_offset_mm"],
             "q1": best["q1"], "q2": best["q2"], "transmission": best["transmission"], "vfom": best["vfom"],
             "zp_asym_rms_mrad": best.get("zp_asym_rms_mrad"), "zp_asym_mean_mrad": best.get("zp_asym_mean_mrad"),
             "z_env_rms_mm": best.get("z_env_rms_mm"), "z_env_mean_mm": best.get("z_env_mean_mm"), "z_asym_rms_mm": best["z_asym_rms_mm"],
             "vert_angle_asym_rms_deg": best["vert_angle_asym_rms_deg"], "z_asym_mean_mm": best["z_asym_mean_mm"],
             "lost_spiral": best["lost_spiral"], "lost_apertures": best["lost_apertures"],
             "best_transmission_on_grid": t_best, "transmission_optimum": {"q1": t_choice["q1"], "q2": t_choice["q2"],
                                                                         "vfom": t_choice["vfom"]},
             "wall_min": (time.time() - t0) / 60.0}
    with open(done_fn, "w") as fh:
        json.dump(point, fh, indent=2)
    stamp("  done {}: T {:.1f} % (grid best {:.1f}), vertical angle {:.1f} mrad rms (envelope {:.1f} mm), q {:+.0f}/{:+.0f}, {:.0f} min".format(
        name, 100 * best["transmission"], 100 * t_best, best["vfom"] or float("nan"), best.get("z_env_mean_mm") or float("nan"),
        best["q1"], best["q2"], point["wall_min"]))
    return point


rows, failed = [], []
for gam in args.gamma:
    for ang in args.angling:
        for sig in args.sigma:
            try:
                rows.append(run_point(sig, ang, gam))
            except Exception as exc:
                stamp("  FAILED {}: {}".format(point_name(sig, ang, gam), exc))
                traceback.print_exc()
                failed.append(point_name(sig, ang, gam))

# ---------------------------------------------------------------- ranking, report, figure
rows = [r for r in rows if r.get("vfom") is not None]
top_t = max(r["transmission"] for r in rows) if rows else 0.0
cands = sorted([r for r in rows if r["transmission"] >= top_t - args.rank_tol], key=lambda r: (r["vfom"], -r["transmission"]))
rest = sorted([r for r in rows if r["transmission"] < top_t - args.rank_tol], key=lambda r: r["vfom"])
ranked = cands + rest
by_t = sorted(rows, key=lambda r: -r["transmission"])
base = next((r for r in rows if abs(r["sigma_mm"] - 1e3 * base_knobs["sigma"]) < 1e-6 and r["angling"] == base_knobs["angling"]
             and r["gamma"] == base_knobs["gamma"]), None)

L = ["# Vertical-focusing shape scan: {}".format(args.name), "",
     "Base run `{}` (knobs except sigma/angling/gamma, its rotated field for the design-orbit solves; R re-solved per point). "
     "Exit truncation {}, dz and voltage from the design-orbit optimizer ({} iterations). Basis fields at {:.0f} mm, "
     "quad grid {}x{} at {} particles. Ranking: the rms vertical angle vz/v_long of the transmitted bunch at the asymptotic "
     "state (55 mm past the exit) among the points within {:.0f} points of the best transmission; the size envelope "
     "(mean rms height over 0..150 mm past the exit) is listed alongside.".format(
         os.path.basename(args.base_run), "free" if TRUNC[1] is None else "fixed at {:.2f} deg".format(TRUNC[1]), args.maxiter,
         1e3 * args.res, int(args.q1[2]), int(args.q2[2]), args.n_scan, 100 * args.rank_tol), "",
     "| rank | point | sigma [mm] | angling | gamma | exit trunc [deg] | V [V] | dz [mm] | R [deg] | q1/q2 [V] | transmission | z' rms [mrad] | envelope [mm rms] | z rms 55 mm | spiral loss |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in enumerate(ranked):
    L.append("| {} | {} | {:.1f} | {:g} | {:g} | {} | {:.0f} | {:+.2f} | {:+.2f} | {:+.0f}/{:+.0f} | {:.1f} % | **{:.2f}** | {:.2f} | {:.2f} | {:.1f} % |".format(
        i + 1, r["name"] + (" (base)" if r is base else ""), r["sigma_mm"], r["angling"], r["gamma"],
        "{:.2f}".format(r["truncations_deg"][1]) if r.get("truncations_deg") else "-", r["voltage_V"], r["dz_mm"], r["R_deg"],
        r["q1"], r["q2"], 100 * r["transmission"], r["vfom"], r.get("z_env_mean_mm") or float("nan"), r["z_asym_rms_mm"], 100 * r["lost_spiral"]))
if ranked:
    w = ranked[0]
    L += ["", "**Winner: {}** -- vertical angle {:.1f} mrad rms at {:.1f} % transmission (best transmission on the grid: {} at {:.1f} %, "
          "{:.1f} mrad).".format(w["name"], w["vfom"], 100 * w["transmission"], by_t[0]["name"], 100 * by_t[0]["transmission"],
                                 by_t[0]["vfom"])]
    if base:
        L.append("Base point {}: {:.1f} mrad rms at {:.1f} %.".format(base["name"], base["vfom"], 100 * base["transmission"]))
if failed:
    L += ["", "Failed points: " + ", ".join(failed)]
L += ["", "![scan](scan.png)", "", "Wall time {:.0f} min.".format((time.time() - T0) / 60)]
with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
with open(os.path.join(OUT, "summary.json"), "w") as fh:
    json.dump({"args": vars(args), "rows": rows, "ranked": [r["name"] for r in ranked], "failed": failed,
               "winner": ranked[0] if ranked else None}, fh, indent=2)
if ranked:
    w = ranked[0]
    with open(os.path.join(OUT, "winner_knobs.json"), "w") as fh:
        json.dump({"sigma": 1e-3 * w["sigma_mm"], "angling": w["angling"], "gamma": w["gamma"]}, fh, indent=2)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sigs, angs, gams = sorted(set(r["sigma_mm"] for r in rows)), sorted(set(r["angling"] for r in rows)), sorted(set(r["gamma"] for r in rows))
    fig, axes = plt.subplots(2, len(gams), figsize=(4.2 * len(gams), 7.5), squeeze=False)
    for j, gam in enumerate(gams):
        for i, (key, label, cmap) in enumerate((("vfom", "rms vertical angle [mrad]", "viridis_r"), ("transmission", "transmission [%]", "viridis"))):
            M = np.full((len(angs), len(sigs)), np.nan)
            for r in rows:
                if r["gamma"] == gam:
                    M[angs.index(r["angling"]), sigs.index(r["sigma_mm"])] = 100 * r[key] if key == "transmission" else r[key]
            ax = axes[i, j]
            im = ax.imshow(M, origin="lower", cmap=cmap, aspect="auto")
            ax.set_xticks(range(len(sigs))); ax.set_xticklabels(["{:.1f}".format(s) for s in sigs])
            ax.set_yticks(range(len(angs))); ax.set_yticklabels(["{:g}".format(a) for a in angs])
            ax.set_xlabel("sigma [mm]"); ax.set_ylabel("angling [deg]")
            ax.set_title("gamma {:g} deg: {}".format(gam, label), fontsize=10)
            for (yy, xx), v in np.ndenumerate(M):
                if np.isfinite(v):
                    ax.text(xx, yy, "{:.1f}".format(v), ha="center", va="center", fontsize=8, color="w")
            fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("{}: rms vertical angle (top) and transmission (bottom) per shape, each at its angle-ranked quad setting".format(args.name))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "scan.png"), dpi=120)
except Exception as exc:
    stamp("figure failed: {}".format(exc))

stamp("scan done: {} points, {} failed; winner {}".format(len(rows), len(failed), ranked[0]["name"] if ranked else "-"))

# ---------------------------------------------------------------- the final push for the winner
if args.final and ranked:
    fname = args.final_name or "{}_final".format(args.name)
    stamp("final push '{}' for the winner {}".format(fname, ranked[0]["name"]))
    run([os.path.join(SCRIPTS, "RunFinalPush.py"), "--name", fname, "--bfield", args.bfield, "--particles", args.particles,
         "--knobs-json", os.path.join(OUT, "winner_knobs.json"), "--free-exit-truncation", "--rank", "vert", "--rank-tol", args.rank_tol,
         "--no-q2-exit-plate", "--fix-truncations", args.fix_entrance, 0.77], os.path.join(OUT, "log_final.txt"))
    stamp("final push done -> {}".format(os.path.join(DECK, "Results", "final", fname)))
stamp("all done ({:.0f} min)".format((time.time() - T0) / 60))
