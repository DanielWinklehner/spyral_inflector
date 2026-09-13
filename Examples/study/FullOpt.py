"""Full DFO-LS optimization of the doublet system with 8 mA space charge in the loop, seeded at a
final-push run's best setting. Knobs (10): k' tilt, angling, gamma, V-depth sigma, entrance and exit
truncation, axial shift dz, quad voltages q1 / q2, spiral gap (the design voltage scales with the gap,
the operating voltage stays the inner gap-centering solve). Each evaluation rebuilds the geometry at
the run's rotation R (design orbit with both truncations and dz fixed: only the centering solve runs),
re-solves the spiral basis field on the run's basis grid (the quad bases are reused), tracks --n core
particles with PyAMG space charge at the quad setting by superposition, and returns the residual vector

    [ rms vertical angle / 3 mrad,  max(0, T_floor - T) / 1 point,  centroid height / 0.3 mm,  centroid angle / 2 mrad ]

at 55 mm past the exit (T_floor = the seed's transmission - 5 points unless --t-floor is given).
DFO-LS (noise-aware, bound-scaled) minimizes the sum of squares. Every evaluation is cached in
<run>/full_opt/evals.json (resumable); the big basis files of each evaluation are deleted after the
measurement. Writes <run>/full_opt.md and <run>/full_opt/best_knobs.json (knobs for RunFinalPush
--knobs-json plus _truncations_deg, _dz_m, _q1, _q2).

    python FullOpt.py --run <deck>\\Results\\final\\vfocus1_level --best-json basis\\scan2_coarse.json --maxfun 80
"""
import argparse
import hashlib
import json
import os
import pickle
import shutil
import subprocess
import sys
import time

import numpy as np

DECK = os.environ.get("SI_DECK", r"D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron\Spiral_inflector")
SCRIPTS = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
sys.path.insert(0, os.path.dirname(os.path.dirname(SCRIPTS)))

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run", required=True, help="final-push run folder (geometry/, steps/, basis/, exit_plane.json, bfield_rotated.pickle)")
p.add_argument("--base-run", default=None, help="folder with bfield_rotated.pickle (default: --run)")
p.add_argument("--best-json", default="basis\\scan2_coarse.json", help="scan json whose 'best' holds the seed q1/q2")
p.add_argument("--quads", type=float, nargs=2, default=[16.0, 24.0], help="basis rotation of the quads on top of R [deg]")
p.add_argument("--n", type=int, default=8000)
p.add_argument("--h", type=float, default=0.004)
p.add_argument("--resolve-every", type=int, default=16)
p.add_argument("--current-ma", type=float, default=8.0)
p.add_argument("--geo-h", type=float, default=0.005)
p.add_argument("--maxfun", type=int, default=80)
p.add_argument("--rhobeg", type=float, default=0.12, help="initial trust region as a fraction of the bound ranges")
p.add_argument("--rhoend", type=float, default=0.01)
p.add_argument("--t-floor", type=float, default=None, help="transmission floor (fraction) of the penalty; default seed - 0.05")
p.add_argument("--bounds", type=float, nargs=20, default=[20, 42, 8, 22, 3, 13, 0.2, 1.6, 0.1, 0.8, 0.3, 1.6, 2.0, 7.0, 5500, 7500, -7700, -5700, 16, 24],
               metavar="B", help="lo hi per knob: tilt[deg] angling[deg] gamma[deg] sigma[mm] t_ent[deg] t_exit[deg] dz[mm] q1[V] q2[V] gap[mm]")
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle"))
p.add_argument("--keep-basis", action="store_true", help="keep every evaluation's basis field (large)")
a = p.parse_args()

NAMES = ["tilt", "angling", "gamma", "sigma_mm", "t_ent", "t_exit", "dz_mm", "q1", "q2", "gap_mm"]
UNITS = ["deg", "deg", "deg", "mm", "deg", "deg", "mm", "V", "V", "mm"]
LO = np.array(a.bounds[0::2], dtype=float)
HI = np.array(a.bounds[1::2], dtype=float)

RUN = os.path.abspath(a.run)
BASE = os.path.abspath(a.base_run) if a.base_run else RUN
BASIS0, STEPS0 = os.path.join(RUN, "basis"), os.path.join(RUN, "steps")
OUT = os.path.join(RUN, "full_opt")
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def jload(fn):
    with open(fn) as fh:
        return json.load(fh)


def run_cmd(cmd, log_path):
    with open(log_path, "w", encoding="utf-8") as fh:
        r = subprocess.run([PY, "-u"] + [str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS)
    if r.returncode != 0:
        raise RuntimeError("{} failed (exit {}), see {}".format(cmd[0], r.returncode, log_path))


from spyral_inflector.tracking.deck import particle_rows  # noqa: E402
from spyral_inflector.tracking.geometry_point import build_geometry, rmtree_retry  # noqa: E402
from PyPATools.field import Field  # noqa: E402

# ---------------------------------------------------------------- the seed
geo0 = jload(os.path.join(RUN, "geometry", "summary.json"))
knobs0 = dict(geo0["knobs"])
opt0 = geo0["optimizer"]
t_ent0, t_exit0 = (opt0.get("truncations_deg") or [0.34, 0.77])
dz0 = 1e-3 * opt0["dz_mm"]                              # the geometry is always rebuilt at dz0; dz acts as a tracking-time shift
lvl = os.path.join(RUN, "sc_level", "sc_level.json")
t_exit_seed, dz_seed = t_exit0, dz0
if os.path.exists(lvl):
    d = jload(lvl)
    t_exit_seed, dz_seed = float(d["t_exit_deg"]), float(d["dz_m"])
R = jload(os.path.join(RUN, "exit_plane.json"))["rotation_deg"]
PHI = 90.0 + R
ROT_BF = os.path.join(BASE, "bfield_rotated.pickle")
ENERGY = float(particle_rows(a.particles, core=True)[:, 8].mean())
bj = a.best_json if os.path.isabs(a.best_json) else os.path.join(RUN, a.best_json)
best0 = jload(bj)["best"]
QROT = list(a.quads)
_g = Field.from_file(os.path.join(BASIS0, "ef_itp_spiral.pickle")).grid
GRID_RES = float(np.asarray(_g["x"], dtype=float)[1] - np.asarray(_g["x"], dtype=float)[0])
GRID_BOX = [float(np.asarray(_g[c], dtype=float)[i]) for c in "xyz" for i in (0, -1)]
VOLT_PER_M = float(knobs0["volt"]) / float(knobs0["gap"])   # the design voltage scales with the gap (same field)

x0 = np.array([knobs0["tilt"], knobs0["angling"], knobs0["gamma"], 1e3 * knobs0["sigma"], t_ent0, t_exit_seed, 1e3 * dz_seed,
               best0["q1"], best0["q2"], 1e3 * knobs0["gap"]], dtype=float)
x0 = np.clip(x0, LO, HI)
stamp("FULL OPT of {}: R {:+.3f} deg fixed, quads rotated {}, {} mA, {} particles, {:.0f} mm cells, basis grid {:.1f} mm; DFO-LS maxfun {}, rhobeg {} of range".format(
    os.path.basename(RUN), R, QROT, a.current_ma, a.n, 1e3 * a.h, 1e3 * GRID_RES, a.maxfun, a.rhobeg))
stamp("  seed: " + ", ".join("{} {:g} {}".format(n, v, u) for n, v, u in zip(NAMES, x0, UNITS)))
stamp("  bounds: " + ", ".join("{} [{:g}, {:g}]".format(n, lo, hi) for n, lo, hi in zip(NAMES, LO, HI)))

# ---------------------------------------------------------------- the evaluation cache
CACHE = os.path.join(OUT, "evals.json")
evals = jload(CACHE) if os.path.exists(CACHE) else []


def key_of(x):
    return hashlib.md5(json.dumps([round(float(v), 4) for v in x]).encode()).hexdigest()[:10]


def save_cache():
    with open(CACHE, "w") as fh:
        json.dump(evals, fh, indent=1)


def measure(out_dir, tag):
    s = jload(os.path.join(out_dir, "bunch_{}.json".format(tag)))
    d = np.load(os.path.join(out_dir, "bunch_{}.npz".format(tag)), allow_pickle=True)
    st = d["asym_state"]
    ok = d["crossed"] & np.all(np.isfinite(st), axis=1)
    if "is_tail" in d.files:
        ok &= ~d["is_tail"]
    z, v = st[ok, 2], st[ok, 3:6]
    zp = v[:, 2] / np.hypot(v[:, 0], v[:, 1])
    return {"z": float(z.mean()), "zp": float(zp.mean()), "z_rms": float(z.std()), "zp_rms": float(zp.std()),
            "transmission": float(s["transmission"]), "n": int(ok.sum()), "voltage_V": s.get("voltages", {}).get("SI_Anode")}


T_FLOOR = a.t_floor          # set after the seed evaluation if None


def residuals(m):
    if m is None:
        return np.array([10.0, 10.0, 10.0, 10.0])
    return np.array([1e3 * m["zp_rms"] / 3.0, max(0.0, (T_FLOOR or 0.0) - m["transmission"]) / 0.01, 1e3 * m["z"] / 0.3, 1e3 * m["zp"] / 2.0])


def evaluate(x):
    x = np.asarray(x, dtype=float)
    key = key_of(x)
    for e in evals:
        if e["key"] == key:
            return e
    idx = len(evals)
    tag = "e{:03d}".format(idx)
    edir = os.path.join(OUT, tag)
    gdir, sdir, bdir = (os.path.join(edir, d) for d in ("geometry", "steps", "basis"))
    tilt, angling, gamma, sigma_mm, t_ent, t_exit, dz_mm, q1, q2, gap_mm = [float(v) for v in x]
    knobs = dict(knobs0)
    knobs.update(tilt=tilt, angling=angling, gamma=gamma, sigma=1e-3 * sigma_mm, gap=1e-3 * gap_mm, volt=VOLT_PER_M * 1e-3 * gap_mm)
    t0 = time.time()
    m, err = None, None
    try:
        rmtree_retry(edir, log=stamp)
        os.makedirs(bdir, exist_ok=True)
        build_geometry(gdir, sdir, ROT_BF, ENERGY, knobs=knobs, fix_truncations=(t_ent, t_exit), fix_dz=dz0, maxiter=0,
                       res=0.005, h=a.geo_h, log=lambda s: None)
        run_cmd(["-m", "spyral_inflector.tracking.bem_reload", "--tag", "spiral", "--step-dir", sdir, "--voltages", os.path.join(gdir, "voltages.csv"),
                 "--state", os.path.join(gdir, "state.pickle"), "--out-dir", bdir, "--res", GRID_RES, "--box"] + GRID_BOX +
                ["--h", a.geo_h, "--bfield", a.bfield, "--rotate-all", R, "--rotate-quads"] + QROT + ["--quad-voltages", 0, 0, "--energy-mev", ENERGY],
                os.path.join(edir, "log_basis_spiral.txt"))
        for fn in ("ef_itp_q1.pickle", "ef_itp_q2.pickle"):
            shutil.copyfile(os.path.join(BASIS0, fn), os.path.join(bdir, fn))
        run_cmd(["-m", "spyral_inflector.tracking.bunch", "--tag", "spiral", "--reload-dir", bdir, "--step-dir", sdir,
                 "--particles", a.particles, "--bfield", a.bfield, "--phi", PHI, "--out-dir", edir, "--out-tag", tag,
                 "--basis-dir", bdir, "--superpose", q1, 0.0, q2, 0.0, "--current-ma", a.current_ma, "--handoff-frame", "machine",
                 "--handoff-distance", 0.0, "--n", a.n, "--sc", "--h", a.h, "--resolve-every", a.resolve_every,
                 "--shift-z", "{:.10f}".format(1e-3 * dz_mm - dz0), "--no-plot"], os.path.join(edir, "log_bunch.txt"))
        m = measure(edir, tag)
    except Exception as exc:               # a collision / failed build / failed solve counts as a bad point, not a crash
        err = str(exc).splitlines()[-1][:200] if str(exc) else repr(exc)
    if not a.keep_basis and os.path.isdir(bdir):
        rmtree_retry(bdir)
    npz = os.path.join(edir, "bunch_{}.npz".format(tag))
    if os.path.exists(npz):
        os.remove(npz)
    r = residuals(m)
    e = {"key": key, "index": idx, "tag": tag, "x": [float(v) for v in x], "measure": m, "error": err,
         "residuals": [float(v) for v in r], "objective": float(np.sum(r ** 2)), "wall_s": time.time() - t0}
    evals.append(e)
    save_cache()
    if m is None:
        stamp("  {} FAILED ({:.0f} s): {} | {}".format(tag, e["wall_s"], err, " ".join("{}={:g}".format(n, v) for n, v in zip(NAMES, x))))
    else:
        stamp("  {} ({:.0f} s): angle {:.1f} mrad rms, T {:.1f} %, centroid {:+.2f} mm / {:+.1f} mrad, V {} -> f {:.3f} | {}".format(
            tag, e["wall_s"], 1e3 * m["zp_rms"], 100 * m["transmission"], 1e3 * m["z"], 1e3 * m["zp"],
            "-" if m["voltage_V"] is None else "{:.0f}".format(m["voltage_V"]), e["objective"],
            " ".join("{}={:g}".format(n, v) for n, v in zip(NAMES, x))))
    return e


# ---------------------------------------------------------------- seed, floor, DFO-LS
seed = evaluate(x0)
if seed["measure"] is None:
    raise SystemExit("the seed evaluation failed: " + str(seed["error"]))
if T_FLOOR is None:
    T_FLOOR = seed["measure"]["transmission"] - 0.05
    for e in evals:                                      # residuals of cached points with the final floor
        e["residuals"] = [float(v) for v in residuals(e["measure"])]
        e["objective"] = float(np.sum(np.array(e["residuals"]) ** 2))
    save_cache()
stamp("  transmission floor {:.1f} % (loop scale); seed objective {:.3f}".format(100 * T_FLOOR, seed["objective"]))

import dfols  # noqa: E402


def objfun(x):
    return np.array(evaluate(x)["residuals"], dtype=float)


stamp("DFO-LS start: {} knobs, maxfun {}".format(len(x0), a.maxfun))
soln = dfols.solve(objfun, x0, bounds=(LO, HI), rhobeg=a.rhobeg, rhoend=a.rhoend, maxfun=a.maxfun,
                   objfun_has_noise=True, scaling_within_bounds=True)
stamp("DFO-LS done: flag {} ({}), {} evaluations, f {:.3f}".format(soln.flag, soln.msg, soln.nf, float(np.sum(soln.resid ** 2))))

# ---------------------------------------------------------------- report and best knobs (from the cache, best objective)
good = [e for e in evals if e["measure"] is not None]
best = min(good, key=lambda e: e["objective"])
xb = best["x"]
bk = dict(knobs0)
bk.update(tilt=xb[0], angling=xb[1], gamma=xb[2], sigma=1e-3 * xb[3], gap=1e-3 * xb[9], volt=VOLT_PER_M * 1e-3 * xb[9])
bk.update({"_truncations_deg": [xb[4], xb[5]], "_dz_m": 1e-3 * xb[6], "_q1": xb[7], "_q2": xb[8], "_objective": best["objective"],
           "_measure": best["measure"], "_seed_objective": seed["objective"], "_R_deg": R, "_tag": best["tag"]})
with open(os.path.join(OUT, "best_knobs.json"), "w") as fh:
    json.dump(bk, fh, indent=2)
L = ["# Full DFO-LS optimization of {} with 8 mA space charge in the loop".format(os.path.basename(RUN)), "",
     "Knobs: k' tilt, angling, gamma, V-depth, entrance and exit truncation, dz, q1, q2, gap (the design voltage scales with the gap; "
     "the operating voltage is the inner centering solve). Rotation R = {:+.3f} deg fixed; spiral basis re-solved per evaluation on the "
     "run's {:.1f} mm grid, quad bases reused; {} core particles, {} mA, {:.0f} mm cells. Residuals: rms vertical angle / 3 mrad, "
     "(floor - T) / 1 point below the floor of {:.1f} %, centroid height / 0.3 mm, centroid angle / 2 mrad. DFO-LS: {} evaluations, "
     "flag {} ({}).".format(R, 1e3 * GRID_RES, a.n, a.current_ma, 1e3 * a.h, 100 * T_FLOOR, soln.nf, soln.flag, soln.msg), "",
     "| | " + " | ".join(NAMES) + " | angle [mrad] | T | centroid [mm / mrad] | f |", "|---|" + "---|" * (len(NAMES) + 4)]
for label, e in (("seed", seed), ("best", best)):
    m = e["measure"]
    L.append("| {} | ".format(label) + " | ".join("{:g}".format(v) for v in e["x"]) +
             " | {:.1f} | {:.1f} % | {:+.2f} / {:+.1f} | {:.3f} |".format(1e3 * m["zp_rms"], 100 * m["transmission"], 1e3 * m["z"], 1e3 * m["zp"], e["objective"]))
L += ["", "## All evaluations (sorted by objective)", "", "| tag | " + " | ".join(NAMES) + " | angle | T | height | z' | f |", "|---|" + "---|" * (len(NAMES) + 5)]
for e in sorted(evals, key=lambda e: e["objective"]):
    m = e["measure"]
    L.append("| {} | ".format(e["tag"]) + " | ".join("{:g}".format(v) for v in e["x"]) +
             (" | {:.1f} | {:.1f} % | {:+.2f} | {:+.1f} | {:.3f} |".format(1e3 * m["zp_rms"], 100 * m["transmission"], 1e3 * m["z"], 1e3 * m["zp"], e["objective"])
              if m else " | FAILED | | | | {} |".format(e["error"])))
L += ["", "Best knobs for the final push: `full_opt/best_knobs.json` (RunFinalPush --knobs-json it --fix-truncations {:.3f} {:.3f} --fix-dz {:.6f} "
      "--q1 {:.0f} {:.0f} 1 --q2 {:.0f} {:.0f} 1 --fine-step 0), then the plate-rule R iterates to its fixed point there.".format(
          xb[4], xb[5], 1e-3 * xb[6], xb[7], xb[7], xb[8], xb[8])]
with open(os.path.join(RUN, "full_opt.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
stamp("report -> {}; best {} f {:.3f} vs seed {:.3f}".format(os.path.join(RUN, "full_opt.md"), best["tag"], best["objective"], seed["objective"]))
