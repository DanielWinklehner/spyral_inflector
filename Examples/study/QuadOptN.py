"""QuadOptN.py - DFO-LS over the quadrupole VOLTAGES AND ROTATIONS of an N-quad TripletScan run,
with space charge in the loop, on the run's own basis fields (normal + skew per quad, superposed;
no BEM solve per evaluation). The inflector is not touched.

Each evaluation is a bunch run like SCGridN's (core particles, --n of them, at --current-ma) with
quad i at voltage q_i and turned by alpha_i on top of the basis orientation (cos 2a normal +
sin 2a skew -- exact for the pure quadrupole part). Measured at the hand-off plane exactly like
FullOpt. Residuals (FullOpt's units):
    w_angle * rms vertical angle at the hand-off / 3 mrad
    (floor - T) / 1 point, one-sided, floor = seed transmission - --t-tol
    w_z100 * rms ballistic height 100 mm past the hand-off / 5 mm            (--w-z100 > 0)
    w_hv * fraction of the bunch lost on HV electrodes (spiral + quad poles) / 1 point   (--w-hv > 0)
No centroid terms: the centroid is set by the levelling (SCLevel) afterwards.

    python QuadOptN.py --run <TripletScan run> --bfield <pickle> --quads 16 24 24 --seed-q 8000 -8000 -2000
           [--seed-alpha 0 0 0] [--q-span 2500] [--alpha-span 15] [--maxfun 60] [--n 4000] [--current-ma 10]

Writes <run>/<out-name>/ (evals.json, log.txt, best.json in the best-json form SCLevel reads,
report.md). Resumable: cached evaluations are reused.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
PY = sys.executable
sys.path.insert(0, SCRIPTS)

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run", required=True, help="TripletScan run folder (geometry/, steps/, basis/ with normal AND skew bases, exit_plane.json)")
p.add_argument("--bfield", required=True)
p.add_argument("--quads", type=float, nargs="+", required=True, help="basis orientation of the quads on top of R [deg] (TripletScan --quads)")
p.add_argument("--seed-q", type=float, nargs="+", required=True, help="start voltages [V], one per quad")
p.add_argument("--seed-alpha", type=float, nargs="+", default=None, help="start extra rotations [deg] (default 0)")
p.add_argument("--q-span", type=float, default=2500.0, help="bounds: seed +- this [V]")
p.add_argument("--alpha-span", type=float, default=15.0, help="bounds: seed +- this [deg]")
p.add_argument("--maxfun", type=int, default=60)
p.add_argument("--rhobeg", type=float, default=0.15, help="initial trust region, fraction of the bound ranges")
p.add_argument("--rhoend", type=float, default=0.01)
p.add_argument("--t-tol", type=float, default=0.03, help="transmission floor = seed - this (fraction)")
p.add_argument("--w-angle", type=float, default=1.0)
p.add_argument("--w-z100", type=float, default=0.0)
p.add_argument("--w-hv", type=float, default=0.0,
               help="weight of the fraction lost on HV electrodes, spiral electrodes AND quad poles (1 point = 1 unit): Daniel wants the "
                    "losses on grounded apertures, so pushing them from the spiral onto a quad pole must not count as a gain")
p.add_argument("--n", type=int, default=4000)
p.add_argument("--h", type=float, default=0.004)
p.add_argument("--resolve-every", type=int, default=10)
p.add_argument("--current-ma", type=float, default=10.0)
p.add_argument("--vscale", type=float, default=1.0)
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--out-name", default="quad_opt")
p.add_argument("--keep", action="store_true", help="keep every evaluation's npz (large)")
a = p.parse_args()

RUN = os.path.abspath(a.run)
BDIR, SDIR = os.path.join(RUN, "basis"), os.path.join(RUN, "steps")
OUT = os.path.join(RUN, a.out_name)
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")
SI_ELECTRODES = ("SI_Anode", "SI_Cathode")


def is_hv(name):
    """The spiral electrodes and the quad poles (D0, D1, ...); plates, housing and the sleeve are grounded."""
    return name in SI_ELECTRODES or (name.startswith("D") and name[1:].isdigit())


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


geo0 = jload(os.path.join(RUN, "geometry", "summary.json"))
NQ = int(geo0["knobs"].get("n_quads", 2) or 2)
for name, v in (("--quads", a.quads), ("--seed-q", a.seed_q)):
    if len(v) != NQ:
        raise SystemExit("{} needs {} values for this {}-quad run".format(name, NQ, NQ))
seed_alpha = list(a.seed_alpha) if a.seed_alpha else [0.0] * NQ
need = ["ef_itp_spiral.pickle"] + ["ef_itp_q{}{}.pickle".format(i + 1, s) for i in range(NQ) for s in ("", "skew")]
missing = [f for f in need if not os.path.exists(os.path.join(BDIR, f))]
if missing:
    raise SystemExit("missing basis fields in {}: {} (TripletScan --skew-bases solves the skew ones)".format(BDIR, ", ".join(missing)))
R = jload(os.path.join(RUN, "exit_plane.json"))["rotation_deg"]
PHI = 90.0 + R
NAMES = ["q{}".format(i + 1) for i in range(NQ)] + ["alpha{}".format(i + 1) for i in range(NQ)]
x0 = np.array(list(a.seed_q) + seed_alpha, dtype=float)
LO = np.array([q - a.q_span for q in a.seed_q] + [s - a.alpha_span for s in seed_alpha])
HI = np.array([q + a.q_span for q in a.seed_q] + [s + a.alpha_span for s in seed_alpha])
stamp("QUAD OPT of {}: {} quads at basis rotation {} on top of R {:+.3f} deg; {} mA, {} core particles, {:.0f} mm cells; "
      "knobs {}; seed {}; bounds q +-{:.0f} V, alpha +-{:.0f} deg; weights angle {:g}, z100 {:g}, HV loss {:g}; DFO-LS maxfun {}".format(
          os.path.basename(RUN), NQ, list(a.quads), R, a.current_ma, a.n, 1e3 * a.h, NAMES,
          ", ".join("{} {:g}".format(n, v) for n, v in zip(NAMES, x0)), a.q_span, a.alpha_span, a.w_angle, a.w_z100, a.w_hv, a.maxfun))

CACHE = os.path.join(OUT, "evals.json")
evals = jload(CACHE) if os.path.exists(CACHE) else []


def key_of(x):
    return hashlib.md5(json.dumps([round(float(v), 3) for v in x]).encode()).hexdigest()[:10]


def save():
    with open(CACHE, "w") as fh:
        json.dump(evals, fh, indent=1)


def measure(out_dir, tag):
    """FullOpt.measure (as in SCGridN), plus the loss split by electrode."""
    s = jload(os.path.join(out_dir, "bunch_{}.json".format(tag)))
    d = np.load(os.path.join(out_dir, "bunch_{}.npz".format(tag)), allow_pickle=True)
    st = d["asym_state"]
    ok = d["crossed"] & np.all(np.isfinite(st), axis=1)
    if "is_tail" in d.files:
        ok &= ~d["is_tail"]
    z, v = st[ok, 2], st[ok, 3:6]
    zp = v[:, 2] / np.hypot(v[:, 0], v[:, 1])
    hs = d["handoff_state"] if "handoff_state" in d.files else None
    okh = ok & np.all(np.isfinite(hs), axis=1) if hs is not None else np.zeros(len(ok), bool)
    M_H2P_KG = 3.3462e-27
    dT_rms_keV, f_low, dt_rms_deg = 0.0, 0.0, 0.0
    if okh.sum() >= 100:
        zh, vh = hs[okh, 2], hs[okh, 3:6]
        zph = vh[:, 2] / np.hypot(vh[:, 0], vh[:, 1])
        zp_h_rms, z100 = float(zph.std()), zh + 0.1 * zph
        T_keV = 0.5 * M_H2P_KG * np.sum(vh ** 2, axis=1) / 1.602176634e-19 / 1e3
        dT_rms_keV, f_low = float(T_keV.std()), float(np.mean(T_keV < 0.965 * T_keV.mean()))
        if "handoff_time" in d.files:
            th = np.asarray(d["handoff_time"], dtype=float)[okh]
            th = th[np.isfinite(th)]
            f_rf = float(s.get("rf_mhz") or s.get("rf_frequency_mhz") or 32.8) * 1e6
            if len(th) >= 100:
                dt_rms_deg = float(360.0 * f_rf * th.std())
    else:
        zp_h_rms, z100 = float(zp.std()), z + 0.1 * zp
    n_all = float(s["n_particles"])
    lbe = s.get("losses_by_electrode", {})
    return {"z": float(z.mean()), "zp": float(zp.mean()), "z_rms": float(z.std()), "zp_rms": float(zp.std()), "zp_h_rms": zp_h_rms,
            "z100_rms": float(z100.std()), "dT_rms_keV": dT_rms_keV, "f_low_energy": f_low, "dt_rms_deg": dt_rms_deg,
            "transmission": float(s["transmission"]), "n": int(ok.sum()), "voltage_V": s.get("voltages", {}).get("SI_Anode"),
            "f_lost_si": sum(lbe.get(k, 0) for k in SI_ELECTRODES) / n_all,
            "f_lost_hv": sum(v_ for k, v_ in lbe.items() if is_hv(k)) / n_all,
            "losses_by_electrode": {k: v_ / n_all for k, v_ in lbe.items()}}


T_FLOOR = None


def residuals(m):
    nres = 2 + (a.w_z100 > 0) + (a.w_hv > 0)
    if m is None:
        return np.full(nres, 10.0)
    r = [a.w_angle * 1e3 * m["zp_h_rms"] / 3.0, max(0.0, (T_FLOOR or 0.0) - m["transmission"]) / 0.01]
    if a.w_z100 > 0:
        r.append(a.w_z100 * 1e3 * m["z100_rms"] / 5.0)
    if a.w_hv > 0:
        r.append(a.w_hv * 100.0 * m["f_lost_hv"])
    return np.array(r)


def evaluate(x):
    x = np.asarray(x, dtype=float)
    key = key_of(x)
    for e in evals:
        if e["key"] == key and e["measure"] is not None:
            return e
    idx = len(evals)
    tag = "q{:03d}".format(idx)
    edir = os.path.join(OUT, tag)
    SUPER = [v for i in range(NQ) for v in (x[i], x[NQ + i])]
    t0 = time.time()
    m, err = None, None
    try:
        if os.path.isdir(edir):
            shutil.rmtree(edir, ignore_errors=True)
        os.makedirs(edir, exist_ok=True)
        run_cmd(["-m", "spyral_inflector.tracking.bunch", "--tag", "spiral", "--reload-dir", BDIR, "--step-dir", SDIR,
                 "--particles", a.particles, "--bfield", a.bfield, "--phi", PHI, "--out-dir", edir, "--out-tag", tag,
                 "--basis-dir", BDIR, "--superpose"] + SUPER + ["--current-ma", a.current_ma, "--handoff-frame", "machine",
                 "--handoff-distance", 0.0, "--n", a.n] + (["--sc"] if a.current_ma > 0 else []) +
                ["--h", a.h, "--resolve-every", a.resolve_every, "--shift-z", "0.0000000000", "--vscale", "{:.6f}".format(a.vscale),
                 "--no-plot"], os.path.join(edir, "log_bunch.txt"))
        m = measure(edir, tag)
    except Exception as exc:
        err = str(exc).splitlines()[-1][:200] if str(exc) else repr(exc)
    if not a.keep:
        npz = os.path.join(edir, "bunch_{}.npz".format(tag))
        if os.path.exists(npz):
            os.remove(npz)
    r = residuals(m)
    e = {"key": key, "index": idx, "tag": tag, "x": [float(v) for v in x], "measure": m, "error": err,
         "residuals": [float(v) for v in r], "objective": float(np.sum(r ** 2)), "wall_s": time.time() - t0}
    evals.append(e)
    save()
    knobs = " ".join("{}={:.6g}".format(n, v) for n, v in zip(NAMES, x))
    if m is None:
        stamp("  {} FAILED ({:.0f} s): {} | {}".format(tag, e["wall_s"], err, knobs))
    else:
        stamp("  {} ({:.0f} s): angle {:.1f} mrad at the hand-off ({:.1f} at 55 mm), T {:.1f} %, lost on SI electrodes {:.1f} % / on all HV {:.1f} %, centroid "
              "{:+.2f} mm / {:+.1f} mrad, proxy {:.1f} mm, dE {:.2f} keV, dt {:.1f} deg -> f {:.3f} | {}".format(
                  tag, e["wall_s"], 1e3 * m["zp_h_rms"], 1e3 * m["zp_rms"], 100 * m["transmission"], 100 * m["f_lost_si"], 100 * m["f_lost_hv"],
                  1e3 * m["z"], 1e3 * m["zp"], 1e3 * m["z100_rms"], m["dT_rms_keV"], m["dt_rms_deg"], e["objective"], knobs))
    return e


seed = evaluate(x0)
if seed["measure"] is None:
    raise SystemExit("the seed evaluation failed: " + str(seed["error"]))
T_FLOOR = seed["measure"]["transmission"] - a.t_tol
for e in evals:
    e["residuals"] = [float(v) for v in residuals(e["measure"])]
    e["objective"] = float(np.sum(np.array(e["residuals"]) ** 2))
save()
stamp("  transmission floor {:.1f} %; seed objective {:.3f}".format(100 * T_FLOOR, seed["objective"]))

import dfols  # noqa: E402

soln = dfols.solve(lambda xx: np.array(evaluate(xx)["residuals"]), x0, bounds=(LO, HI), rhobeg=a.rhobeg, rhoend=a.rhoend,
                   maxfun=a.maxfun, objfun_has_noise=True, scaling_within_bounds=True)
stamp("DFO-LS done: flag {} ({}), {} evaluations".format(soln.flag, soln.msg, soln.nf))

good = sorted([e for e in evals if e["measure"] is not None], key=lambda e: e["objective"])
best = good[0]
xb = best["x"]
bj = {"best": {**{"q{}".format(i + 1): xb[i] for i in range(NQ)}, **{"alpha{}".format(i + 1): xb[NQ + i] for i in range(NQ)},
               "transmission": best["measure"]["transmission"], "vfom_h": 1e3 * best["measure"]["zp_h_rms"]},
      "quad_rotation_deg": [float(q) + xb[NQ + i] for i, q in enumerate(a.quads)], "basis_rotation_deg": list(a.quads),
      "measure": best["measure"], "tag": best["tag"], "objective": best["objective"]}
with open(os.path.join(OUT, "best.json"), "w") as fh:
    json.dump(bj, fh, indent=2)
L = ["# Quad voltages and rotations of {} ({} mA in the loop)".format(os.path.basename(RUN), a.current_ma), "",
     "Basis orientation {} deg on top of R = {:+.3f} deg; the total rotation of quad i is basis + alpha_i. {} core particles per "
     "evaluation (the tail is not in the loop; every-particle runs confirm the result). DFO-LS: {} evaluations, flag {} ({}).".format(
         list(a.quads), R, a.n, soln.nf, soln.flag, soln.msg), "",
     "| | " + " | ".join(NAMES) + " | angle hand-off [mrad] | T | lost on SI electrodes | f |", "|---|" + "---|" * (len(NAMES) + 4)]
for label, e in (("seed", seed), ("best", best)):
    m = e["measure"]
    L.append("| {} | ".format(label) + " | ".join("{:.4g}".format(v) for v in e["x"]) +
             " | {:.1f} | {:.1f} % | {:.1f} % | {:.3f} |".format(1e3 * m["zp_h_rms"], 100 * m["transmission"], 100 * m["f_lost_si"], e["objective"]))
L += ["", "Total quad rotations at the best point: " + ", ".join("q{} {:+.2f} deg".format(i + 1, r) for i, r in enumerate(bj["quad_rotation_deg"])),
      "", "## All evaluations, by objective", "", "| tag | " + " | ".join(NAMES) + " | angle | T | SI loss | f |", "|---|" + "---|" * (len(NAMES) + 4)]
for e in good:
    m = e["measure"]
    L.append("| {} | ".format(e["tag"]) + " | ".join("{:.4g}".format(v) for v in e["x"]) +
             " | {:.1f} | {:.1f} % | {:.1f} % | {:.3f} |".format(1e3 * m["zp_h_rms"], 100 * m["transmission"], 100 * m["f_lost_si"], e["objective"]))
with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
stamp("best {} f {:.3f} (seed {:.3f}); rotations {}; report -> {}".format(best["tag"], best["objective"], seed["objective"],
                                                                          [round(r, 2) for r in bj["quad_rotation_deg"]], os.path.join(OUT, "report.md")))
