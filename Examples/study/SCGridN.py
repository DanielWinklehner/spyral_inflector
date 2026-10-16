"""SCGridN.py - a space-charge quad grid for an N-quad TripletScan run, ranked three ways.

Every (q1, q2, q3[, vscale]) point is a full bunch run WITH space charge on the run's own basis fields (superposed, no
BEM solve per point), measured at the hand-off plane exactly like FullOpt: rms vertical angle z', transmission, centroid
height / angle, aperture proxy (rms ballistic height 100 mm on), energy spread, and the bunch length as the rms arrival
time in RF degrees. Then three picks within --tol of the best transmission (Daniel, 2026-09-15 19:55: "2-3 versions
emphasizing a) min vertical divergence, b) min phase width, c) balanced; all should include transmission as well, but
not as the main goal"):
  a  minimum rms vertical angle at the hand-off
  b  minimum rms arrival time (RF degrees)
  c  minimum of angle / 3 mrad + bunch length / 10 deg + energy spread / 1 keV + proxy / 5 mm (FullOpt's residual units)
Each pick is written as <run>/sc_grid/pick_<x>.json in the best-json form SCLevel reads ({"best": {"q1", "q2", "q3"}},
plus "vscale" and the measure), the whole grid as sc_grid/points.json (resume-safe: finished points are skipped) and
sc_grid/report.md.

  python SCGridN.py --run <TripletScan run> --base-run <RunFinalPush run> --bfield <pickle> --quads 16 24 24
         --q1 7500 8100 8700 --q2 -8500 -7900 -7300 --q3 -5200 -4600 -4000 [--vscale 1.0] [--n 4000] [--current-ma 10]
         [--h 0.004] [--resolve-every 10] [--tol 0.05] [--particles <dst>]
"""
import argparse
import itertools
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

DECK = os.environ.get("SI_DECK", r"D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron\Spiral_inflector")
SCRIPTS = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

p = argparse.ArgumentParser()
p.add_argument("--run", required=True, help="TripletScan run folder (geometry/, steps/, basis/ with ef_itp_*.pickle, exit_plane.json)")
p.add_argument("--base-run", default=None, help="folder with bfield_rotated.pickle (default: --run)")
p.add_argument("--quads", type=float, nargs="+", default=[16.0, 24.0, 24.0], help="basis rotation of the quads on top of R [deg]")
p.add_argument("--q1", type=float, nargs="+", required=True, help="quad 1 voltages [V]")
p.add_argument("--q2", type=float, nargs="+", required=True)
p.add_argument("--q3", type=float, nargs="+", default=None, help="quad 3 voltages [V] (3-quad runs)")
p.add_argument("--vscale", type=float, nargs="+", default=[1.0], help="operating spiral voltage / centering value")
p.add_argument("--n", type=int, default=4000)
p.add_argument("--h", type=float, default=0.004)
p.add_argument("--resolve-every", type=int, default=10)   # 10 since 2026-09-15 23:50 (was 16)
p.add_argument("--current-ma", type=float, default=10.0)
p.add_argument("--tol", type=float, default=0.05, help="transmission tolerance of the rankings (fraction below the best)")
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", required=True)
p.add_argument("--out-name", default="sc_grid")
p.add_argument("--keep", action="store_true", help="keep every point's npz (large)")
a = p.parse_args()

RUN = os.path.abspath(a.run)
BASE = os.path.abspath(a.base_run) if a.base_run else RUN
BDIR, SDIR = os.path.join(RUN, "basis"), os.path.join(RUN, "steps")
OUT = os.path.join(RUN, a.out_name)
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

geo0 = jload(os.path.join(RUN, "geometry", "summary.json"))
knobs0 = dict(geo0["knobs"])
NQ = int(knobs0.get("n_quads", 2) or 2)
if NQ >= 3 and not a.q3:
    raise SystemExit("this is a {}-quad run: give --q3".format(NQ))
QROT = list(a.quads)
if len(QROT) != NQ:
    raise SystemExit("--quads needs {} angles, got {}".format(NQ, len(QROT)))
R = jload(os.path.join(RUN, "exit_plane.json"))["rotation_deg"]
PHI = 90.0 + R
ENERGY = float(particle_rows(a.particles, core=True)[:, 8].mean())
for fn in ["ef_itp_spiral.pickle"] + ["ef_itp_q{}.pickle".format(i + 1) for i in range(NQ)]:
    if not os.path.exists(os.path.join(BDIR, fn)):
        raise SystemExit("missing basis field {} in {}".format(fn, BDIR))
axes = [a.q1, a.q2] + ([a.q3] if NQ >= 3 else []) + [a.vscale]
points = [tuple(float(v) for v in t) for t in itertools.product(*axes)]
# 2026-09-15 23:00 (Daniel): hold here while the space-charge fixes from the audits go in. The chain launches this
# script with a fixed command line, so the pause is a file: <RUN>/HOLD.txt or <DECK>/Results/HOLD_SCGRID.txt.
_holds = [os.path.join(RUN, "HOLD.txt"), os.path.join(DECK, "Results", "HOLD_SCGRID.txt")]
_said = False
while any(os.path.exists(h) for h in _holds):
    if not _said:
        stamp("HOLD: waiting while one of these exists: " + ", ".join(_holds))
        _said = True
    time.sleep(60)
if _said:
    stamp("hold released")
stamp("SC GRID of {}: R {:+.3f} deg, quads rotated {}, {} mA, {} particles, {:.0f} mm cells, {} points ({} x {}{} x {} vscale), tol {:.0f} points".format(
    os.path.basename(RUN), R, QROT, a.current_ma, a.n, 1e3 * a.h, len(points), len(a.q1), len(a.q2), (" x {}".format(len(a.q3)) if NQ >= 3 else ""),
    len(a.vscale), 100 * a.tol))

PTS = os.path.join(OUT, "points.json")
done = jload(PTS) if os.path.exists(PTS) else []


def key_of(pt):
    return "|".join("{:.3f}".format(v) for v in pt)


def save():
    with open(PTS, "w") as fh:
        json.dump(done, fh, indent=1)


def measure(out_dir, tag):
    """FullOpt.measure with the bunch-length term (kept in step with Examples/study/FullOpt.py)."""
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
    return {"z": float(z.mean()), "zp": float(zp.mean()), "z_rms": float(z.std()), "zp_rms": float(zp.std()), "zp_h_rms": zp_h_rms,
            "z100_rms": float(z100.std()), "dT_rms_keV": dT_rms_keV, "f_low_energy": f_low, "dt_rms_deg": dt_rms_deg,
            "transmission": float(s["transmission"]), "n": int(ok.sum()), "voltage_V": s.get("voltages", {}).get("SI_Anode")}


for i, pt in enumerate(points):
    key = key_of(pt)
    if any(e["key"] == key and e.get("measure") for e in done):
        continue
    tag = "p{:03d}".format(i)
    edir = os.path.join(OUT, tag)
    qv, vscale = list(pt[:NQ]), float(pt[NQ])
    SUPER = [v for q in qv for v in (q, 0.0)]
    t0 = time.time()
    m, err = None, None
    try:
        if os.path.isdir(edir):
            shutil.rmtree(edir, ignore_errors=True)
        os.makedirs(edir, exist_ok=True)
        run_cmd(["-m", "spyral_inflector.tracking.bunch", "--tag", "spiral", "--reload-dir", BDIR, "--step-dir", SDIR,
                 "--particles", a.particles, "--bfield", a.bfield, "--phi", PHI, "--out-dir", edir, "--out-tag", tag,
                 "--basis-dir", BDIR, "--superpose"] + SUPER + ["--current-ma", a.current_ma, "--handoff-frame", "machine",
                 "--handoff-distance", 0.0, "--n", a.n, "--sc", "--h", a.h, "--resolve-every", a.resolve_every,
                 "--shift-z", "0.0000000000", "--vscale", "{:.6f}".format(vscale), "--no-plot"], os.path.join(edir, "log_bunch.txt"))
        m = measure(edir, tag)
    except Exception as exc:
        err = str(exc).splitlines()[-1][:200] if str(exc) else repr(exc)
    if not a.keep:
        npz = os.path.join(edir, "bunch_{}.npz".format(tag))
        if os.path.exists(npz):
            os.remove(npz)
    e = {"key": key, "tag": tag, "q": qv, "vscale": vscale, "measure": m, "error": err, "wall_s": time.time() - t0}
    done = [d_ for d_ in done if d_["key"] != key] + [e]
    save()
    if m is None:
        stamp("  {} FAILED ({:.0f} s): {} | q {} vscale {:g}".format(tag, e["wall_s"], err, [round(q) for q in qv], vscale))
    else:
        stamp("  {} ({:.0f} s): angle {:.1f} mrad at the hand-off ({:.1f} at 55 mm), T {:.1f} %, centroid {:+.2f} mm / {:+.1f} mrad, proxy {:.1f} mm, "
              "dE {:.2f} keV ({:.0f} % low), dt {:.1f} deg, V {} | q {} vscale {:g}".format(
                  tag, e["wall_s"], 1e3 * m["zp_h_rms"], 1e3 * m["zp_rms"], 100 * m["transmission"], 1e3 * m["z"], 1e3 * m["zp"], 1e3 * m["z100_rms"],
                  m["dT_rms_keV"], 100 * m["f_low_energy"], m["dt_rms_deg"], "-" if m["voltage_V"] is None else "{:.0f}".format(m["voltage_V"]),
                  [round(q) for q in qv], vscale))

# ---------------------------------------------------------------- rankings
good = [e for e in done if e.get("measure")]
if not good:
    raise SystemExit("no successful points")
t_best = max(e["measure"]["transmission"] for e in good)
floor = t_best - a.tol
pool = [e for e in good if e["measure"]["transmission"] >= floor]


def score_c(m):
    return 1e3 * m["zp_h_rms"] / 3.0 + m["dt_rms_deg"] / 10.0 + m["dT_rms_keV"] / 1.0 + 1e3 * m["z100_rms"] / 5.0


picks = {"a": min(pool, key=lambda e: e["measure"]["zp_h_rms"]),
         "b": min(pool, key=lambda e: e["measure"]["dt_rms_deg"]),
         "c": min(pool, key=lambda e: score_c(e["measure"]))}
labels = {"a": "minimum vertical angle at the hand-off", "b": "minimum bunch length (rms arrival time, RF degrees)",
          "c": "balanced: angle / 3 mrad + bunch length / 10 deg + energy spread / 1 keV + proxy / 5 mm"}
for k, e in picks.items():
    with open(os.path.join(OUT, "pick_{}.json".format(k)), "w") as fh:
        json.dump({"best": {"q{}".format(i + 1): e["q"][i] for i in range(NQ)}, "vscale": e["vscale"], "tag": e["tag"], "rank": k,
                   "label": labels[k], "measure": e["measure"], "transmission_floor": floor, "score_c": score_c(e["measure"])}, fh, indent=2)
    m = e["measure"]
    stamp("pick {} ({}): {} q {} vscale {:g}: angle {:.1f} mrad, dt {:.1f} deg, dE {:.2f} keV, proxy {:.1f} mm, T {:.1f} %".format(
        k, labels[k], e["tag"], [round(q) for q in e["q"]], e["vscale"], 1e3 * m["zp_h_rms"], m["dt_rms_deg"], m["dT_rms_keV"], 1e3 * m["z100_rms"], 100 * m["transmission"]))
qn = ["q{}".format(i + 1) for i in range(NQ)]
L = ["# Space-charge quad grid of {} ({} mA, {} particles, {:.0f} mm cells, {} points)".format(os.path.basename(RUN), a.current_ma, a.n, 1e3 * a.h, len(good)), "",
     "R = {:+.3f} deg, quads rotated {}, bases of the run superposed (no BEM solve per point), unlevelled geometry (design orbit). "
     "Measured at the hand-off plane. Rankings within {:.0f} points of the best transmission ({:.1f} %, floor {:.1f} %).".format(R, QROT, 100 * a.tol, 100 * t_best, 100 * floor), "",
     "| pick | rule | tag | " + " | ".join(qn) + " | vscale | angle [mrad] | dt [deg] | dE [keV] | proxy [mm] | T |", "|---|---|---|" + "---|" * (NQ + 6)]
for k, e in picks.items():
    m = e["measure"]
    L.append("| {} | {} | {} | ".format(k, labels[k], e["tag"]) + " | ".join("{:.0f}".format(q) for q in e["q"]) +
             " | {:g} | {:.1f} | {:.1f} | {:.2f} | {:.1f} | {:.1f} % |".format(e["vscale"], 1e3 * m["zp_h_rms"], m["dt_rms_deg"], m["dT_rms_keV"], 1e3 * m["z100_rms"], 100 * m["transmission"]))
L += ["", "## All points (sorted by the balanced score)", "",
      "| tag | " + " | ".join(qn) + " | vscale | angle hand-off | angle 55 mm | T | height | z' | proxy | dE | low E | dt [deg] | score c |", "|---|" + "---|" * (NQ + 11)]
for e in sorted(good, key=lambda e: score_c(e["measure"])):
    m = e["measure"]
    L.append("| {} | ".format(e["tag"]) + " | ".join("{:.0f}".format(q) for q in e["q"]) +
             " | {:g} | {:.1f} | {:.1f} | {:.1f} % | {:+.2f} | {:+.1f} | {:.1f} | {:.2f} | {:.0f} % | {:.1f} | {:.2f} |".format(
                 e["vscale"], 1e3 * m["zp_h_rms"], 1e3 * m["zp_rms"], 100 * m["transmission"], 1e3 * m["z"], 1e3 * m["zp"], 1e3 * m["z100_rms"],
                 m["dT_rms_keV"], 100 * m["f_low_energy"], m["dt_rms_deg"], score_c(m)))
failed = [e for e in done if not e.get("measure")]
if failed:
    L += ["", "Failed: " + ", ".join("{} ({})".format(e["tag"], e["error"]) for e in failed)]
with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
stamp("report -> {}".format(os.path.join(OUT, "report.md")))
