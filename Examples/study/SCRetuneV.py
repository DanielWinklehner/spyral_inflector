"""Space charge in the loop: re-tune the quads of a finished final-push run WITH space charge,
ranked by the rms vertical angle, then re-run the final space-charge run at that setting.

The final push tunes the quads in vacuum (superposed basis fields, 1500 particles) and only
evaluates the 8 mA case once. Here a grid of quad voltages around the vacuum choice is
tracked with PyAMG space charge (--current-ma, --h cells, --n core particles) on the same
superposed fields; each point's transmission and the rms vertical angle vz/v_long of the
transmitted bunch at the asymptotic state are taken from the run; the best point is the
smallest angle within --tol of the best transmission. That setting is then run with every
particle of the RFQ file (core + injected tail) with space charge at the final cell size,
and once in vacuum for reference, both with hand-off files.

    python SCRetuneV.py --run ..\\Results\\final\\vfocus1_final
Writes <run>/sc_retune/ (one bunch run per grid point, bunch_sc_all_scretune.*, handoff_sc_all_scretune.h5)
and <run>/sc_retune.md.
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
PY = sys.executable
C = 299792458.0

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run", required=True, help="finished RunFinalPush folder")
p.add_argument("--tol", type=float, default=0.05)
p.add_argument("--dq1", type=float, nargs="+", default=[-500, 0, 500], help="q1 offsets from the vacuum choice [V]")
p.add_argument("--dq2", type=float, nargs="+", default=[-500, -250, 0, 250, 500], help="q2 offsets from the vacuum choice [V]")
p.add_argument("--n", type=int, default=5000, help="core particles per grid point")
p.add_argument("--h", type=float, default=0.004, help="Poisson cell size of the grid runs [m]")
p.add_argument("--resolve-every", type=int, default=16)
p.add_argument("--final-h", type=float, default=0.002)
p.add_argument("--final-resolve-every", type=int, default=8)
p.add_argument("--current-ma", type=float, default=8.0)
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=None, help="default: the run's field (final1_baseline_1mm.pickle)")
p.add_argument("--skip-final", action="store_true")
a = p.parse_args()

RUN = os.path.abspath(a.run)
BASIS, STEPS = os.path.join(RUN, "basis"), os.path.join(RUN, "steps")
OUT = os.path.join(RUN, "sc_retune")
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")
BF = a.bfield or os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle")


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def jload(fn):
    with open(fn) as fh:
        return json.load(fh)


R = jload(os.path.join(RUN, "exit_plane.json"))["rotation_deg"]
PHI = 90.0 + R
best_vac = jload(os.path.join(BASIS, "scan2_fine.json"))["best"]
q1_0, q2_0 = best_vac["q1"], best_vac["q2"]
stamp("SC RETUNE of {}: vacuum choice q1 {:+.0f} / q2 {:+.0f} V ({:.1f} %, {} mrad); grid {} x {} at {} mA, {} particles, {:.0f} mm cells".format(
    os.path.basename(RUN), q1_0, q2_0, 100 * best_vac["transmission"],
    "-" if best_vac.get("vfom") is None else "{:.1f}".format(best_vac["vfom"]), len(a.dq1), len(a.dq2), a.current_ma, a.n, 1e3 * a.h))


def bunch(tag, extra, log_name):
    cmd = [PY, "-u", "-m", "spyral_inflector.tracking.bunch", "--tag", "fine_best", "--reload-dir", BASIS, "--step-dir", STEPS,
           "--particles", a.particles, "--bfield", BF, "--phi", str(PHI), "--out-dir", OUT, "--out-tag", tag,
           "--basis-dir", BASIS, "--current-ma", str(a.current_ma), "--handoff-frame", "machine", "--handoff-distance", "0.0"] + [str(x) for x in extra]
    with open(os.path.join(OUT, log_name), "w", encoding="utf-8") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS)
    if r.returncode != 0:
        raise RuntimeError("bunch run {} failed (exit {}), see {}".format(tag, r.returncode, os.path.join(OUT, log_name)))


def metrics(tag):
    """transmission, rms vertical angle [mrad] and rms height [mm] at the asymptotic state of the transmitted bunch."""
    s = jload(os.path.join(OUT, "bunch_{}.json".format(tag)))
    d = np.load(os.path.join(OUT, "bunch_{}.npz".format(tag)), allow_pickle=True)
    st = d["asym_state"]
    ok = d["crossed"] & np.all(np.isfinite(st), axis=1)
    if "is_tail" in d.files:
        ok &= ~d["is_tail"]
    v = st[ok, 3:6]
    zp = v[:, 2] / np.hypot(v[:, 0], v[:, 1])
    return {"transmission": s["transmission"], "zp_rms_mrad": 1e3 * float(zp.std()), "zp_mean_mrad": 1e3 * float(zp.mean()),
            "z_asym_rms_mm": 1e3 * float(st[ok, 2].std()), "z_asym_mean_mm": 1e3 * float(st[ok, 2].mean()),
            "n": int(ok.sum()), "tail": s.get("tail")}


# ---------------------------------------------------------------- the grid with space charge
rows = []
for dq1, dq2 in itertools.product(a.dq1, a.dq2):
    q1, q2 = q1_0 + dq1, q2_0 + dq2
    tag = "scq_{:+.0f}_{:+.0f}".format(q1, q2)
    if not os.path.exists(os.path.join(OUT, "bunch_{}.json".format(tag))):
        t0 = time.time()
        bunch(tag, ["--n", a.n, "--sc", "--h", a.h, "--resolve-every", a.resolve_every, "--superpose", q1, 0, q2, 0, "--no-plot"],
              "log_{}.txt".format(tag))
        stamp("  q1 {:+.0f} q2 {:+.0f}: {:.0f} s".format(q1, q2, time.time() - t0))
    m = metrics(tag)
    m.update(q1=q1, q2=q2, tag=tag)
    rows.append(m)
    stamp("  q1 {:+.0f} q2 {:+.0f}: T {:.1f} %, vertical angle {:.1f} mrad rms, z rms {:.1f} mm".format(q1, q2, 100 * m["transmission"], m["zp_rms_mrad"], m["z_asym_rms_mm"]))

t_max = max(r["transmission"] for r in rows)
cands = [r for r in rows if r["transmission"] >= t_max - a.tol]
best = min(cands, key=lambda r: (r["zp_rms_mrad"], -r["transmission"]))
vac_pt = next((r for r in rows if r["q1"] == q1_0 and r["q2"] == q2_0), None)
stamp("best with space charge: q1 {:+.0f} / q2 {:+.0f} V, {:.1f} %, {:.1f} mrad (vacuum choice with SC: {})".format(
    best["q1"], best["q2"], 100 * best["transmission"], best["zp_rms_mrad"],
    "-" if vac_pt is None else "{:.1f} %, {:.1f} mrad".format(100 * vac_pt["transmission"], vac_pt["zp_rms_mrad"])))
with open(os.path.join(OUT, "sc_grid.json"), "w") as fh:
    json.dump({"rows": rows, "best": best, "vacuum_choice": {"q1": q1_0, "q2": q2_0}, "tol": a.tol}, fh, indent=2)

# ---------------------------------------------------------------- the final runs at the SC-optimal setting
final = {}
if not a.skip_final:
    for tag, extra in (("sc_all_scretune", ["--sc", "--h", a.final_h, "--resolve-every", a.final_resolve_every, "--sc-min-particles", 200]),
                       ("vac_all_scretune", [])):
        if not os.path.exists(os.path.join(OUT, "bunch_{}.json".format(tag))):
            stamp("final run {} at q1 {:+.0f} / q2 {:+.0f}: every particle (core + injected tail)".format(tag, best["q1"], best["q2"]))
            t0 = time.time()
            bunch(tag, ["--n", 10 ** 7, "--stragglers", "--superpose", best["q1"], 0, best["q2"], 0,
                        "--save-openpmd", os.path.join(OUT, "handoff_{}.h5".format(tag))] + extra, "log_{}.txt".format(tag))
            stamp("  done in {:.0f} min".format((time.time() - t0) / 60))
        final[tag] = metrics(tag)
        final[tag]["summary"] = jload(os.path.join(OUT, "bunch_{}.json".format(tag)))

# ---------------------------------------------------------------- report
ref = {}
for tag in ("vac_all", "sc_all"):
    fn = os.path.join(BASIS, "bunch_{}.json".format(tag))
    if os.path.exists(fn):
        s = jload(fn)
        d = np.load(os.path.join(BASIS, "bunch_{}.npz".format(tag)), allow_pickle=True)
        st = d["asym_state"]
        ok = d["crossed"] & np.all(np.isfinite(st), axis=1) & ~d["is_tail"]
        v = st[ok, 3:6]
        ref[tag] = {"transmission": s["transmission"], "zp_rms_mrad": 1e3 * float((v[:, 2] / np.hypot(v[:, 0], v[:, 1])).std()),
                    "tail": s.get("tail")}
L = ["# Space-charge quad retune of {} (vertical ranking)".format(os.path.basename(RUN)), "",
     "Grid around the vacuum choice q1 {:+.0f} / q2 {:+.0f} V, {} core particles per point, PyAMG space charge {:g} mA at {:.0f} mm cells, "
     "superposed basis fields of the run; best = smallest rms vertical angle (vz/v_long at 55 mm past the exit) within {:.0f} points of the "
     "best transmission.".format(q1_0, q2_0, a.n, a.current_ma, 1e3 * a.h, 100 * a.tol), "",
     "| q1 [V] | q2 [V] | transmission (core) | vertical angle [mrad rms] | z rms at 55 mm [mm] |", "|---|---|---|---|---|"]
for r in sorted(rows, key=lambda r: (r["q1"], r["q2"])):
    mark = " **(SC best)**" if r is best else (" (vacuum choice)" if r is vac_pt else "")
    L.append("| {:+.0f} | {:+.0f} | {:.1f} % | {:.1f} | {:.1f} |{}".format(r["q1"], r["q2"], 100 * r["transmission"], r["zp_rms_mrad"], r["z_asym_rms_mm"], mark))
if final:
    L += ["", "## Every particle of the RFQ file at the SC-optimal setting", "",
          "| run | transmitted | core | tail | vertical angle, core [mrad rms] |", "|---|---|---|---|---|"]
    for tag, label in (("vac_all", "final push, vacuum choice, vacuum"), ("sc_all", "final push, vacuum choice, {:g} mA".format(a.current_ma)),
                       ("vac_all_scretune", "SC retune choice, vacuum"), ("sc_all_scretune", "SC retune choice, {:g} mA".format(a.current_ma))):
        m = final.get(tag) or ref.get(tag)
        if not m:
            continue
        t = m.get("tail") or {}
        L.append("| {} | {:.2f} % | {} | {} | {:.1f} |".format(
            label, 100 * m["transmission"],
            "{:.2f} %".format(100 * t["core_transmission"]) if t else "-", "{:.2f} %".format(100 * t["tail_transmission"]) if t else "-",
            m["zp_rms_mrad"]))
    L += ["", "Hand-off files: `sc_retune/handoff_sc_all_scretune.h5`, `sc_retune/handoff_vac_all_scretune.h5` (machine frame, at the electrode exit)."]
with open(os.path.join(RUN, "sc_retune.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
stamp("report -> {}".format(os.path.join(RUN, "sc_retune.md")))
