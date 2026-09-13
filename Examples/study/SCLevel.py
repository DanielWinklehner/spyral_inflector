"""Space charge in the loop: level the BUNCH centroid -- its mean height and mean vertical
angle at the asymptotic state (55 mm past the exit) -- with the EXIT TRUNCATION and the AXIAL
SHIFT as knobs, the spiral voltage staying the gap-centering knob (the inner solve of every
rebuild), on a run whose basis fields exist.

Each evaluation is one space-charge bunch run of --n core particles at the quad setting of
--best-json. A truncation change rebuilds the geometry (design-orbit build with both
truncations and dz fixed: only the inner voltage/centering solve runs), re-solves the spiral
basis field of the new electrodes (the quad basis fields are reused) and tracks; a dz change
is an exact translation of the assembly's field, design orbit and collision geometry at
tracking time (bunch --shift-z), no solve. Newton iterations on a measured 2x2 Jacobian
(finite differences --dt [deg], --ddz [m]), damped to --max-dt / --max-ddz per step, until
|height| < --tol-z and |angle| < --tol-zp. With --final the levelled setting is run with every
RFQ particle (core + injected tail) with space charge at --final-h cells and once in vacuum,
both with hand-off files.

    python SCLevel.py --run <deck>\\Results\\triplet\\triplet1 --final
    python SCLevel.py --run <deck>\\Results\\final\\vfocus1_final --best-json basis\\scan2_fine.json --quads 16 24 --final

Writes <run>/sc_level/ (geometries, basis fields and bunch runs per evaluation, the full runs,
hand-offs) and <run>/sc_level.md.
"""
import argparse
import json
import os
import pickle
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
p.add_argument("--run", required=True, help="run folder with geometry/summary.json, exit_plane.json, basis/ and steps/")
p.add_argument("--best-json", default=None, help="scan json with the quad setting (default: basis/scan2_scan.json, else scan2_fine.json)")
p.add_argument("--quads", type=float, nargs="+", default=None, help="basis rotation of the quads [deg] (default: from the run's spiral state)")
p.add_argument("--base-run", default=os.path.join(DECK, "Results", "final", "final1"), help="the rotated field the run's geometry was optimized in")
p.add_argument("--n", type=int, default=5000)
p.add_argument("--h", type=float, default=0.004)
p.add_argument("--resolve-every", type=int, default=16)
p.add_argument("--current-ma", type=float, default=8.0)
p.add_argument("--dt", type=float, default=0.3, help="finite-difference step of the exit truncation [deg]")
p.add_argument("--ddz", type=float, default=0.001, help="finite-difference step of the axial shift [m]")
p.add_argument("--max-dt", type=float, default=1.0)
p.add_argument("--max-ddz", type=float, default=0.003)
p.add_argument("--iters", type=int, default=3)
p.add_argument("--tol-z", type=float, default=0.3e-3, help="[m]")
p.add_argument("--tol-zp", type=float, default=2.0e-3, help="[rad]")
p.add_argument("--res", type=float, default=None, help="resolution of the re-solved spiral basis [m]; default: the run's basis grid (and its box)")
p.add_argument("--geo-h", type=float, default=0.005)
p.add_argument("--final", action="store_true")
p.add_argument("--final-h", type=float, default=0.002)
p.add_argument("--final-resolve-every", type=int, default=8)
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle"))
a = p.parse_args()

RUN = os.path.abspath(a.run)
BASIS0, STEPS0 = os.path.join(RUN, "basis"), os.path.join(RUN, "steps")
OUT = os.path.join(RUN, "sc_level")
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
        raise SystemExit("{} failed (exit {}), see {}".format(cmd[0], r.returncode, log_path))


from spyral_inflector.tracking.deck import particle_rows  # noqa: E402
from spyral_inflector.tracking.geometry_point import build_geometry  # noqa: E402

geo0 = jload(os.path.join(RUN, "geometry", "summary.json"))
knobs = dict(geo0["knobs"])
opt0 = geo0["optimizer"]
t_ent, t_exit0 = (opt0.get("truncations_deg") or [0.34, 0.77])
dz0 = 1e-3 * opt0["dz_mm"]
R = jload(os.path.join(RUN, "exit_plane.json"))["rotation_deg"]
PHI = 90.0 + R
ROT_BF = os.path.join(a.base_run, "bfield_rotated.pickle")
ENERGY = float(particle_rows(a.particles, core=True)[:, 8].mean())
with open(os.path.join(BASIS0, "si_state_spiral.pickle"), "rb") as fh:
    st_spiral = pickle.load(fh)
QROT = a.quads or [float(x) for x in (st_spiral.get("quad_rotation_deg") or [16.0, 24.0])]
n_quads = int(knobs.get("n_quads", 2) or 2)
QROT = (QROT + [0.0] * n_quads)[:n_quads]
bj = a.best_json or next((f for f in ("scan2_scan.json", "scan2_fine.json") if os.path.exists(os.path.join(BASIS0, f))), None)
if bj is None:
    raise SystemExit("no scan json with the quad setting found in " + BASIS0)
best = jload(bj if os.path.isabs(bj) else os.path.join(RUN, bj) if os.path.exists(os.path.join(RUN, bj)) else os.path.join(BASIS0, bj))["best"]
quads = [best["q1"], best.get("alpha1", 0.0), best["q2"], best.get("alpha2", 0.0)]
if best.get("q3") is not None:
    quads += [best["q3"], best.get("alpha3", 0.0)]
# the re-solved spiral basis must sit on the run's basis grid (superposition needs identical grids)
from PyPATools.field import Field  # noqa: E402
_g = Field.from_file(os.path.join(BASIS0, "ef_itp_spiral.pickle")).grid
GRID_RES = a.res or float(np.asarray(_g["x"], dtype=float)[1] - np.asarray(_g["x"], dtype=float)[0])
GRID_BOX = [float(np.asarray(_g[c], dtype=float)[i]) for c in "xyz" for i in (0, -1)]
quad_basis_files = ["ef_itp_q{}.pickle".format(i + 1) for i in range(len(quads) // 2)] + \
                   ["ef_itp_q{}skew.pickle".format(i + 1) for i in range(len(quads) // 2) if quads[2 * i + 1] != 0.0]
stamp("SC LEVEL of {}: quads {} (from {}), basis rotation {}, R {:+.3f} deg; start exit truncation {:.3f} deg, dz {:+.2f} mm; "
      "{} mA, {} particles, {:.0f} mm cells; basis grid {:.1f} mm".format(
          os.path.basename(RUN), " ".join("{:+.0f}@{:g}".format(quads[2 * i], quads[2 * i + 1]) for i in range(len(quads) // 2)),
          os.path.basename(bj), QROT, R, t_exit0, 1e3 * dz0, a.current_ma, a.n, 1e3 * a.h, 1e3 * GRID_RES))


def geometry_for(t_exit):
    """The run's geometry rebuilt with the exit truncation t_exit (entrance and dz as in the run,
    the voltage from the inner centering solve) and its spiral basis field; the quad basis
    fields of the run are reused. Returns the basis dir."""
    key = "t{:.3f}".format(t_exit)
    gdir, sdir, bdir = (os.path.join(OUT, key, d) for d in ("geometry", "steps", "basis"))
    if abs(t_exit - t_exit0) < 1e-9:
        gdir, sdir, bdir = os.path.join(RUN, "geometry"), STEPS0, BASIS0
        return bdir
    if not os.path.exists(os.path.join(gdir, "summary.json")):
        if os.path.isdir(gdir):
            shutil.rmtree(gdir)
        stamp("  geometry at exit truncation {:.3f} deg (dz {:+.2f} mm fixed, voltage from the centering solve)".format(t_exit, 1e3 * dz0))
        build_geometry(gdir, sdir, ROT_BF, ENERGY, knobs=knobs, fix_truncations=(t_ent, t_exit), fix_dz=dz0, maxiter=0,
                       res=0.005, h=a.geo_h, log=stamp)
    os.makedirs(bdir, exist_ok=True)
    if not os.path.exists(os.path.join(bdir, "ef_itp_spiral.pickle")):
        stamp("  spiral basis field of that geometry")
        run_cmd(["-m", "spyral_inflector.tracking.bem_reload", "--tag", "spiral", "--step-dir", sdir, "--voltages", os.path.join(gdir, "voltages.csv"),
                 "--state", os.path.join(gdir, "state.pickle"), "--out-dir", bdir, "--res", GRID_RES, "--box"] + GRID_BOX +
                ["--h", a.geo_h, "--bfield", a.bfield, "--rotate-all", R, "--rotate-quads"] + QROT + ["--quad-voltages"] + [0] * n_quads +
                ["--energy-mev", ENERGY],
                os.path.join(bdir, "log_basis_spiral.txt"))
    for fn in quad_basis_files:
        if not os.path.exists(os.path.join(bdir, fn)):
            shutil.copyfile(os.path.join(BASIS0, fn), os.path.join(bdir, fn))
    return bdir


def bunch(bdir, sdir, tag, extra):
    cmd = ["-m", "spyral_inflector.tracking.bunch", "--tag", "spiral", "--reload-dir", bdir, "--step-dir", sdir,
           "--particles", a.particles, "--bfield", a.bfield, "--phi", PHI, "--out-dir", OUT, "--out-tag", tag,
           "--basis-dir", bdir, "--superpose"] + quads + ["--current-ma", a.current_ma, "--handoff-frame", "machine",
           "--handoff-distance", 0.0] + list(extra)
    run_cmd(cmd, os.path.join(OUT, "log_{}.txt".format(tag)))


def measure(tag):
    s = jload(os.path.join(OUT, "bunch_{}.json".format(tag)))
    d = np.load(os.path.join(OUT, "bunch_{}.npz".format(tag)), allow_pickle=True)
    st = d["asym_state"]
    ok = d["crossed"] & np.all(np.isfinite(st), axis=1)
    if "is_tail" in d.files:
        ok &= ~d["is_tail"]
    z = st[ok, 2]
    v = st[ok, 3:6]
    zp = v[:, 2] / np.hypot(v[:, 0], v[:, 1])
    return {"z": float(z.mean()), "zp": float(zp.mean()), "z_rms": float(z.std()), "zp_rms": float(zp.std()),
            "transmission": s["transmission"], "n": int(ok.sum()), "tail": s.get("tail"),
            "voltage_V": s.get("voltages", {}).get("SI_Anode")}


def evaluate(t_exit, dz, label):
    bdir = geometry_for(t_exit)
    sdir = STEPS0 if bdir == BASIS0 else os.path.join(os.path.dirname(bdir), "steps")
    tag = "lvl_{}".format(label)
    if not os.path.exists(os.path.join(OUT, "bunch_{}.json".format(tag))):
        t0 = time.time()
        bunch(bdir, sdir, tag, ["--n", a.n, "--sc", "--h", a.h, "--resolve-every", a.resolve_every, "--shift-z", dz - dz0, "--no-plot"])
        stamp("  {}: {:.0f} s".format(tag, time.time() - t0))
    m = measure(tag)
    m.update(t_exit=t_exit, dz=dz, tag=tag, basis=bdir)
    stamp("  truncation {:.3f} deg, dz {:+.2f} mm -> height {:+.2f} mm, angle {:+.1f} mrad (rms {:.1f} mm / {:.1f} mrad), T {:.1f} %, V {}".format(
        t_exit, 1e3 * dz, 1e3 * m["z"], 1e3 * m["zp"], 1e3 * m["z_rms"], 1e3 * m["zp_rms"], 100 * m["transmission"],
        "-" if m["voltage_V"] is None else "{:.0f} V".format(m["voltage_V"])))
    return m


# ---------------------------------------------------------------- Newton on the measured Jacobian
history = []
t_exit, dz = t_exit0, dz0
cur = evaluate(t_exit, dz, "0")
history.append(cur)
converged = abs(cur["z"]) < a.tol_z and abs(cur["zp"]) < a.tol_zp
for it in range(1, a.iters + 1):
    if converged:
        break
    stamp("iteration {}: Jacobian".format(it))
    mt = evaluate(min(t_exit + a.dt, 15.0), dz, "{}t".format(it))
    mz = evaluate(t_exit, dz + a.ddz, "{}z".format(it))
    dt_eff = mt["t_exit"] - t_exit
    J = np.array([[(mt["z"] - cur["z"]) / dt_eff, (mz["z"] - cur["z"]) / a.ddz],
                  [(mt["zp"] - cur["zp"]) / dt_eff, (mz["zp"] - cur["zp"]) / a.ddz]])
    stamp("  J = [[dz/dt {:+.2e} m/deg, dz/ddz {:+.3f}], [dzp/dt {:+.2e} rad/deg, dzp/ddz {:+.3f} rad/m]]".format(J[0, 0], J[0, 1], J[1, 0], J[1, 1]))
    try:
        step = -np.linalg.solve(J, np.array([cur["z"], cur["zp"]]))
    except np.linalg.LinAlgError:
        stamp("  singular Jacobian; stopping")
        break
    scale = min(1.0, a.max_dt / max(abs(step[0]), 1e-12), a.max_ddz / max(abs(step[1]), 1e-12))
    step *= scale
    t_exit = float(np.clip(t_exit + step[0], 0.0, 15.0))
    dz = dz + float(step[1])
    stamp("  step: dt {:+.3f} deg, ddz {:+.2f} mm (damping {:.2f}) -> truncation {:.3f} deg, dz {:+.2f} mm".format(step[0], 1e3 * step[1], scale, t_exit, 1e3 * dz))
    cur = evaluate(t_exit, dz, str(it))
    history.append(cur)
    converged = abs(cur["z"]) < a.tol_z and abs(cur["zp"]) < a.tol_zp
stamp("levelled: exit truncation {:.3f} deg, dz {:+.2f} mm -> height {:+.2f} mm, angle {:+.1f} mrad ({})".format(
    t_exit, 1e3 * dz, 1e3 * cur["z"], 1e3 * cur["zp"], "converged" if converged else "not converged"))
with open(os.path.join(OUT, "sc_level.json"), "w") as fh:
    json.dump({"quads": quads, "quad_rotation_deg": QROT, "history": history, "t_exit_deg": t_exit, "dz_m": dz, "dz0_m": dz0,
               "shift_z_m": dz - dz0, "converged": converged, "basis": cur["basis"]}, fh, indent=2)

# ---------------------------------------------------------------- final runs at the levelled setting
final = {}
if a.final:
    bdir = cur["basis"]
    sdir = STEPS0 if bdir == BASIS0 else os.path.join(os.path.dirname(bdir), "steps")
    for tag, extra in (("sc_all_level", ["--sc", "--h", a.final_h, "--resolve-every", a.final_resolve_every, "--sc-min-particles", 200]),
                       ("vac_all_level", [])):
        if not os.path.exists(os.path.join(OUT, "bunch_{}.json".format(tag))):
            stamp("final run {}: every particle (core + injected tail) at truncation {:.3f} deg, dz {:+.2f} mm".format(tag, t_exit, 1e3 * dz))
            t0 = time.time()
            bunch(bdir, sdir, tag, ["--n", 10 ** 7, "--stragglers", "--shift-z", dz - dz0,
                                    "--save-openpmd", os.path.join(OUT, "handoff_{}.h5".format(tag))] + extra)
            stamp("  done in {:.0f} min".format((time.time() - t0) / 60))
        final[tag] = measure(tag)

# ---------------------------------------------------------------- report
L = ["# Space-charge levelling of {}".format(os.path.basename(RUN)), "",
     "Quads {} at basis rotation {} deg (from {}); {:g} mA, {} core particles per evaluation at {:.0f} mm cells. Knobs: exit "
     "truncation (geometry rebuilt, spiral voltage from the centering solve, spiral basis re-solved, quad bases reused) and axial "
     "shift of the whole assembly (exact field translation). Targets: centroid height and mean vertical angle vz/v_long of the "
     "transmitted core at 55 mm past the exit. Start: the run's design orbit (truncation {:.3f} deg, dz {:+.2f} mm).".format(
         " ".join("{:+.0f} V @ {:g} deg".format(quads[2 * i], quads[2 * i + 1]) for i in range(len(quads) // 2)), QROT, os.path.basename(bj),
         a.current_ma, a.n, 1e3 * a.h, t_exit0, 1e3 * dz0), "",
     "| evaluation | exit trunc [deg] | dz [mm] | V [V] | height [mm] | angle [mrad] | rms height | rms angle | transmission |",
     "|---|---|---|---|---|---|---|---|---|"]
for m in history:
    L.append("| {} | {:.3f} | {:+.2f} | {} | {:+.2f} | {:+.1f} | {:.1f} | {:.1f} | {:.1f} % |".format(
        m["tag"], m["t_exit"], 1e3 * m["dz"], "-" if m["voltage_V"] is None else "{:.0f}".format(m["voltage_V"]), 1e3 * m["z"], 1e3 * m["zp"],
        1e3 * m["z_rms"], 1e3 * m["zp_rms"], 100 * m["transmission"]))
L += ["", "**Levelled: exit truncation {:.3f} deg, dz {:+.2f} mm (shift {:+.2f} mm from the run); height {:+.2f} mm, angle {:+.1f} mrad ({}).**".format(
    t_exit, 1e3 * dz, 1e3 * (dz - dz0), 1e3 * cur["z"], 1e3 * cur["zp"], "converged" if converged else "not converged")]
if final:
    L += ["", "## Every particle at the levelled setting", "", "| run | transmitted | core | tail | height [mm] | angle [mrad] | rms angle |", "|---|---|---|---|---|---|---|"]
    for tag, label in (("vac_all_level", "vacuum"), ("sc_all_level", "{:g} mA".format(a.current_ma))):
        m = final.get(tag)
        if not m:
            continue
        t = m.get("tail") or {}
        L.append("| {} | {:.2f} % | {} | {} | {:+.2f} | {:+.1f} | {:.1f} |".format(
            label, 100 * m["transmission"], "{:.2f} %".format(100 * t["core_transmission"]) if t else "-",
            "{:.2f} %".format(100 * t["tail_transmission"]) if t else "-", 1e3 * m["z"], 1e3 * m["zp"], 1e3 * m["zp_rms"]))
    L += ["", "Hand-off files: `sc_level/handoff_sc_all_level.h5`, `sc_level/handoff_vac_all_level.h5` (machine frame, electrode exit).",
          "The levelled geometry: the STEP files of `sc_level/{}/steps` (exit truncation {:.3f} deg) shifted by {:+.2f} mm along z; "
          "the exit plate moves with it, so the plate-rule rotation should be re-solved before a Baseline export.".format(
              "t{:.3f}".format(t_exit) if abs(t_exit - t_exit0) > 1e-9 else "..", t_exit, 1e3 * (dz - dz0))]
with open(os.path.join(RUN, "sc_level.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
stamp("report -> {}".format(os.path.join(RUN, "sc_level.md")))
