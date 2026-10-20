"""A/B of a geometry change -- a RECTANGULAR first aperture (the vacuum-sleeve stand-in, ent0), plate holes, or any
knob via --knob (2026-10-20: the round housing exit opening, exit_circle) -- against a levelled reference run.

Daniel, 2026-10-16: "When you have the rectangle figured out, run a single full run with all particles and 10 mA as an
A/B and show me the results before doing any optimization runs."

Everything is the reference's LEVELLED setting except the sleeve: the reference run's knobs with entrance_rect added,
both truncations fixed at the levelled values and dz fixed at the run's design-orbit dz (the voltage from the same
centering solve the levelling used), the same system rotation R, the same quads (voltages AND turns), the same voltage
scale, levelling shift, beam rotation, particles and space-charge settings. Bases re-solved as the reference had them:
the spiral basis without the quad poles (SCLevel's default), the quad bases and their skews with every electrode
(TripletScan's), all on this geometry -- the rectangle changes quad 1's fringe, so the quad bases must be new.
No optimization and no levelling: one every-particle bunch run at the reference's setting.

    python ABRect.py --ref-run E:\\HCHC1.5\\Results\\triplet\\t1nc --name t1nc_rect --entrance-rect 0.017 0.0065 9.4926
    python ABRect.py --ref-run E:\\HCHC1.5\\Results\\triplet\\t1nr_c2 --name t1nr_c2_circ --base-run <ref> --knob exit_circle=[0.01,-0.000956,-0.001214]
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
PY = sys.executable
sys.path.insert(0, SCRIPTS)
p = argparse.ArgumentParser()
p.add_argument("--ref-run", required=True)
p.add_argument("--name", required=True)
p.add_argument("--entrance-rect", type=float, nargs=3, default=None, metavar=("HALF_U", "HALF_W", "ANGLE_DEG"),
               help="half-widths [m] and the turn of the rectangle's long (u) axis from the deck x axis [deg], before R "
                    "(default: the reference's)")
p.add_argument("--knob", action="append", default=[], metavar="KEY=JSON",
               help="set any geometry knob, e.g. exit_circle=[0.01,-0.000956,-0.001214]; repeatable")
p.add_argument("--base-run", default=r"E:\HCHC1.5\Results\final\d1", help="folder with bfield_rotated.pickle (the design-orbit field)")
p.add_argument("--bfield", default=r"E:\HCHC1.5\Fields\hchc15_injax_1mm.pickle")
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--res", type=float, default=0.0025)
p.add_argument("--h", type=float, default=0.005)
p.add_argument("--geo-res", type=float, default=0.005)
p.add_argument("--record", type=int, default=48843)
p.add_argument("--exit-rect", type=float, nargs=4, action="append", default=None, metavar=("K", "HALF_U", "HALF_W", "ANGLE_DEG"),
               help="rectangular hole in the plate AFTER quad K (1 = the plate between quads 1 and 2); repeatable")
p.add_argument("--bases-only", action="store_true", help="build the geometry and solve the bases, no bunch run (for QuadOptN)")
p.add_argument("--quads-json", default=None, help="take the quads from this QuadOptN best.json instead of the reference's")
p.add_argument("--tag", default="sc_all_ab")
a = p.parse_args()

from spyral_inflector.tracking.deck import particle_rows  # noqa: E402
from spyral_inflector.tracking.geometry_point import build_geometry  # noqa: E402

RESULTS = os.environ.get("SI_RESULTS") or os.path.join(DECK, "Results")
OUT = os.path.join(RESULTS, "triplet", a.name)
GEO, STEPS, BASIS = (os.path.join(OUT, d) for d in ("geometry", "steps", "basis"))
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def jload(path):
    with open(path) as fh:
        return json.load(fh)


def run(cmd, log_path):
    with open(log_path, "w", encoding="utf-8") as fh:
        r = subprocess.run([PY, "-u"] + [str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS)
    if r.returncode != 0:
        stamp("FAILED (exit {}): see {}".format(r.returncode, log_path))
        raise SystemExit(r.returncode)


ref = a.ref_run
s0 = jload(os.path.join(ref, "geometry", "summary.json"))
knobs = dict(s0["knobs"])
t_ent = float(s0["optimizer"]["truncations_deg"][0])
lv = jload(os.path.join(ref, "sc_level", "sc_level.json"))
t_exit, dz0, shift = float(lv["t_exit_deg"]), float(lv["dz0_m"]), float(lv["shift_z_m"])
quads, qrot, vscale = [float(x) for x in lv["quads"]], [float(x) for x in lv["quad_rotation_deg"]], float(lv["vscale"])
bj = jload(os.path.join(ref, "sc_level", "bunch_sc_all_level.json"))
R, phi = float(bj["assembly_rotation_deg"]), float(bj["phi_deg"])
if a.entrance_rect is not None:
    knobs["entrance_rect"] = list(a.entrance_rect)
for _kv in a.knob:
    _key, _, _val = _kv.partition("=")
    knobs[_key.strip()] = json.loads(_val)
if a.exit_rect:
    _n = int(knobs.get("n_quads", 3) or 3)
    _xr = [None] * _n
    for k_, hu, hw, ang in a.exit_rect:
        _xr[int(k_) - 1] = [hu, hw, ang]
    knobs["exit_rects"] = _xr
if a.quads_json:
    _b = jload(a.quads_json)["best"]
    quads = [float(_b["q1"]), float(_b.get("alpha1", 0.0)), float(_b["q2"]), float(_b.get("alpha2", 0.0)),
             float(_b["q3"]), float(_b.get("alpha3", 0.0))]
ENERGY = float(particle_rows(a.particles, core=True)[:, 8].mean())
ROT_BF = os.path.join(a.base_run, "bfield_rotated.pickle")
_changed = {k_: v_ for k_, v_ in knobs.items() if s0["knobs"].get(k_) != v_}
stamp("A/B {} vs {}: knobs changed {}; everything else the reference's levelled setting: truncations {:.4f} / {:.4f} deg, "
      "dz {:+.3f} mm + shift {:+.3f} mm, R {:+.4f} deg, quads {} at basis turns {}, vscale {}".format(
          a.name, os.path.basename(ref), _changed, t_ent, t_exit, 1e3 * dz0, 1e3 * shift, R, quads, qrot, vscale))

# 1. geometry at the reference's levelled truncations and dz (voltage from the centering solve, as SCLevel builds it)
if not os.path.exists(os.path.join(GEO, "summary.json")):
    build_geometry(GEO, STEPS, ROT_BF, ENERGY, knobs=knobs, fix_truncations=(t_ent, t_exit), fix_dz=dz0, maxiter=0,
                   res=a.geo_res, h=a.h, log=stamp)
else:
    # an existing folder is reused only if it IS this geometry (2026-10-16: a stale folder of the same name from a
    # stopped chain was silently reused with a different sleeve and no quad1-quad2 slot)
    _s = jload(os.path.join(GEO, "summary.json"))
    _k = _s.get("knobs", {})
    _same = (json.loads(json.dumps(_k)) == json.loads(json.dumps(knobs))
             and np.allclose(_s.get("fixed_truncations_deg") or [np.nan, np.nan], [t_ent, t_exit])
             and _s.get("fix_dz_m") is not None and abs(float(_s["fix_dz_m"]) - dz0) < 1e-9)
    if not _same:
        raise SystemExit("{} holds a DIFFERENT geometry (knobs differing {}, truncations {}, fix_dz {}); use another "
                         "--name".format(GEO, {k_: (_k.get(k_), v_) for k_, v_ in knobs.items() if _k.get(k_) != v_},
                                         _s.get("fixed_truncations_deg"), _s.get("fix_dz_m")))
    stamp("  reusing the geometry already in {} (same knobs, truncations and dz)".format(GEO))
volt_ab = jload(os.path.join(GEO, "summary.json"))["voltage"]
ref_geo = os.path.join(ref, "sc_level", "t{:.3f}".format(t_exit), "geometry", "summary.json")
volt_ref = jload(ref_geo)["voltage"] if os.path.exists(ref_geo) else float("nan")
stamp("  spiral voltage {:.1f} V (reference levelled geometry {:.1f} V)".format(volt_ab, volt_ref))

# 2. bases: spiral without the quad poles (SCLevel), quads + skews with every electrode (TripletScan)
os.makedirs(BASIS, exist_ok=True)
a1, a2, a3 = qrot
common = ["--step-dir", STEPS, "--voltages", os.path.join(GEO, "voltages.csv"), "--state", os.path.join(GEO, "state.pickle"),
          "--out-dir", BASIS, "--res", a.res, "--h", a.h, "--bfield", a.bfield, "--rotate-all", R, "--energy-mev", ENERGY]
bases = (("spiral", ["--quad-voltages", 0, 0, 0, "--rotate-quads", a1, a2, a3, "--drop", "D*"]),
         ("q1", ["--spiral-voltage", 0, "--quad-voltages", 3500, 0, 0, "--rotate-quads", a1, a2, a3, "--no-test"]),
         ("q1skew", ["--spiral-voltage", 0, "--quad-voltages", 3500, 0, 0, "--rotate-quads", a1 + 45, a2, a3, "--no-test"]),
         ("q2", ["--spiral-voltage", 0, "--quad-voltages", 0, 3500, 0, "--rotate-quads", a1, a2, a3, "--no-test"]),
         ("q2skew", ["--spiral-voltage", 0, "--quad-voltages", 0, 3500, 0, "--rotate-quads", a1, a2 + 45, a3, "--no-test"]),
         ("q3", ["--spiral-voltage", 0, "--quad-voltages", 0, 0, 3500, "--rotate-quads", a1, a2, a3, "--no-test"]),
         ("q3skew", ["--spiral-voltage", 0, "--quad-voltages", 0, 0, 3500, "--rotate-quads", a1, a2, a3 + 45, "--no-test"]))
for tag, extra in bases:
    if os.path.exists(os.path.join(BASIS, "ef_itp_{}.pickle".format(tag))):
        continue
    t0 = time.time()
    run(["-m", "spyral_inflector.tracking.bem_reload", "--tag", tag] + extra + common, os.path.join(OUT, "log_basis_{}.txt".format(tag)))
    stamp("  basis {} ({:.0f} s)".format(tag, time.time() - t0))

# QuadOptN reads the system rotation from the run folder: the same R as the reference (the geometry is turned by it)
_ep = os.path.join(OUT, "exit_plane.json")
if not os.path.exists(_ep):
    with open(_ep, "w") as fh:
        json.dump(dict(jload(os.path.join(ref, "exit_plane.json")), note="copied from the reference run by ABRect.py"), fh, indent=1)
if a.bases_only:
    stamp("bases only: geometry and the seven bases are in {}".format(OUT))
    raise SystemExit(0)

# 3. one every-particle run at 10 mA, the reference's setting (or the --quads-json quads)
tag = a.tag
if not os.path.exists(os.path.join(OUT, "bunch_{}.json".format(tag))):
    stamp("every particle at {:g} mA, the reference's quads and shift".format(bj["sc"]["current_mA"]))
    t0 = time.time()
    run(["-m", "spyral_inflector.tracking.bunch", "--tag", "spiral", "--reload-dir", BASIS, "--step-dir", STEPS, "--particles", a.particles,
         "--bfield", a.bfield, "--phi", "{:.10f}".format(phi), "--out-dir", OUT, "--out-tag", tag, "--basis-dir", BASIS, "--superpose"] +
        ["{:.10f}".format(x) for x in quads] +
        ["--current-ma", bj["sc"]["current_mA"], "--vscale", "{:.6f}".format(vscale), "--handoff-frame", "machine", "--handoff-distance", 0.0,
         "--n", 10 ** 7, "--stragglers", "--shift-z={:.10f}".format(shift), "--sc", "--h", bj["sc"]["h_m"], "--resolve-every", 8,
         "--sc-min-particles", 200, "--record", a.record], os.path.join(OUT, "log_{}.txt".format(tag)))
    stamp("  done in {:.0f} min".format((time.time() - t0) / 60))
r = jload(os.path.join(OUT, "bunch_{}.json".format(tag)))
stamp("A/B result: transmitted {:.2f} % (reference {:.2f} %)".format(100 * r["transmission"], 100 * bj["transmission"]))
