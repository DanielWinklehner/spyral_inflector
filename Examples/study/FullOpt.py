"""Full DFO-LS optimization of the doublet or triplet system with space charge in the loop, seeded at a
final-push / TripletScan / levelled run's setting. Knobs (10 for a doublet, 11 for a triplet): k' tilt, angling, gamma,
V-depth sigma, entrance and exit truncation, axial shift dz, quad voltages q1..qN (the run says N), spiral gap,
and vscale = the operating spiral voltage as a multiple of the gap-centering value (the 2026-09-14 laptop result: 1.03
buys exit transmission; the angle term and the aperture proxy --w-z100 see what it costs downstream) (the design voltage scales with the gap,
the operating voltage stays the inner gap-centering solve). Each evaluation rebuilds the geometry at
the run's rotation R (design orbit with both truncations and dz fixed: only the centering solve runs),
re-solves the spiral basis field on the run's basis grid (the quad bases are reused), tracks --n
particles with PyAMG space charge at the quad setting by superposition, and returns the residual vector

    [ rms vertical angle / 3 mrad,  max(0, T_floor - T) / 1 point,  centroid height / 0.3 mm,  centroid angle / 2 mrad ]

at 55 mm past the exit (T_floor = the seed's transmission - 5 points unless --t-floor is given).
DFO-LS (noise-aware, bound-scaled) minimizes the sum of squares. Every evaluation is cached in
<run>/full_opt/evals.json (resumable); the big basis files of each evaluation are deleted after the
measurement. Writes <run>/full_opt.md and <run>/full_opt/best_knobs.json (knobs for RunFinalPush
--knobs-json plus _truncations_deg, _dz_m, _q1, _q2).

2026-10-21 (HCHC-1.5, after Barnard et al., PRAB 24, 023501 (2021), whose inflector C3 shaped the two field-gradient
profiles Q1(s), Q2(s) piecewise linearly on 3 nodes each and cut the rf phase spread from +-39 to +-22 deg):
  --profiles angling sigma   frees angling (Q2-type, the angled electrode faces) and/or the V-depth sigma (Q1-type) as
                             3-node PROFILES along the electrode, entrance / middle / exit of the truncated electrode,
                             linearly interpolated (knobs ang_ent ang_mid ang_exit [deg], sig_ent_mm sig_mid_mm sig_exit_mm
                             [mm]; geometry knobs angling_profile / sigma_profile). The seed is the scalar's legacy meaning:
                             angling a -> [a, 0, -a] (the old +a -> -a ramp), sigma s -> [s, s, s]. Without --profiles the
                             scalars work as before and old caches stay valid.
  --free-aspect              electrode width / gap as a knob (2026-10-16: the anode-cathode gap closes between the electrode
                             EDGES, so width, not angling, is the lever under the 8 mm rule --min-gap 8).
  whole-file loop (default)  the loop draws --n particles at random from the WHOLE RFQ file, core and tail (bunch.py
                             --stragglers, the tail injected at its arrival time), and every term is over all of them
                             (Daniel: no core-only numbers anywhere). --core-only restores the old core loop.
  --w-dt                     the rms RF phase at the hand-off plane, WRAPPED to +-180 deg around the core's median (a tail
                             particle N periods late sees the same phase; beam_report.py's definition) / 10 deg.
  --w-loss, --w-hv           (2026-10-24) QuadOptN's loss terms: w * 100 * (1 - T) and w * points lost on HV electrodes, so a
                             shape round sees transmission (sp1's objective did not).
  --free-height              (2026-11-10) the electric radius A of the design orbit (= the inflector height; A = T / (q E / 2),
                             E = 2 V / gap) as a knob height_mm: the design voltage becomes T gap / A. Without it V/gap is the run's
                             and A is fixed. With B and T fixed, the hand-off radius AND the exit direction need A and k' together
                             (v6: k' alone trades 71 mm / +14.9 deg for 62 mm / +6.6 deg).
                             CAUTION: the builder runs the LAST quad down to the entrance plate on the inflector top (z = -A), so A
                             changes Q3's length (v6r60_c3: 25 -> 38 mm at A 106 -> 93) while every evaluation reuses the RUN's quad
                             bases. Only consistent for small A steps; to move A, build a seed at the new A (level stage) and hold it.
  --r-target-mm, --w-r       (2026-11-10) w * (r - target) / 1 mm, r = |mean position| of the transmitted beam at the hand-off plane.
  --alpha-target-deg, --w-alpha  w * (alpha - target) / 1 deg, alpha = angle of the MEAN momentum at the hand-off from the local
                             tangent (direction of circulation), + = outward. Both always measured and logged; the central
                             region wants the beam turned outward (Daniel, 2026-11-10: ~60 mm and ~7 deg more than v6t40o_c3's 7.1).
  seed quads                 from --best-json if given, else from <run>/sc_level/sc_level.json (voltages AND the skew turns
                             alpha, which stay fixed; their skew bases are copied with the normal ones), else
                             basis/scan2_coarse.json; --quads (basis rotations) default to the levelling's.
  --hold NAME..., --set NAME=VALUE, --bound NAME LO HI, --seed-only, --dry-run: see --help.

    python FullOpt.py --run <deck>\\Results\\final\\vfocus1_level --best-json basis\\scan2_coarse.json --maxfun 80
    python FullOpt.py --run E:\\HCHC1.5\\Results\\triplet\\t1nr_circ_c2 --out-name prof1 --profiles angling sigma --free-aspect
                      --min-gap 8 --w-dt 3 --current-ma 10 --n 4000 --maxfun 60
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
p.add_argument("--best-json", default=None,
               help="scan json whose 'best' holds the seed q1..qN (and optionally alpha1..alphaN). Default: the levelled quads of "
                    "<run>/sc_level/sc_level.json if it has them, else basis\\scan2_coarse.json")
p.add_argument("--quads", type=float, nargs="+", default=None,
               help="basis rotation of the quads on top of R [deg], one per quad (the run says how many); default: the levelling's "
                    "quad_rotation_deg, else 16 24")
p.add_argument("--n", type=int, default=8000)
p.add_argument("--h", type=float, default=0.004)
p.add_argument("--resolve-every", type=int, default=10)   # 10 since 2026-09-15 23:50 (was 16)
p.add_argument("--current-ma", type=float, default=8.0)
p.add_argument("--geo-h", type=float, default=0.005)
p.add_argument("--hold-quads", action="store_true",
               help="take the quad voltages OUT of the optimizer's knob vector and hold them at the "
                    "seed. Pinning them by bounds instead leaves them in the vector, so DFO-LS still "
                    "spends one of its n+1 model-building evaluations on each -- three ~150 s "
                    "evaluations perturbing a voltage across a 2 V bound. The full vector is still "
                    "what is evaluated, cached and reported, so caches stay interchangeable")
p.add_argument("--hold", nargs="+", default=[], metavar="NAME",
               help="take these knobs out of the optimizer's vector and hold them at the seed (generalizes --hold-quads; the full "
                    "vector is still evaluated and cached)")
p.add_argument("--geo-res", type=float, default=0.005,
               help="potential grid of the geometry build's design-orbit solve [m]. Measured "
                    "2026-09-19: 5 mm and 2.5 mm put the converged design point 0.263 deg of exit "
                    "truncation and 0.130 mm of dz apart, while BOTH report residuals inside "
                    "tolerance -- so this is an accuracy knob, not a speed one")
p.add_argument("--keep-quad-poles", dest="drop_quad_poles", action="store_false",
               help="leave the grounded quad POLES in the per-evaluation spiral solve. They are "
                    "dropped by default since 2026-09-19: measured against a solve that keeps them, "
                    "the spiral field moves 0.003 %% of peak (the bounding plates, which do the "
                    "shielding, stay), while the mesh drops ~26 %% -- at geo-h 4 mm that is 43.3k "
                    "elements and a 15.0 GB dense block, right at this card's 16 GiB, against ~32k "
                    "and 8.3 GB with them out")
p.add_argument("--maxfun", type=int, default=80)
p.add_argument("--rhobeg", type=float, default=0.12, help="initial trust region as a fraction of the bound ranges")
p.add_argument("--rhoend", type=float, default=0.01)
p.add_argument("--t-floor", type=float, default=None, help="transmission floor (fraction) of the penalty; default seed - 0.05")
p.add_argument("--bounds", type=float, nargs="+", default=None,
               metavar="B", help="lo hi per knob in the vector's order (see the log's 'bounds' line), e.g. for the legacy vector: "
                    "tilt[deg] angling[deg] gamma[deg] sigma[mm] t_ent[deg] t_exit[deg] dz[mm] q1..qN[V] gap[mm] vscale; "
                    "default: the doublet values, or +-35 %% around the seed voltages for three or more quads")
p.add_argument("--bound", nargs=3, action="append", default=[], metavar=("NAME", "LO", "HI"),
               help="override one knob's bounds, e.g. --bound ang_exit -28 0; repeatable")
p.add_argument("--profiles", nargs="+", default=[], choices=["angling", "sigma"],
               help="free these as 3-node profiles along the electrode (entrance / middle / exit) instead of scalars: angling -> "
                    "ang_ent ang_mid ang_exit [deg], sigma -> sig_ent_mm sig_mid_mm sig_exit_mm [mm] (Barnard et al. 2021's Q2(s), Q1(s))")
p.add_argument("--profile-nodes", type=int, default=3,
               help="nodes per profile: 3 (entrance / middle / exit) or 5 (+ the quarter points, knobs *_q1 / *_q3; Daniel 2026-10-23). "
                    "A 3-node seed is carried onto 5 nodes exactly (its nodes are a subset)")
p.add_argument("--free-aspect", action="store_true",
               help="electrode width / gap (geometry knob aspect) as a knob; with --min-gap the lever on the anode-cathode end gap")
p.add_argument("--set", action="append", default=[], metavar="NAME=VALUE",
               help="override a seed value (knob units as logged), e.g. --set ang_mid=6; repeatable. For smoke tests: the "
                    "levelled reference is the right seed for an optimization")
p.add_argument("--particles", default=None, help="RFQ file (default: the levelling's, else the deck's ext3_exit.dst)")
p.add_argument("--bfield", default=None, help="magnet map (default: the levelling's, else the deck's final1_baseline_1mm.pickle)")
p.add_argument("--core-only", action="store_true",
               help="the OLD loop: --n particles of the RFQ CORE only and core-only terms. Default since 2026-10-21: --n drawn at "
                    "random from the whole file (core + tail, bunch.py --stragglers) and every term over all of them")
p.add_argument("--sc-min-particles", type=int, default=200,
               help="whole-file loop: stop re-solving the space charge below this many depositing particles (the late tail)")
p.add_argument("--seed-only", action="store_true", help="evaluate the seed, write the report, no DFO-LS (smoke tests)")
p.add_argument("--final", action="store_true",
               help="after DFO-LS, re-run the WINNER with EVERY particle of the RFQ file (core + tail, --stragglers) at --current-ma: "
                    "its spiral basis re-solved, SCLevel --final's space-charge settings (--final-h, --final-resolve-every), hand-off "
                    "file; -> <out>/final/ (bunch_sc_all_final.*, final.json) and a section in the report. Daniel 2026-10-21: a core-only "
                    "loop, the winner confirmed with the full distribution")
p.add_argument("--final-only", action="store_true", help="no DFO-LS: take the best evaluation of the cache and do the --final run")
p.add_argument("--interim", action="store_true",
               help="while an optimization is still RUNNING on this cache: re-run its current best evaluation with every particle into "
                    "<out>/interim_<tag>/ (report, best_knobs.json and the final run there); nothing of the running optimizer is "
                    "written (no cache save, no report, no final/). Daniel 2026-10-23: an intermediate result")
p.add_argument("--final-n", type=int, default=10 ** 7, help="particles of the final run (default: every particle of the file; less only to test)")
p.add_argument("--final-h", type=float, default=0.002, help="space-charge cell of the final every-particle run [m] (SCLevel --final)")
p.add_argument("--final-resolve-every", type=int, default=8, help="Poisson solve cadence of the final every-particle run (SCLevel --final)")
p.add_argument("--dry-run", action="store_true",
               help="print the knob vector, seed, bounds, residual layout and the seed evaluation's commands; build nothing")
p.add_argument("--keep-basis", action="store_true", help="keep every evaluation's basis field (large)")
p.add_argument("--keep-npz", action="store_true",
               help="keep each evaluation's hand-off npz (~1 MB at n=4000) instead of deleting it. The npz carries "
                    "handoff_time / handoff_state, so any hand-off quantity can be recomputed later without re-running. "
                    "2026-09-16: their absence is why bunch length was unrecoverable for the 135 evaluations of 09-13/15, "
                    "which measured only what the objective of the day asked for")
p.add_argument("--out-name", default="full_opt",
               help="name of the output / cache folder inside the run (default full_opt). A run with a different objective "
                    "(e.g. --w-de) must use its own folder: cached measures of an older objective lack the new terms")
p.add_argument("--w-de", type=float, default=0.0,
               help="weight of the hand-off energy-spread residual (rms kinetic energy at the hand-off plane / 1 keV). 0 (default) = off. "
                    "2026-09-15: the central region's bar ends take the particles below ~65 keV (29 %% of the doublet hand-off, 3.3 keV "
                    "rms, 71 %% correlated with height); a 3 keV beam at weight 3 weighs like the 33 mrad angle residual")
p.add_argument("--w-angle", type=float, default=1.0,
               help="weight of the hand-off vertical-angle residual (rms z' at the hand-off plane / 3 mrad); 1 = the campaign default. "
                    "2026-09-15: lowered for the 'minimum phase width' variant so the bunch-length term can lead")
p.add_argument("--robust-angle", action="store_true",
               help="the angle residual takes the ROBUST width (p84 - p16) / 2 of z' at the hand-off instead of the rms. 2026-10-21: "
                    "in a 4000-particle whole-file draw a few strays carry the rms (t1nr_circ_c2 seed: 38.1 mrad rms, 25.9 robust; "
                    "every particle 33.1 / 27.5), so the rms term can jump when one stray comes or goes. Off = the reported metric")
p.add_argument("--w-dt", type=float, default=0.0,
               help="weight of the hand-off bunch-length residual: rms RF phase at the hand-off plane, wrapped to +-180 deg around "
                    "the core's median, in RF degrees / 10 deg. 0 (default) = off. 2026-09-15: the central region's transmission "
                    "follows the bunch length at the hand-off (52 deg rms -> 39.8 %%, 26 deg -> 66 %%, about one point per degree). "
                    "DFO-LS squares it: at 33 mrad / 48 deg, --w-dt 2.3 weighs the phase like the angle. Needs its own --out-name")
p.add_argument("--w-loss", type=float, default=0.0,
               help="weight of a residual on ALL particles lost, w * 100 * (1 - T) (one point = one unit), QuadOptN's term. 2026-10-24: "
                    "the sp1 shape round (no such term) was blind to the 14-point space-charge loss on the grounded SI entrance aperture. "
                    "0 = off; needs its own --out-name")
p.add_argument("--w-hv", type=float, default=0.0,
               help="weight of the fraction lost on HV electrodes (spiral electrodes AND quad poles) / 1 point, QuadOptN's term (Daniel "
                    "wants losses on grounded apertures). 0 = off; needs its own --out-name")
p.add_argument("--min-gap", type=float, default=0.0,
               help="minimum anode-cathode distance [mm] (Daniel, 2026-10-15: the electrodes nearly touched at the exit end - a sure "
                    "spark). Measured exactly on the evaluation's solids; 0 = off")
p.add_argument("--w-gap", type=float, default=10.0, help="weight of the one-sided gap residual: w * (min_gap - gap) per mm short")
p.add_argument("--w-fit", type=float, default=0.0, help="weight of the pole-steel residual: conductor mesh vertices where |B| > 1.6 T in the magnet map, after R and dz (HCHC-1.5: the housing top rim sat 8 mm inside a hill at the HCHC-60 seed). 0 = off; 1 makes one vertex one unit")
p.add_argument("--w-z100", type=float, default=0.0,
               help="weight of an aperture proxy residual: rms ballistic height 100 mm past the measurement plane (z + 0.1 z') / 5 mm. "
                    "0 (default) = off. The exit metrics alone do not rank designs for the cyclotron (2026-09-14: a variant gained 2.5 "
                    "core points at the exit and lost 11 points of vertical survival in the machine); this proxy penalizes a beam that "
                    "has spread vertically by the first-turn apertures")
p.add_argument("--free-height", action="store_true",
               help="the design orbit's electric radius A [mm] (the inflector height) as a knob height_mm; the design voltage becomes "
                    "T gap / A. Off: the run's V/gap, A fixed. Needs its own --out-name")
p.add_argument("--r-target-mm", type=float, default=None, help="hand-off radius target [mm] (mean position of the transmitted beam)")
p.add_argument("--w-r", type=float, default=0.0, help="weight of (hand-off r - target) / 1 mm; 0 = off")
p.add_argument("--alpha-target-deg", type=float, default=None,
               help="target angle of the mean momentum at the hand-off from the local tangent, + = outward [deg]")
p.add_argument("--w-alpha", type=float, default=0.0, help="weight of (alpha - target) / 1 deg; 0 = off")
a = p.parse_args()
if (a.w_r > 0 and a.r_target_mm is None) or (a.w_alpha > 0 and a.alpha_target_deg is None):
    raise SystemExit("--w-r needs --r-target-mm and --w-alpha needs --alpha-target-deg")

# NAMES / UNITS / LO / HI are set after the seed is read: the run says how many quads there are

RUN = os.path.abspath(a.run)
BASE = os.path.abspath(a.base_run) if a.base_run else RUN
BASIS0, STEPS0 = os.path.join(RUN, "basis"), os.path.join(RUN, "steps")
OUT = os.path.join(RUN, a.out_name)
WHOLE = not a.core_only
if not a.dry_run:
    os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8") if not a.dry_run else None


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    if LOG is not None:
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
VSCALE0 = 1.0                                          # operating voltage / centering value; the levelling may record one
lvl_d = {}
if os.path.exists(lvl):
    lvl_d = jload(lvl)
    t_exit_seed, dz_seed = float(lvl_d["t_exit_deg"]), float(lvl_d["dz_m"])
    VSCALE0 = float(lvl_d.get("vscale", 1.0))
PARTICLES = a.particles or lvl_d.get("particles") or os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst")
BFIELD = a.bfield or lvl_d.get("bfield") or os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle")
R = jload(os.path.join(RUN, "exit_plane.json"))["rotation_deg"]
PHI = 90.0 + R
ROT_BF = os.path.join(BASE, "bfield_rotated.pickle")
ENERGY = float(particle_rows(PARTICLES, core=True)[:, 8].mean())   # the DESIGN energy (core mean), not a beam number
NQ = int(knobs0.get("n_quads", 2) or 2)
QN = ["q{}".format(i + 1) for i in range(NQ)]
# quad seed: an explicit scan json, else the levelled setting (voltages AND skew turns), else the old default scan json
if a.best_json:
    bj = a.best_json if os.path.isabs(a.best_json) else os.path.join(RUN, a.best_json)
    best0 = jload(bj)["best"]
    qseed = [float(best0[q]) for q in QN]
    QALPHA = [float(best0.get("alpha{}".format(i + 1), 0.0)) for i in range(NQ)]
    QSRC = bj
elif lvl_d.get("quads") and len(lvl_d["quads"]) == 2 * NQ:
    qseed = [float(v) for v in lvl_d["quads"][0::2]]
    QALPHA = [float(v) for v in lvl_d["quads"][1::2]]
    QSRC = lvl
else:
    bj = os.path.join(RUN, "basis", "scan2_coarse.json")
    best0 = jload(bj)["best"]
    qseed = [float(best0[q]) for q in QN]
    QALPHA = [0.0] * NQ
    QSRC = bj
QROT = list(a.quads) if a.quads else [float(v) for v in (lvl_d.get("quad_rotation_deg") or [16.0, 24.0])]
if len(QROT) != NQ:
    raise SystemExit("--quads needs {} angles for this {}-quad run, got {}".format(NQ, NQ, len(QROT)))

# ---------------------------------------------------------------- the knob vector
for _key in ("angling", "sigma"):
    if knobs0.get(_key + "_profile") is not None and _key not in a.profiles:
        raise SystemExit("this run's geometry carries a {0}_profile, which would override the scalar {0} knob: "
                         "run with --profiles {0}".format(_key))
# profile nodes spread evenly from the entrance (xi -1) to the exit (+1); 5 adds the quarter points (2026-10-23, Daniel)
NODE_TAGS = {3: ["ent", "mid", "exit"], 5: ["ent", "q1", "mid", "q3", "exit"]}
if a.profile_nodes not in NODE_TAGS:
    raise SystemExit("--profile-nodes: 3 or 5, got {}".format(a.profile_nodes))
TAGS = NODE_TAGS[a.profile_nodes]
ANG = ["ang_" + t for t in TAGS] if "angling" in a.profiles else ["angling"]
SIG = ["sig_{}_mm".format(t) for t in TAGS] if "sigma" in a.profiles else ["sigma_mm"]
NAMES = (["tilt"] + ANG + ["gamma"] + SIG + ["t_ent", "t_exit", "dz_mm"] + QN + ["gap_mm", "vscale"] + (["aspect"] if a.free_aspect else [])
         + (["height_mm"] if a.free_height else []))
UNIT_OF = dict(tilt="deg", angling="deg", gamma="deg", sigma_mm="mm", t_ent="deg", t_exit="deg", dz_mm="mm", gap_mm="mm", vscale="x",
               aspect="x", height_mm="mm", **{"ang_" + t: "deg" for t in NODE_TAGS[5]}, **{"sig_{}_mm".format(t): "mm" for t in NODE_TAGS[5]},
               **{q: "V" for q in QN})
UNITS = [UNIT_OF[n] for n in NAMES]

seed_val = {"tilt": knobs0["tilt"], "angling": knobs0["angling"], "gamma": knobs0["gamma"], "sigma_mm": 1e3 * knobs0["sigma"],
            "t_ent": t_ent0, "t_exit": t_exit_seed, "dz_mm": 1e3 * dz_seed, "gap_mm": 1e3 * knobs0["gap"], "vscale": VSCALE0,
            "aspect": knobs0.get("aspect", 2.4), **dict(zip(QN, qseed)),
            # A = T gap / V for q = 1 (SpiralInflector.initialize's height), T the design energy of build_geometry
            "height_mm": 1e3 * 1e6 * ENERGY * float(knobs0["gap"]) / float(knobs0["volt"])}
# profile seeds: the run's own profile if it has one, else the scalar's legacy meaning (angling a: the ramp +a -> -a,
# i.e. [a, 0, -a]; a <= 0 meant "no angling"; sigma s: constant)
_ap = knobs0.get("angling_profile")
if _ap is None:
    _a = float(knobs0["angling"])
    _ap = [_a, 0.0, -_a] if _a > 0.0 else [0.0, 0.0, 0.0]
_sp = knobs0.get("sigma_profile")
_sp = [1e3 * float(v) for v in _sp] if _sp is not None else [1e3 * float(knobs0["sigma"])] * 3
_xi = np.linspace(-1.0, 1.0, a.profile_nodes)
for _name, _prof in (("angling", _ap), ("sigma", _sp)):
    _old = np.linspace(-1.0, 1.0, len(_prof))
    if not np.all(np.isclose(_old[:, None], _xi[None, :]).any(axis=1)):
        raise SystemExit("the run's {}-node {} profile is not a subset of {} nodes: the seed would change the geometry".format(
            len(_prof), _name, a.profile_nodes))
# the run's profile carried onto the node set by linear interpolation: exact while its nodes are a subset (3 -> 5)
seed_val.update({"ang_" + t: float(v) for t, v in zip(TAGS, np.interp(_xi, np.linspace(-1.0, 1.0, len(_ap)), _ap))})
seed_val.update({"sig_{}_mm".format(t): float(v) for t, v in zip(TAGS, np.interp(_xi, np.linspace(-1.0, 1.0, len(_sp)), _sp))})
for item in a.set:
    _k, _v = item.split("=", 1)
    if _k not in NAMES:
        raise SystemExit("--set {}: not a knob of this vector ({})".format(_k, ", ".join(NAMES)))
    seed_val[_k] = float(_v)

# default bounds (the legacy ones for the legacy knobs)
BOUND0 = {"tilt": (20, 42), "angling": (8, 22), "gamma": (3, 13), "sigma_mm": (0.2, 1.6), "t_ent": (0.1, 0.8), "t_exit": (0.3, 1.6),
          "dz_mm": (2.0, 7.0), "gap_mm": (16, 24), "vscale": (0.96, 1.06),
          # 2026-10-21 profiles: signed, wider than the scalars (HCHC-1.5 brief: open the parameter space); the ends are bound
          # in practice by --min-gap (the end gap closes ~0.7 mm per degree of angling at aspect 2.4, less when narrow)
          "ang_ent": (-5, 30), "ang_q1": (-15, 30), "ang_mid": (-15, 15), "ang_q3": (-30, 15), "ang_exit": (-30, 5),
          "sig_ent_mm": (-1.0, 2.5), "sig_q1_mm": (-1.0, 2.5), "sig_mid_mm": (-1.0, 2.5), "sig_q3_mm": (-1.0, 2.5),
          "sig_exit_mm": (-1.0, 2.5),
          "aspect": (1.2, 2.4), "height_mm": (85.0, 110.0)}
if NQ == 2:
    BOUND0.update(q1=(5500, 7500), q2=(-7700, -5700))
else:                                     # +-35 % around each seed voltage, sign preserved
    BOUND0.update({q: tuple(sorted([0.65 * v, 1.35 * v])) for q, v in zip(QN, qseed)})
if a.bounds is not None:
    if len(a.bounds) != 2 * len(NAMES):
        raise SystemExit("--bounds needs {} numbers (lo hi per knob: {}), got {}".format(2 * len(NAMES), " ".join(NAMES), len(a.bounds)))
    BOUND0.update({n: (a.bounds[2 * i], a.bounds[2 * i + 1]) for i, n in enumerate(NAMES)})
for _n, _lo, _hi in a.bound:
    if _n not in NAMES:
        raise SystemExit("--bound {}: not a knob of this vector ({})".format(_n, ", ".join(NAMES)))
    BOUND0[_n] = (float(_lo), float(_hi))
LO = np.array([BOUND0[n][0] for n in NAMES], dtype=float)
HI = np.array([BOUND0[n][1] for n in NAMES], dtype=float)
if np.any(HI <= LO):
    raise SystemExit("empty bounds for " + ", ".join(n for n, lo, hi in zip(NAMES, LO, HI) if hi <= lo))
_sb = Field.from_file(os.path.join(BASIS0, "ef_itp_spiral.pickle")) if not a.dry_run else None
if _sb is not None:
    _g = _sb.grid
    GRID_RES = float(np.asarray(_g["x"], dtype=float)[1] - np.asarray(_g["x"], dtype=float)[0])
    GRID_BOX = [float(np.asarray(_g[c], dtype=float)[i]) for c in "xyz" for i in (0, -1)]
    del _sb
else:                                     # the dry run reads the grid from the reload record instead of the 1 GB field
    _rj = jload(os.path.join(BASIS0, "reload_spiral.json"))
    GRID_RES, GRID_BOX = float(_rj.get("res_m", 0.0025)), ["<box of basis/ef_itp_spiral.pickle>"]
VOLT_PER_M = float(knobs0["volt"]) / float(knobs0["gap"])   # the design voltage scales with the gap (same field)

x0_raw = np.array([seed_val[n] for n in NAMES], dtype=float)
x0 = np.clip(x0_raw, LO, HI)
_clipped = [n for n, u, c in zip(NAMES, x0_raw, x0) if u != c]
# The optimizer sees only the free knobs; evaluate(), the cache and the report keep working on the
# full vector, so an --out-name folder is shared between a held and an unheld run.
HOLD = set(a.hold) | (set(QN) if a.hold_quads else set())
if HOLD - set(NAMES):
    raise SystemExit("--hold: not knobs of this vector: " + ", ".join(sorted(HOLD - set(NAMES))))
FREE = [i for i, n in enumerate(NAMES) if n not in HOLD]
HELD = [n for n in NAMES if n in HOLD]


def expand(xf):
    x = np.array(x0, dtype=float)
    x[FREE] = xf
    return x


def geometry_knobs(vals):
    """The build_geometry knobs of one full knob vector (as a name -> value dict)."""
    k = dict(knobs0)
    # the design voltage: the run's field V/gap (A fixed), or with --free-height the field of the electric radius height_mm
    volt = (1e6 * ENERGY * vals["gap_mm"] / vals["height_mm"]) if "height_mm" in vals else VOLT_PER_M * 1e-3 * vals["gap_mm"]
    k.update(tilt=vals["tilt"], gamma=vals["gamma"], gap=1e-3 * vals["gap_mm"], volt=volt)
    if "angling" in vals:
        k["angling"] = vals["angling"]
    else:
        prof = [vals[n] for n in ANG]
        # the scalar is ignored by the build once a profile is given; it keeps the legacy-equivalent amplitude for tables
        k.update(angling_profile=prof, angling=0.5 * (prof[0] - prof[-1]))
    if "sigma_mm" in vals:
        k["sigma"] = 1e-3 * vals["sigma_mm"]
    else:
        prof = [1e-3 * vals[n] for n in SIG]
        k.update(sigma_profile=prof, sigma=float(np.mean(prof)))
    if "aspect" in vals:
        k["aspect"] = vals["aspect"]
    return k


SKEW = {i: os.path.exists(os.path.join(BASIS0, "ef_itp_q{}skew.pickle".format(i + 1))) for i in range(NQ)}
for i, al in enumerate(QALPHA):
    if al != 0.0 and not SKEW[i] and not os.path.exists(os.path.join(BASIS0, "basis_meta.json")):
        raise SystemExit("quad {} is turned by alpha {:+.4f} deg but {} has no ef_itp_q{}skew.pickle".format(i + 1, al, BASIS0, i + 1))

stamp("FULL OPT of {}: R {:+.3f} deg fixed, quads rotated {} (+ alpha {} held), {} mA, {} particles {}, {:.0f} mm cells, basis grid {:.1f} mm, "
      "mesh {:.1f} mm, design-orbit grid {:.1f} mm, quad poles {} of the spiral solve; DFO-LS maxfun {}, rhobeg {} of range".format(
          os.path.basename(RUN), R, QROT, QALPHA, a.current_ma, a.n,
          "drawn from the WHOLE RFQ file (core + tail, every term over all of them)" if WHOLE else "of the CORE only (old loop)",
          1e3 * a.h, 1e3 * GRID_RES, 1e3 * a.geo_h, 1e3 * a.geo_res, "OUT" if a.drop_quad_poles else "in", a.maxfun, a.rhobeg))
stamp("  quad seed from {}; particles {}; B map {}".format(QSRC, PARTICLES, BFIELD))
if a.profiles or a.free_aspect:
    stamp("  shape knobs: {}{} (profile nodes: entrance / middle / exit section of the truncated electrode, linear between)".format(
        ", ".join(a.profiles) + " as {}-node profiles".format(a.profile_nodes) if a.profiles else "scalars", " + aspect" if a.free_aspect else ""))
stamp("  seed: " + ", ".join("{} {:g} {}".format(n, v, u) for n, v, u in zip(NAMES, x0, UNITS)))
if _clipped:
    stamp("  WARNING: seed clipped into the bounds for " + ", ".join(_clipped))
stamp("  bounds: " + ", ".join("{} [{:g}, {:g}]".format(n, lo, hi) for n, lo, hi in zip(NAMES, LO, HI)))
_res_names = (["{}angle/3mrad x{:g}".format("robust " if a.robust_angle else "", a.w_angle), "floor", "height/0.3mm", "angle centroid/2mrad"] + (["z100/5mm"] if a.w_z100 > 0 else [])
              + (["dE/1keV x{:g}".format(a.w_de)] if a.w_de > 0 else []) + (["phase/10deg x{:g}".format(a.w_dt)] if a.w_dt > 0 else [])
              + (["lost/1pt x{:g}".format(a.w_loss)] if a.w_loss > 0 else []) + (["lost on HV/1pt x{:g}".format(a.w_hv)] if a.w_hv > 0 else [])
              + (["steel"] if a.w_fit > 0 else []) + (["gap<{:g}mm x{:g}/mm".format(a.min_gap, a.w_gap)] if a.min_gap > 0 else [])
              + (["hand-off r-{:g}mm/1mm x{:g}".format(a.r_target_mm, a.w_r)] if a.w_r > 0 else [])
              + (["alpha-{:g}deg/1deg x{:g}".format(a.alpha_target_deg, a.w_alpha)] if a.w_alpha > 0 else []))
stamp("  residuals: {}; {} free knobs -> DFO-LS's initial model takes {} evaluations".format(" | ".join(_res_names), len(FREE), len(FREE) + 1))

if a.dry_run:
    k_seed = geometry_knobs({n: float(v) for n, v in zip(NAMES, x0)})
    stamp("DRY RUN, nothing built. Seed geometry knobs that differ from the run's: " + ", ".join(
        "{} {} -> {}".format(kk, knobs0.get(kk), vv) for kk, vv in k_seed.items() if knobs0.get(kk) != vv))
    _super = [v for q, al in zip(qseed, QALPHA) for v in (q, al)]
    stamp("  per evaluation: build_geometry(fix_truncations=(t_ent, t_exit), fix_dz={:.6f} m, maxiter=0, res={}, h={}); then".format(dz0, a.geo_res, a.geo_h))
    stamp("  python -m spyral_inflector.tracking.bem_reload --tag spiral --step-dir <e>/steps --res {} --box <run's> --h {} --bfield {} "
          "--rotate-all {} --rotate-quads {} --quad-voltages {} --energy-mev {:.6f}{}".format(
              GRID_RES, a.geo_h, BFIELD, R, " ".join(map(str, QROT)), " ".join(["0"] * NQ), ENERGY, " --drop D*" if a.drop_quad_poles else ""))
    stamp("  python -m spyral_inflector.tracking.bunch --tag spiral --particles {} --bfield {} --phi {} --superpose {} --current-ma {} "
          "--handoff-frame machine --handoff-distance 0.0 --n {}{}{} --h {} --resolve-every {} --shift-z <dz - dz0> --vscale <vscale> --no-plot".format(
              PARTICLES, BFIELD, PHI, " ".join("{:g}".format(v) for v in _super), a.current_ma, a.n,
              " --stragglers --sc-min-particles {}".format(a.sc_min_particles) if WHOLE else "", " --sc" if a.current_ma > 0 else "", a.h,
              a.resolve_every))
    stamp("  copied per evaluation: " + ", ".join(["ef_itp_q{}.pickle".format(i + 1) for i in range(NQ)] +
                                                 ["ef_itp_q{}skew.pickle".format(i + 1) for i in range(NQ) if QALPHA[i] != 0.0 and SKEW[i]]))
    stamp("  output folder {}".format(OUT))
    sys.exit(0)

# ---------------------------------------------------------------- the evaluation cache
CACHE = os.path.join(OUT, "evals.json")
META = os.path.join(OUT, "meta.json")
evals = jload(CACHE) if os.path.exists(CACHE) else []
meta_now = {"names": NAMES, "whole_file": WHOLE, "n": a.n, "current_ma": a.current_ma, "alphas": QALPHA, "particles": PARTICLES}
if evals:
    # an older folder (no meta.json) is the legacy loop: core only, the scalar vector
    meta_old = jload(META) if os.path.exists(META) else {"names": None, "whole_file": False}
    if meta_old.get("names") is None and len(evals[0]["x"]) != len(NAMES):
        meta_old["names"] = "a {}-knob vector".format(len(evals[0]["x"]))
    clash = [k for k in ("whole_file", "names", "n", "current_ma", "alphas")
             if k in meta_old and meta_old[k] is not None and meta_old[k] != meta_now[k]]
    if clash:
        raise SystemExit("{} holds evaluations of another loop ({}: {} vs now {}); use another --out-name (a folder from before "
                         "2026-10-21, without meta.json, is the core-only loop: resume it with --core-only)".format(
                             CACHE, ", ".join(clash), [meta_old[k] for k in clash], [meta_now[k] for k in clash]))
with open(META, "w") as fh:
    json.dump(meta_now, fh, indent=1)


def key_of(x):
    return hashlib.md5(json.dumps([round(float(v), 4) for v in x]).encode()).hexdigest()[:10]


def save_cache():
    with open(CACHE, "w") as fh:
        json.dump(evals, fh, indent=1)


SI_ELECTRODES = ("SI_Anode", "SI_Cathode")


def measure(out_dir, tag, core_only=None):
    core_only = a.core_only if core_only is None else core_only
    s = jload(os.path.join(out_dir, "bunch_{}.json".format(tag)))
    d = np.load(os.path.join(out_dir, "bunch_{}.npz".format(tag)), allow_pickle=True)
    st = d["asym_state"]
    crossed = np.asarray(d["crossed"], dtype=bool)         # through the housing exit opening (post-exit hits excluded)
    tail = np.asarray(d["is_tail"], dtype=bool) if "is_tail" in d.files else np.zeros(len(crossed), bool)
    sel = crossed & ~tail if core_only else crossed        # whole-file loop: every transmitted particle, tail included
    ok = sel & np.all(np.isfinite(st), axis=1)
    z, v = st[ok, 2], st[ok, 3:6]
    zp = v[:, 2] / np.hypot(v[:, 0], v[:, 1])
    # the hand-off plane (electrode exit, machine frame) is where the cyclotron sets its injection cut: the laptop's
    # 2026-09-14 A/Bs put the survival cliff at ~40 mrad rms of z' THERE, while the 55 mm plane compresses beams that
    # differ by 11 points of survival into 33-36 mrad. Angle and proxy are therefore measured at the hand-off plane.
    hs = d["handoff_state"] if "handoff_state" in d.files else None
    okh = sel & np.all(np.isfinite(hs), axis=1) if hs is not None else np.zeros(len(sel), bool)   # beam_report.py's selection
    # 2026-09-15: the energy spread at the hand-off is the longitudinal channel of the central region (the bar ends take
    # the particles below ~65 keV: 29 % of the doublet's hand-off, 3.3 keV rms, 71 % correlated with the height because
    # the spiral's exit potential varies with z). Reported always, in the objective only with --w-dE > 0.
    M_H2P_KG = 3.3462e-27                                             # 2 protons + 1 electron
    dT_rms_keV, f_low, dt_rms_deg, dt_robust_deg, zp_h_robust = 0.0, 0.0, 0.0, 0.0, 0.0
    r_h_mm, alpha_deg = None, None
    if okh.sum() >= 100:
        # 2026-11-10: where the beam is handed to the central region and which way it points there (x-y plane, machine frame):
        # the radius of the mean position and the angle of the MEAN momentum from the local tangent, + = outward
        pm, vm = hs[okh, 0:2].mean(axis=0), hs[okh, 3:5].mean(axis=0)
        r_h_mm = float(1e3 * np.linalg.norm(pm))
        rhat, vhat = pm / np.linalg.norm(pm), vm / np.linalg.norm(vm)
        that = np.sign(rhat[0] * vhat[1] - rhat[1] * vhat[0]) * np.array([-rhat[1], rhat[0]])   # along the circulation
        alpha_deg = float(np.degrees(np.arctan2(vhat @ rhat, vhat @ that)))
        zh, vh = hs[okh, 2], hs[okh, 3:6]
        zph = vh[:, 2] / np.hypot(vh[:, 0], vh[:, 1])
        zp_h_rms, z100 = float(zph.std()), zh + 0.1 * zph          # ballistic height 100 mm past the hand-off plane
        zp_h_robust = float(0.5 * np.subtract(*np.percentile(zph, [84.135, 15.865])))
        T_keV = 0.5 * M_H2P_KG * np.sum(vh ** 2, axis=1) / 1.602176634e-19 / 1e3
        dT_rms_keV, f_low = float(T_keV.std()), float(np.mean(T_keV < 0.965 * T_keV.mean()))
        # 2026-09-15: the bunch length at the hand-off in RF degrees (the central region's longitudinal channel);
        # the arrival times are on the tracker's clock (the RFQ file's own phase spread included). 2026-10-21: WRAPPED to
        # +-180 deg around the core's median, as beam_report.py reports it -- once the tail is in the loop, a straggler
        # N periods late sees the same RF phase as one on time
        if "handoff_time" in d.files:
            th = np.asarray(d["handoff_time"], dtype=float)[okh]
            fin = np.isfinite(th)
            th, core_h = th[fin], ~tail[okh][fin]
            f_rf = float(s.get("rf_mhz") or s.get("rf_frequency_mhz") or 32.8) * 1e6
            if len(th) >= 100:
                ph = 360.0 * f_rf * th
                ref = np.median(ph[core_h]) if core_h.any() else np.median(ph)
                phw = (ph - ref + 180.0) % 360.0 - 180.0
                dt_rms_deg = float(phw.std())
                dt_robust_deg = float(0.5 * np.subtract(*np.percentile(phw, [84.135, 15.865])))
    else:
        zp_h_rms, z100 = float(zp.std()), z + 0.1 * zp
    n_all = float(s.get("n_particles") or len(crossed))
    lbe = s.get("losses_by_electrode", {})
    return {"z": float(z.mean()), "zp": float(zp.mean()), "z_rms": float(z.std()), "zp_rms": float(zp.std()), "zp_h_rms": zp_h_rms,
            "zp_h_robust": zp_h_robust, "z100_rms": float(z100.std()), "dT_rms_keV": dT_rms_keV, "f_low_energy": f_low,
            "dt_rms_deg": dt_rms_deg, "dt_robust_deg": dt_robust_deg, "r_h_mm": r_h_mm, "alpha_deg": alpha_deg,
            "transmission": float(crossed[~tail].mean()) if core_only and tail.any() else float(s["transmission"]),
            "n": int(ok.sum()), "n_particles": int(n_all), "n_tail": int(tail.sum()), "n_tail_through": int((crossed & tail).sum()),
            "whole_file": not core_only, "f_lost_si": sum(lbe.get(k_, 0) for k_ in SI_ELECTRODES) / n_all,
            "f_lost_hv": sum(v_ for k_, v_ in lbe.items() if k_ in SI_ELECTRODES or (k_.startswith("D") and k_[1:].isdigit())) / n_all,
            "voltage_V": s.get("voltages", {}).get("SI_Anode")}


T_FLOOR = a.t_floor          # set after the seed evaluation if None


def residuals(m):
    if m is None:
        return np.array([10.0] * (4 + (1 if a.w_z100 > 0 else 0) + (1 if a.w_de > 0 else 0) + (1 if a.w_dt > 0 else 0)
                                  + (1 if a.w_loss > 0 else 0) + (1 if a.w_hv > 0 else 0)
                                  + (1 if a.w_fit > 0 else 0) + (1 if a.min_gap > 0 else 0)
                                  + (1 if a.w_r > 0 else 0) + (1 if a.w_alpha > 0 else 0)))
    ang = m.get("zp_h_rms", m["zp_rms"])
    if a.robust_angle:
        ang = m.get("zp_h_robust") or ang
    r = [a.w_angle * 1e3 * ang / 3.0, max(0.0, (T_FLOOR or 0.0) - m["transmission"]) / 0.01, 1e3 * m["z"] / 0.3, 1e3 * m["zp"] / 2.0]
    if a.w_z100 > 0:
        r.append(a.w_z100 * 1e3 * m.get("z100_rms", 0.0) / 5.0)
    if a.w_de > 0:
        r.append(a.w_de * m.get("dT_rms_keV", 0.0) / 1.0)            # 1 keV rms of hand-off energy spread = one unit
    if a.w_dt > 0:
        r.append(a.w_dt * m.get("dt_rms_deg", 0.0) / 10.0)           # 10 RF degrees rms of hand-off arrival phase = one unit
    if a.w_loss > 0:
        r.append(a.w_loss * 100.0 * (1.0 - m["transmission"]))       # one point of transmission lost = one unit (QuadOptN's term)
    if a.w_hv > 0:
        r.append(a.w_hv * 100.0 * m.get("f_lost_hv", m.get("f_lost_si", 0.0)))   # one point lost on HV electrodes = one unit
    if a.w_fit > 0:
        r.append(a.w_fit * float(m.get("n_steel", 0)) / 1.0)          # one conductor mesh vertex in pole steel = one unit
    if a.min_gap > 0:
        g = m.get("electrode_gap_mm")
        r.append(a.w_gap * max(0.0, a.min_gap - (g if g is not None else 0.0)))   # one-sided: nothing above the minimum
    if a.w_r > 0:                                                    # 1 mm off the hand-off radius target = one unit
        rh_ = m.get("r_h_mm")
        r.append(a.w_r * ((rh_ - a.r_target_mm) if rh_ is not None else 10.0))
    if a.w_alpha > 0:                                                # 1 deg off the exit-direction target = one unit
        al_ = m.get("alpha_deg")
        r.append(a.w_alpha * ((al_ - a.alpha_target_deg) if al_ is not None else 10.0))
    return np.array(r)


# 2026-10-14 (HCHC-1.5): the BEM and the tracker know nothing about the magnet's iron, so an
# optimizer left to itself will happily put the housing into a pole. This counts the conductor
# mesh vertices (the state pickle carries the mesh) that sit where |B| > B_IRON in the UNROTATED
# map, after the system rotation R and the axial shift of this evaluation. The e071 seed on
# HCHC-1.5 had 4 such vertices in the housing's top rim; its FullOpt must see them.
B_IRON = 1.6
_BF_FIT = None


def steel_vertices(state_pickle, R_deg, dz_m):
    global _BF_FIT
    if _BF_FIT is None:
        from spyral_inflector.tracking.deck import load_bfield
        _BF_FIT = load_bfield(BFIELD)
    with open(state_pickle, "rb") as fh:
        mesh = pickle.load(fh)["mesh"]
    V = np.asarray(mesh["verts"], dtype=float)
    E = np.asarray(mesh["elems"], dtype=int)
    if V.shape[0] == 3 and V.shape[1] != 3:
        V = V.T
    if E.shape[0] == 3 and E.shape[1] != 3:
        E = E.T
    V = V[np.unique(E)]
    c, s = np.cos(np.radians(R_deg)), np.sin(np.radians(R_deg))
    v = np.empty_like(V)
    v[:, 0] = c * V[:, 0] - s * V[:, 1]
    v[:, 1] = s * V[:, 0] + c * V[:, 1]
    v[:, 2] = V[:, 2] + dz_m
    return int((np.linalg.norm(_BF_FIT(v), axis=1) > B_IRON).sum())


def spiral_basis(gdir, sdir, bdir, log_path):
    """Re-solve the spiral basis of one built geometry on the run's basis grid and copy the run's quad bases next to it."""
    run_cmd(["-m", "spyral_inflector.tracking.bem_reload", "--tag", "spiral", "--step-dir", sdir, "--voltages", os.path.join(gdir, "voltages.csv"),
             "--state", os.path.join(gdir, "state.pickle"), "--out-dir", bdir, "--res", GRID_RES, "--box"] + GRID_BOX +
            ["--h", a.geo_h, "--bfield", BFIELD, "--rotate-all", R, "--rotate-quads"] + QROT + ["--quad-voltages"] + [0] * NQ +
            ["--energy-mev", ENERGY] + (["--drop", "D*"] if a.drop_quad_poles else []), log_path)
    for fn in (["ef_itp_q{}.pickle".format(i + 1) for i in range(NQ)] +
               ["ef_itp_q{}skew.pickle".format(i + 1) for i in range(NQ) if QALPHA[i] != 0.0 and SKEW[i]] +
               (["basis_meta.json"] if os.path.exists(os.path.join(BASIS0, "basis_meta.json")) else [])):
        shutil.copyfile(os.path.join(BASIS0, fn), os.path.join(bdir, fn))


def bunch_cmd(bdir, sdir, out_dir, tag, superpose, dz_mm, vscale, n, whole, h, resolve_every):
    """bunch.py arguments of one run of a built geometry at a quad setting; whole = draw from the whole RFQ file (core + tail)."""
    return (["-m", "spyral_inflector.tracking.bunch", "--tag", "spiral", "--reload-dir", bdir, "--step-dir", sdir,
             "--particles", PARTICLES, "--bfield", BFIELD, "--phi", PHI, "--out-dir", out_dir, "--out-tag", tag,
             "--basis-dir", bdir, "--superpose"] + list(superpose) + ["--current-ma", a.current_ma, "--handoff-frame", "machine",
             "--handoff-distance", 0.0, "--n", n] + (["--stragglers", "--sc-min-particles", a.sc_min_particles] if whole else []) +
            (["--sc"] if a.current_ma > 0 else []) + ["--h", h, "--resolve-every", resolve_every,
             "--shift-z", "{:.10f}".format(1e-3 * dz_mm - dz0), "--vscale", "{:.6f}".format(vscale)])


def evaluate(x):
    x = np.asarray(x, dtype=float)
    key = key_of(x)
    for e in evals:
        if e["key"] == key and e["measure"] is not None:
            # a cached FAILURE is re-evaluated on resume: 2026-09-15 a memory squeeze on the desktop poisoned the GPU
            # host-memory pool and evaluations e041.. failed in seconds; replaying those as failures would be pointless
            return e
    idx = len(evals)
    tag = "e{:03d}".format(idx)
    edir = os.path.join(OUT, tag)
    gdir, sdir, bdir = (os.path.join(edir, d) for d in ("geometry", "steps", "basis"))
    vals = dict(zip(NAMES, [float(v) for v in x]))
    t_ent, t_exit, dz_mm, vscale = (vals[k] for k in ("t_ent", "t_exit", "dz_mm", "vscale"))
    qv = [vals[q] for q in QN]
    SUPER = [v for q, al in zip(qv, QALPHA) for v in (q, al)]   # voltage, turn alpha per quad (alpha held at the seed)
    knobs = geometry_knobs(vals)
    t0 = time.time()
    m, err, t_geo, t_basis = None, None, None, None
    try:
        rmtree_retry(edir, log=stamp)
        os.makedirs(bdir, exist_ok=True)
        build_geometry(gdir, sdir, ROT_BF, ENERGY, knobs=knobs, fix_truncations=(t_ent, t_exit), fix_dz=dz0, maxiter=0,
                       res=a.geo_res, h=a.geo_h, log=lambda s: None)
        t_geo = time.time() - t0
        # the geometry is built at dz0 and the tracking shifts it by (dz - dz0): the conductors sit at dz
        n_steel = steel_vertices(os.path.join(gdir, "state.pickle"), R, 1e-3 * dz_mm) if a.w_fit > 0 else 0
        _g = (jload(os.path.join(gdir, "summary.json")).get("electrode_gap") or {}).get("min_mm")
        spiral_basis(gdir, sdir, bdir, os.path.join(edir, "log_basis_spiral.txt"))
        t_basis = time.time() - t0 - t_geo
        run_cmd(bunch_cmd(bdir, sdir, edir, tag, SUPER, dz_mm, vscale, a.n, WHOLE, a.h, a.resolve_every) + ["--no-plot"],
                os.path.join(edir, "log_bunch.txt"))
        m = measure(edir, tag)
        m["n_steel"] = n_steel
        m["electrode_gap_mm"] = _g
    except Exception as exc:               # a collision / failed build / failed solve counts as a bad point, not a crash
        err = str(exc).splitlines()[-1][:200] if str(exc) else repr(exc)
    if not a.keep_basis and os.path.isdir(bdir):
        rmtree_retry(bdir)
    npz = os.path.join(edir, "bunch_{}.npz".format(tag))
    if os.path.exists(npz) and not a.keep_npz:
        os.remove(npz)
    r = residuals(m)
    wall = time.time() - t0
    e = {"key": key, "index": idx, "tag": tag, "x": [float(v) for v in x], "measure": m, "error": err,
         "residuals": [float(v) for v in r], "objective": float(np.sum(r ** 2)), "wall_s": wall,
         "wall_split_s": {"geometry": t_geo, "basis": t_basis,
                          "tracking": None if (t_geo is None or t_basis is None) else wall - t_geo - t_basis}}
    evals.append(e)
    save_cache()
    if m is None:
        stamp("  {} FAILED ({:.0f} s): {} | {}".format(tag, e["wall_s"], err, " ".join("{}={:g}".format(n, v) for n, v in zip(NAMES, x))))
    else:
        stamp("  {} ({:.0f} s: geometry {:.0f}, basis {:.0f}, tracking {:.0f}): [{}] angle {:.1f} mrad rms at the hand-off (robust {:.1f}; {:.1f} at 55 mm), "
              "T {:.1f} %, on SI electrodes {:.2f} %, centroid {:+.2f} mm / {:+.1f} mrad, proxy {:.1f} mm, dE {:.2f} keV ({:.0f} % low), "
              "phase {:.1f} deg rms (robust {:.1f}), hand-off r {} mm / alpha {} deg, steel {}, electrode gap {} mm, V {} -> f {:.3f} | {}".format(
            tag, e["wall_s"], t_geo, t_basis, e["wall_split_s"]["tracking"],
            "all {:,d} particles incl. {:,d} tail, {:,d} tail through".format(m["n_particles"], m["n_tail"], m["n_tail_through"]) if WHOLE
            else "CORE ONLY - not for reporting",
            1e3 * m.get("zp_h_rms", m["zp_rms"]), 1e3 * m.get("zp_h_robust", 0.0), 1e3 * m["zp_rms"], 100 * m["transmission"],
            100 * m.get("f_lost_si", 0.0), 1e3 * m["z"], 1e3 * m["zp"], 1e3 * m.get("z100_rms", 0.0),
            m.get("dT_rms_keV", 0.0), 100 * m.get("f_low_energy", 0.0), m.get("dt_rms_deg", 0.0), m.get("dt_robust_deg", 0.0),
            "-" if m.get("r_h_mm") is None else "{:.1f}".format(m["r_h_mm"]),
            "-" if m.get("alpha_deg") is None else "{:+.1f}".format(m["alpha_deg"]), m.get("n_steel", "-"),
            "-" if m.get("electrode_gap_mm") is None else "{:.2f}".format(m["electrode_gap_mm"]),
            "-" if m["voltage_V"] is None else "{:.0f}".format(m["voltage_V"]), e["objective"],
            " ".join("{}={:g}".format(n, v) for n, v in zip(NAMES, x))))
    return e


# ---------------------------------------------------------------- seed, floor, DFO-LS
seed = evaluate(x0)
if seed["measure"] is None:
    raise SystemExit("the seed evaluation failed: " + str(seed["error"]))
if T_FLOOR is None and a.interim:
    T_FLOOR = seed["measure"]["transmission"] - 0.05     # the cached objectives are the running optimizer's; nothing is saved
elif T_FLOOR is None:
    T_FLOOR = seed["measure"]["transmission"] - 0.05
    for e in evals:                                      # residuals of cached points with the final floor
        e["residuals"] = [float(v) for v in residuals(e["measure"])]
        e["objective"] = float(np.sum(np.array(e["residuals"]) ** 2))
    save_cache()
stamp("  transmission floor {:.1f} % (loop scale); seed objective {:.3f} = {}".format(
    100 * T_FLOOR, seed["objective"], " + ".join("{:.1f} ({})".format(r ** 2, nm) for r, nm in zip(seed["residuals"], _res_names))))

soln = None
if a.seed_only or a.final_only or a.interim:
    stamp("--{}: no DFO-LS".format("seed-only" if a.seed_only else "final-only" if a.final_only else "interim"))
else:
    import dfols  # noqa: E402

    def objfun(xf):
        return np.array(evaluate(expand(xf))["residuals"], dtype=float)

    stamp("DFO-LS start: {} free knobs ({}), maxfun {}{}".format(
        len(FREE), ", ".join(NAMES[i] for i in FREE), a.maxfun,
        "; held at the seed: " + ", ".join("{} {:g} {}".format(n, x0[NAMES.index(n)], UNIT_OF[n]) for n in HELD) if HELD else ""))
    soln = dfols.solve(objfun, x0[FREE], bounds=(LO[FREE], HI[FREE]), rhobeg=a.rhobeg, rhoend=a.rhoend,
                       maxfun=a.maxfun, objfun_has_noise=True, scaling_within_bounds=True)
    stamp("DFO-LS done: flag {} ({}), {} evaluations, f {:.3f}".format(soln.flag, soln.msg, soln.nf, float(np.sum(soln.resid ** 2))))

# ---------------------------------------------------------------- report and best knobs (from the cache, best objective)
good = [e for e in evals if e["measure"] is not None]
best = min(good, key=lambda e: e["objective"])
# --interim writes only into <out>/interim_<tag>/: the running optimizer keeps its own report, best_knobs.json and final/
REP_DIR = os.path.join(OUT, "interim_" + best["tag"]) if a.interim else OUT
REP_MD = os.path.join(REP_DIR, "report.md") if a.interim else os.path.join(RUN, a.out_name + ".md")
os.makedirs(REP_DIR, exist_ok=True)
xb = best["x"]
vb = dict(zip(NAMES, xb))
bk = geometry_knobs(vb)
bk.update({"_" + q: vb[q] for q in QN})
bk.update({"_alpha{}".format(i + 1): al for i, al in enumerate(QALPHA)})
bk["_vscale"] = vb["vscale"]
bk.update({"_truncations_deg": [vb["t_ent"], vb["t_exit"]], "_dz_m": 1e-3 * vb["dz_mm"], "_objective": best["objective"],
           "_measure": best["measure"], "_seed_objective": seed["objective"], "_R_deg": R, "_tag": best["tag"], "_names": NAMES,
           "_whole_file": WHOLE})
with open(os.path.join(REP_DIR, "best_knobs.json"), "w") as fh:
    json.dump(bk, fh, indent=2)
with open(os.path.join(REP_DIR, "best_quads.json"), "w") as fh:     # the quad setting in SCLevel's / ABRect's --best-json format
    json.dump({"best": dict({q: vb[q] for q in QN}, **{"alpha{}".format(i + 1): al for i, al in enumerate(QALPHA)}),
               "note": "FullOpt {} winner {} (quad turns alpha held at the seed's)".format(a.out_name, best["tag"])}, fh, indent=1)
who = ("{} particles drawn from the WHOLE RFQ file (core + tail; every number over all of them)".format(a.n) if WHOLE
       else "{} CORE particles (old loop; confirm with an every-particle run before reporting)".format(a.n))
L = ["# Full DFO-LS optimization of {} {}".format(
    os.path.basename(RUN),
    "with {} mA space charge in the loop".format(a.current_ma) if a.current_ma > 0 else "with no space charge"), "",
     "Knobs: " + ", ".join(NAMES) + " (profile nodes ent/mid/exit: entrance, middle and exit section of the truncated electrode, linear "
     "between; the design voltage scales with the gap; the operating voltage is vscale x the inner centering solve). Rotation R = {:+.3f} deg "
     "fixed; quad turns alpha {} held; spiral basis re-solved per evaluation on the run's {:.1f} mm grid, quad bases reused; {}, {} mA, "
     "{:.0f} mm cells. Residuals: {}; floor {:.1f} %. {}.".format(
         R, QALPHA, 1e3 * GRID_RES, who, a.current_ma, 1e3 * a.h, " | ".join(_res_names), 100 * T_FLOOR,
         "DFO-LS: {} evaluations, flag {} ({})".format(soln.nf, soln.flag, soln.msg) if soln is not None
         else ("Seed only (--seed-only)" if a.seed_only else "No DFO-LS in this call (--final-only)")), "",
     "| | " + " | ".join(NAMES) + " | angle at hand-off [mrad rms (robust)] | angle at 55 mm | T | phase [deg rms] | electrode gap [mm] | centroid [mm / mrad] | f |",
     "|---|" + "---|" * (len(NAMES) + 7)]


def _gap_txt(m):
    return "-" if m.get("electrode_gap_mm") is None else "{:.2f}".format(m["electrode_gap_mm"])


for label, e in (("seed", seed), ("best", best)):
    m = e["measure"]
    L.append("| {} | ".format(label) + " | ".join("{:g}".format(v) for v in e["x"]) +
             " | {:.1f} ({:.1f}) | {:.1f} | {:.1f} % | {:.1f} | {} | {:+.2f} / {:+.1f} | {:.3f} |".format(
                 1e3 * m.get("zp_h_rms", m["zp_rms"]), 1e3 * m.get("zp_h_robust", 0.0), 1e3 * m["zp_rms"], 100 * m["transmission"],
                 m.get("dt_rms_deg", 0.0), _gap_txt(m), 1e3 * m["z"], 1e3 * m["zp"], e["objective"]))
L += ["", "## All evaluations (sorted by objective)", "",
      "| tag | " + " | ".join(NAMES) + " | angle hand-off rms (robust) | angle 55 mm | T | phase | gap | height | z' | f | wall [s] |",
      "|---|" + "---|" * (len(NAMES) + 9)]
for e in sorted(evals, key=lambda e: e["objective"]):
    m = e["measure"]
    L.append("| {} | ".format(e["tag"]) + " | ".join("{:g}".format(v) for v in e["x"]) +
             (" | {:.1f} ({:.1f}) | {:.1f} | {:.1f} % | {:.1f} | {} | {:+.2f} | {:+.1f} | {:.3f} | {:.0f} |".format(
                 1e3 * m.get("zp_h_rms", m["zp_rms"]), 1e3 * m.get("zp_h_robust", 0.0), 1e3 * m["zp_rms"], 100 * m["transmission"],
                 m.get("dt_rms_deg", 0.0), _gap_txt(m), 1e3 * m["z"], 1e3 * m["zp"], e["objective"], e.get("wall_s", 0.0))
              if m else " | FAILED | | | | | | | {} | |".format(e["error"])))
L += ["", "The design exit radius is not in the objective; it follows k' (about -0.6 mm per degree). "
      "`FullOptRadius.py --run <this run>` tabulates it per evaluation -> `full_opt_radius.md`."]
L += ["", "Best knobs for the final push: `{}/best_knobs.json` (RunFinalPush --knobs-json it --fix-truncations {:.3f} {:.3f} --fix-dz {:.6f} "
      "{} --fine-step 0; quad turns alpha {} on top of the basis rotations {}), then the plate-rule R iterates to its fixed point there; "
      "level with SCLevel --vscale {:.4f}.".format(
          a.out_name, vb["t_ent"], vb["t_exit"], 1e-3 * vb["dz_mm"], " ".join("--{} {:.0f} {:.0f} 1".format(q, vb[q], vb[q]) for q in QN),
          QALPHA, QROT, vb["vscale"])]
with open(REP_MD, "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
stamp("report -> {}; best {} f {:.3f} vs seed {:.3f}".format(REP_MD, best["tag"], best["objective"], seed["objective"]))

# ---------------------------------------------------------------- the winner with EVERY particle
# Daniel, 2026-10-21: "only inject the core particles during optimization, re-run the winner with the full dist at the end".
# The winner's geometry and STEP files are kept in its evaluation folder; its spiral basis is re-solved the same way, and one
# bunch run tracks all particles of the RFQ file (core + tail injected at its arrival time) with SCLevel --final's space-charge
# settings. Not re-levelled and R not re-solved: the loop's dz / t_exit knobs and its centroid residuals are the levelling here.
if a.final or a.final_only or a.interim:
    ftag = "sc_all_final" if a.current_ma > 0 else "vac_all_final"
    edir = os.path.join(OUT, best["tag"])
    gdir, sdir = os.path.join(edir, "geometry"), os.path.join(edir, "steps")
    fdir = os.path.join(REP_DIR, "final")
    fbdir = os.path.join(fdir, "basis")
    os.makedirs(os.path.join(fdir, "geometry"), exist_ok=True)
    shutil.copyfile(os.path.join(gdir, "voltages.csv"), os.path.join(fdir, "geometry", "voltages.csv"))   # beam_report's pole signs
    stamp("FINAL: the winner {} with EVERY particle of the RFQ file (core + tail) at {:g} mA, {:.1f} mm cells, Poisson solve every {} "
          "steps -> {}".format(best["tag"], a.current_ma, 1e3 * a.final_h, a.final_resolve_every, fdir))
    t0 = time.time()
    if not os.path.exists(os.path.join(fbdir, "ef_itp_spiral.pickle")):
        os.makedirs(fbdir, exist_ok=True)
        spiral_basis(gdir, sdir, fbdir, os.path.join(fdir, "log_basis_spiral.txt"))
    stamp("  spiral basis of {} re-solved ({:.0f} s)".format(best["tag"], time.time() - t0))
    superpose = [v for q, al in zip([vb[q] for q in QN], QALPHA) for v in (q, al)]
    run_cmd(bunch_cmd(fbdir, sdir, fdir, ftag, superpose, vb["dz_mm"], vb["vscale"], a.final_n, True, a.final_h, a.final_resolve_every) +
            ["--save-openpmd", os.path.join(fdir, "handoff_{}.h5".format(ftag))], os.path.join(fdir, "log_{}.txt".format(ftag)))
    mf = measure(fdir, ftag, core_only=False)
    mf["electrode_gap_mm"] = best["measure"].get("electrode_gap_mm")
    ref_tag = "{}_all_level".format("sc" if a.current_ma > 0 else "vac")
    ref_dir = os.path.join(RUN, "sc_level")
    mref = None
    if os.path.exists(os.path.join(ref_dir, "bunch_{}.npz".format(ref_tag))):
        mref = measure(ref_dir, ref_tag, core_only=False)
        mref["electrode_gap_mm"] = (geo0.get("electrode_gap") or {}).get("min_mm")
    with open(os.path.join(fdir, "final.json"), "w") as fh:
        json.dump({"tag": best["tag"], "names": NAMES, "x": best["x"], "measure": mf, "reference": mref, "reference_npz":
                   os.path.join(ref_dir, "bunch_{}.npz".format(ref_tag)), "npz": os.path.join(fdir, "bunch_{}.npz".format(ftag)),
                   "wall_s": time.time() - t0}, fh, indent=1)
    if not a.keep_basis:
        rmtree_retry(fbdir)

    def _row(label, m):
        return "| {} | {:.2f} % | {:.1f} ({:.1f}) | {:.1f} | {:.1f} ({:.1f}) | {:.2f} | {:.2f} % | {:+.2f} / {:+.1f} | {} |".format(
            label, 100 * m["transmission"], 1e3 * m["zp_h_rms"], 1e3 * m["zp_h_robust"], 1e3 * m["zp_rms"], m["dt_rms_deg"],
            m["dt_robust_deg"], m["dT_rms_keV"], 100 * m["f_lost_si"], 1e3 * m["z"], 1e3 * m["zp"], _gap_txt(m))
    L += ["", "## The winner {} with every particle (core + tail), {:g} mA".format(best["tag"], a.current_ma), "",
          "All {:,d} particles of the RFQ file ({:,d} tail, injected at their arrival time), {:.1f} mm cells, Poisson solve every {} steps; "
          "the winner's geometry at the reference's R, not re-levelled. Files: `{}/final/`.".format(
              mf["n_particles"], mf["n_tail"], 1e3 * a.final_h, a.final_resolve_every, a.out_name), "",
          "| | T | angle at hand-off [mrad rms (robust)] | angle at 55 mm | phase [deg rms (robust)] | dE [keV] | on SI electrodes | "
          "centroid [mm / mrad] | electrode gap [mm] |", "|---|" + "---|" * 8]
    if mref is not None:
        L.append(_row("{} (reference, levelled)".format(os.path.basename(RUN)), mref))
    L.append(_row("winner {}".format(best["tag"]), mf))
    with open(REP_MD, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    stamp("FINAL {} ({:.0f} min), every particle: T {:.2f} %, angle {:.1f} mrad rms (robust {:.1f}) at the hand-off, {:.1f} at 55 mm, "
          "phase {:.1f} deg rms, on SI {:.2f} %, centroid {:+.2f} mm / {:+.1f} mrad, hand-off r {} mm / alpha {} deg{}".format(
              best["tag"], (time.time() - t0) / 60, 100 * mf["transmission"], 1e3 * mf["zp_h_rms"], 1e3 * mf["zp_h_robust"], 1e3 * mf["zp_rms"],
              mf["dt_rms_deg"], 100 * mf["f_lost_si"], 1e3 * mf["z"], 1e3 * mf["zp"],
              "-" if mf.get("r_h_mm") is None else "{:.1f}".format(mf["r_h_mm"]),
              "-" if mf.get("alpha_deg") is None else "{:+.1f}".format(mf["alpha_deg"]),
              "" if mref is None else " | reference: T {:.2f} %, {:.1f} mrad, {:.1f} deg".format(
                  100 * mref["transmission"], 1e3 * mref["zp_h_rms"], mref["dt_rms_deg"])))
