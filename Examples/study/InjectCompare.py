"""Compare the two ways of presenting the RFQ core to the line, on a finished final-push run:

  spatial   the core is a frozen bunch at t = 0, its phase turned into a longitudinal spread
            around the start plane (the convention of every run so far). About a third of the
            core is then born level with or inside the bore of quad 1's entrance plate, and
            the head starts 16 mm downstream of the plane, so it never traverses the field it
            has formally already passed.
  injected  every particle crosses the start plane at its own arrival time t = phase/(360 f),
            the convention of the inflector -> central-region hand-off (bunch --inject-core).
            Nothing is born inside an electrode; the price is that the charge still upstream
            of the plane is missing from the Poisson solve during the injection ramp, which
            mostly removes longitudinal push (see the TODO in tracking/bunch.py).

Runs the injected convention for the run's vacuum and space-charge settings and compares each
with the run's existing spatial result. The vacuum pair isolates the initial-condition
artifact; the 8 mA pair adds the ramp.

    python InjectCompare.py --run <deck>\\Results\\final\\vfocus1_level

Writes <run>/inject_compare/ (the runs) and <run>/inject_compare.md.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

DECK = os.environ.get("SI_DECK", r"D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron\Spiral_inflector")
SCRIPTS = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run", required=True)
p.add_argument("--tag", default="fine_best", help="basis tag of the run's tracking field")
p.add_argument("--n", type=int, default=10 ** 7, help="particles (at or above the file size = all of them)")
p.add_argument("--current-ma", type=float, default=8.0)
p.add_argument("--h", type=float, default=0.002)
p.add_argument("--resolve-every", type=int, default=8)
p.add_argument("--sc-min-particles", type=int, default=200)
p.add_argument("--cases", nargs="+", default=["vac", "sc"], help="which pairs to run")
p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "final1_baseline_1mm.pickle"))
a = p.parse_args()

RUN = os.path.abspath(a.run)
BASIS, STEPS = os.path.join(RUN, "basis"), os.path.join(RUN, "steps")
OUT = os.path.join(RUN, "inject_compare")
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


R = jload(os.path.join(RUN, "exit_plane.json"))["rotation_deg"]
PHI = 90.0 + R
stamp("INJECT COMPARE on {}: R {:+.3f} deg, phi {:+.3f} deg, {} particles, {} mA, {:.0f} mm cells".format(
    os.path.basename(RUN), R, PHI, a.n, a.current_ma, 1e3 * a.h))


def run_injected(case):
    """The injected-convention twin of the run's <case>_all result."""
    tag = "{}_all_inject".format(case)
    if os.path.exists(os.path.join(OUT, "bunch_{}.json".format(tag))):
        stamp("  {} exists".format(tag))
        return tag
    cmd = ["-m", "spyral_inflector.tracking.bunch", "--tag", a.tag, "--reload-dir", BASIS, "--step-dir", STEPS,
           "--particles", a.particles, "--bfield", a.bfield, "--phi", PHI, "--out-dir", OUT, "--out-tag", tag,
           "--current-ma", a.current_ma, "--handoff-frame", "machine", "--handoff-distance", 0.0,
           "--n", a.n, "--stragglers", "--inject-core", "--no-plot"]
    if case == "sc":
        cmd += ["--sc", "--h", a.h, "--resolve-every", a.resolve_every, "--sc-min-particles", a.sc_min_particles]
    t0 = time.time()
    stamp("  running {} (injected convention)".format(tag))
    with open(os.path.join(OUT, "log_{}.txt".format(tag)), "w", encoding="utf-8") as fh:
        rc = subprocess.call([PY, "-u"] + [str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS)
    if rc != 0:
        raise SystemExit("{} failed (exit {}), see {}".format(tag, rc, os.path.join(OUT, "log_{}.txt".format(tag))))
    stamp("  {} done in {:.0f} min".format(tag, (time.time() - t0) / 60.0))
    return tag


def measure(folder, tag):
    s = jload(os.path.join(folder, "bunch_{}.json".format(tag)))
    d = np.load(os.path.join(folder, "bunch_{}.npz".format(tag)), allow_pickle=True)
    st = d["asym_state"]
    core = ~d["is_tail"] if "is_tail" in d.files else np.ones(len(d["crossed"]), dtype=bool)
    ok = d["crossed"] & np.all(np.isfinite(st), axis=1) & core
    z = -1e3 * st[ok, 2]                               # Baseline height [mm]
    v = st[ok, 3:6]
    zp = -1e3 * v[:, 2] / np.hypot(v[:, 0], v[:, 1])   # Baseline vertical angle [mrad]
    ti = s.get("tail") or {}
    losses = s.get("losses_by_electrode") or {}
    n = s["n_particles"]
    return {"T": s["transmission"], "core": ti.get("core_transmission"), "tail": ti.get("tail_transmission"),
            "z": float(z.mean()), "zp": float(zp.mean()), "z_rms": float(z.std()), "zp_rms": float(zp.std()),
            "ent0": losses.get("ent0", 0) / n, "aper": losses.get("Entrance_Aperture", 0) / n,
            "anode": losses.get("SI_Anode", 0) / n, "n": int(ok.sum())}


rows = []
for case in a.cases:
    base_tag = "{}_all".format(case)
    if not os.path.exists(os.path.join(BASIS, "bunch_{}.json".format(base_tag))):
        stamp("  no baseline {} in {}; skipping {}".format(base_tag, BASIS, case))
        continue
    inj_tag = run_injected(case)
    old, new = measure(BASIS, base_tag), measure(OUT, inj_tag)
    rows.append((case, old, new))
    stamp("  {}: T {:.2f} -> {:.2f} %, angle {:.1f} -> {:.1f} mrad, centroid {:+.2f}/{:+.1f} -> {:+.2f}/{:+.1f}".format(
        case, 100 * old["T"], 100 * new["T"], old["zp_rms"], new["zp_rms"], old["z"], old["zp"], new["z"], new["zp"]))

L = ["# Spatial bunch versus arrival-time injection of the RFQ core: {}".format(os.path.basename(RUN)), "",
     "Same geometry, field, particles and space-charge settings; only the way the core is presented at the start plane "
     "differs. **spatial**: frozen bunch at t = 0, the phase as a longitudinal spread, which puts 31 % of the core level "
     "with or inside the bore of quad 1's entrance plate and its head 16 mm downstream of the plane. **injected**: every "
     "particle crosses the plane at its own arrival time (the hand-off convention), so nothing is born inside an "
     "electrode, at the price of the charge still upstream being absent from the Poisson solve during the ~22 ns ramp.",
     "", "The vacuum pair isolates the initial-condition artifact; the 8 mA pair adds the ramp, which removes mostly "
     "longitudinal push (TODO in tracking/bunch.py: feed an estimate of the uninjected line charge into the AMG solve).",
     "", "| case | convention | transmitted | core | tail | rms vert. angle [mrad] | centroid [mm / mrad] | rms height [mm] | ent0 loss | entrance aperture | SI anode |",
     "|---|---|---|---|---|---|---|---|---|---|---|"]
for case, old, new in rows:
    for label, m in (("spatial", old), ("injected", new)):
        L.append("| {} | {} | {:.2f} % | {} | {} | {:.1f} | {:+.2f} / {:+.1f} | {:.2f} | {:.2f} % | {:.2f} % | {:.2f} % |".format(
            "8 mA" if case == "sc" else "vacuum", label, 100 * m["T"],
            "-" if m["core"] is None else "{:.2f} %".format(100 * m["core"]),
            "-" if m["tail"] is None else "{:.2f} %".format(100 * m["tail"]),
            m["zp_rms"], m["z"], m["zp"], m["z_rms"], 100 * m["ent0"], 100 * m["aper"], 100 * m["anode"]))
if rows:
    L += ["", "## Differences (injected minus spatial)", "",
          "| case | transmission | core | rms vert. angle | centroid height | centroid angle |", "|---|---|---|---|---|---|"]
    for case, old, new in rows:
        L.append("| {} | {:+.2f} points | {} | {:+.1f} mrad | {:+.2f} mm | {:+.1f} mrad |".format(
            "8 mA" if case == "sc" else "vacuum", 100 * (new["T"] - old["T"]),
            "-" if old["core"] is None or new["core"] is None else "{:+.2f} points".format(100 * (new["core"] - old["core"])),
            new["zp_rms"] - old["zp_rms"], new["z"] - old["z"], new["zp"] - old["zp"]))
    L += ["", "Run-to-run noise of these quantities at the full particle count is well below 0.1 points and 0.2 mrad "
          "(the tracker is deterministic; only the convention changed), so any difference here is systematic."]
with open(os.path.join(RUN, "inject_compare.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
stamp("report -> {}".format(os.path.join(RUN, "inject_compare.md")))
