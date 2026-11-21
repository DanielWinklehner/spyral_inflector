"""Particle files at the start plane and halfway through the SI entrance aperture for ANY levelled run (2026-10-21).

RunProbe.ps1 did this for t1nr_circ_c2 with its plane and setting written in by hand (Daniel, 2026-10-20). Here they come
from the run itself:
  * the setting = the run's final every-particle run <run>/<level>/bunch_<tag>.json (basis and STEP folders, phi, levelling
    shift, quads incl. turns, vscale, space-charge cell and cadence, particles) and the levelling's B map;
  * the plane = the mid-plane of the Entrance_Aperture plate of that run's STEP folder (deck frame), moved by the assembly's
    axial shifts the tracking applies (the basis state's shift_all_m + the levelling shift), normal to z.
One every-particle bunch run with --probe-z into <run>/probe, then probe_files.py writes <run>/probe/files/
<run>_initial_particles.txt and <run>_probe_z<plane>mm.txt (+ _index.txt): x y z vx vy vz phase, Baseline frame, m, m/s, rad.

    python probe_run.py <run> [--level sc_level] [--from-tag sc_all_level] [--dry-run]
"""
import argparse
import glob
import json
import os
import pickle
import subprocess
import sys
import time


def main(argv=None):
    """Particle files at the start plane and in the entrance aperture of a levelled run."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run")
    p.add_argument("--level", default="sc_level")
    p.add_argument("--from-tag", default="sc_all_level", help="the final run whose setting is repeated")
    p.add_argument("--tag", default="sc_all_probe")
    p.add_argument("--dry-run", action="store_true", help="print the plane and the bunch command, run nothing")
    p.add_argument("--sc-h", type=float, default=None, help="space-charge cell [m] instead of the --from-tag run's (e.g. 0.002 = SCLevel --final)")
    p.add_argument("--resolve-every", type=int, default=None, help="Poisson cadence instead of the --from-tag run's (8 = SCLevel --final)")
    a = p.parse_args(argv)

    SCRIPTS = os.path.dirname(os.path.abspath(__file__))
    RUN = os.path.abspath(a.run)
    NAME = os.path.basename(RUN)
    LVL = os.path.join(RUN, a.level)
    OUT = os.path.join(RUN, "probe")
    js = json.load(open(os.path.join(LVL, "bunch_{}.json".format(a.from_tag))))
    lv = json.load(open(os.path.join(LVL, "sc_level.json")))
    sp = js["superpose"]
    basis, steps = sp["basis_dir"], js["step_dir"]
    quads = [v for pair in sp["quads"] for v in pair]
    sc = js.get("sc") or {}
    with open(os.path.join(basis, "si_state_spiral.pickle"), "rb") as fh:
        shift_all = float(pickle.load(fh).get("shift_all_m", 0.0) or 0.0)
    shift_z = float(js.get("shift_z_m") or 0.0)

    # the entrance aperture plate's mid-plane, tracking frame
    import gmsh  # noqa: E402
    ap = glob.glob(os.path.join(steps, "*Entrance_Aperture.stp"))
    if len(ap) != 1:
        raise SystemExit("expected one Entrance_Aperture STEP in {}, found {}".format(steps, ap))
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        vols = gmsh.model.occ.importShapes(ap[0])
        gmsh.model.occ.synchronize()
        zlo = min(gmsh.model.getBoundingBox(d_, t_)[2] for d_, t_ in vols)
        zhi = max(gmsh.model.getBoundingBox(d_, t_)[5] for d_, t_ in vols)
    finally:
        gmsh.finalize()
    z_plane = 0.5 * (zlo + zhi) + shift_all + shift_z
    print("{}: entrance aperture plate z {:+.2f} .. {:+.2f} mm in its STEP (deck frame), shifts {:+.3f} (state) {:+.3f} (levelling) mm "
          "-> probe plane deck z {:+.5f} m = Baseline {:+.2f} mm".format(NAME, 1e3 * zlo, 1e3 * zhi, 1e3 * shift_all, 1e3 * shift_z,
                                                                         z_plane, -1e3 * z_plane), flush=True)
    cmd = ["-m", "spyral_inflector.tracking.bunch", "--tag", js.get("tag", "spiral"), "--reload-dir", basis, "--step-dir", steps,
           "--particles", js["particles"], "--bfield", lv["bfield"], "--phi={:.10f}".format(float(js["phi_deg"])), "--out-dir", OUT,
           "--out-tag", a.tag, "--basis-dir", basis, "--superpose"] + ["{:.10g}".format(v) for v in quads] + [
           "--current-ma", "{:g}".format(float(sc.get("current_mA", 0.0))), "--vscale", "{:.6f}".format(float(js.get("vscale") or 1.0)),
           "--handoff-frame", "machine", "--handoff-distance", "0.0", "--n", "10000000", "--stragglers",
           "--shift-z={:.10f}".format(shift_z), "--probe-z={:.6f}".format(z_plane)]
    if js.get("space_charge"):
        cmd += ["--sc", "--h", "{:g}".format(a.sc_h or float(sc["h_m"])), "--resolve-every", str(a.resolve_every or int(sc["resolve_every"])),
                "--sc-min-particles", "200"]
    print("python " + " ".join(cmd), flush=True)
    if a.dry_run:
        sys.exit(0)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(OUT, "log_{}.txt".format(a.tag)), "w", encoding="utf-8") as fh:
        r = subprocess.run([sys.executable, "-u"] + cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=SCRIPTS)
    if r.returncode != 0:
        raise SystemExit("bunch failed (exit {}), see {}".format(r.returncode, os.path.join(OUT, "log_{}.txt".format(a.tag))))
    pj = json.load(open(os.path.join(OUT, "bunch_{}.json".format(a.tag))))
    print("probe run done in {:.0f} min: transmitted {:.2f} % (the final run {:.2f} %)".format(
        (time.time() - t0) / 60, 100 * pj["transmission"], 100 * js["transmission"]), flush=True)
    r = subprocess.run([sys.executable, "-u", os.path.join(SCRIPTS, "probe_files.py"), os.path.join(OUT, "bunch_{}.npz".format(a.tag)),
                        "--out", os.path.join(OUT, "files"), "--prefix", NAME], cwd=SCRIPTS)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
