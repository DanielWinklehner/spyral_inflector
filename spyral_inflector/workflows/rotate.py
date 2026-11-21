"""One pass of the CORNER-rule rotation/levelling fixed point (the e016e procedure, rebuilt 2026-10-19).

Daniel, 2026-10-19: re-export t1nr with the CORNER rule -- "rotation fixed point -> bases -> SCLevel ->
ExportLevelled, repeated until |dR| < 0.2 deg. Nothing else changes."

This script does the first two stages for a new run folder, starting from the LEVELLED setting of --prev:
  1. rotation fixed point (RunFinalPush phases 4-5): rotate the B map by R, rebuild the geometry in it with
     --prev's knobs, both truncations fixed at --prev's levelled values and dz fixed at --prev's levelled dz
     (maxiter 0: only the voltage centering solve runs), re-solve R from the corner rule; repeat until the
     rotation moves by less than --r-tol. Writes <out>/bfield_rotated.pickle and <out>/exit_plane.json.
  2. the seven bases at that R, exactly as ABRect.py solves them: the spiral without the quad poles
     (--drop D*), the quads and their skews with every electrode, all with --rotate-all R on the unrotated map.
SCLevel (with --base-run <out>) and ExportLevelled (--rotation-rule corner ...) follow; corner_pass_run.ps1
chains them and starts another pass while the export's re-solved R differs from the run's by --r-tol or more.

    python corner_pass.py --prev E:\\HCHC1.5\\Results\\triplet\\t1nr --name t1nr_c1 --r0 -91.0
    python corner_pass.py --prev E:\\HCHC1.5\\Results\\triplet\\t1nr_c1 --name t1nr_c2      (R0 from t1nr_c1's export)
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time


def main(argv=None):
    """One corner-rule rotation pass: the rotation fixed point and the bases of a new run (command line: see the module docstring)."""
    SCRIPTS = os.path.dirname(os.path.abspath(__file__))
    DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
    PY = sys.executable
    sys.path.insert(0, SCRIPTS)
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--prev", required=True, help="run whose LEVELLED setting (sc_level/sc_level.json) starts this pass")
    p.add_argument("--name", required=True, help="new run folder under <SI_RESULTS>/triplet")
    p.add_argument("--best-knobs", default=None,
                   help="start from a FullOpt winner instead of --prev's levelled setting: its best_knobs.json gives the knobs (incl. "
                        "angling_profile / sigma_profile / aspect), both truncations (_truncations_deg) and dz (_dz_m) -- FullOpt builds "
                        "at the seed's dz and tracks the shift, so _dz_m is the setting the winner was scored at. --prev still gives "
                        "the quad basis turns and the default R0 (2026-10-21, the profile campaign)")
    p.add_argument("--r0", type=float, default=None,
                   help="starting rotation [deg]; default: --prev's corner-rule export (export_baseline_levelled/exit_plane_spec.json)")
    p.add_argument("--gap-azimuth", type=float, default=24.0)
    p.add_argument("--half-gap", type=float, default=0.00325)
    p.add_argument("--corner", choices=("outer", "inner"), default="outer")
    p.add_argument("--corner-offset-mm", type=float, default=3.25)
    p.add_argument("--r-tol", type=float, default=0.2)
    p.add_argument("--max-iter", type=int, default=4)
    p.add_argument("--quads", type=float, nargs="+", default=None, help="basis turns of the quads [deg] (default: --prev's levelling)")
    p.add_argument("--bfield", default=r"E:\HCHC1.5\Fields\hchc15_injax_1mm.pickle", help="the UNROTATED map")
    p.add_argument("--particles", default=os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit.dst"))
    p.add_argument("--res", type=float, default=0.0025, help="basis grid [m]")
    p.add_argument("--h", type=float, default=0.005)
    p.add_argument("--geo-res", type=float, default=0.005)
    p.add_argument("--no-bases", action="store_true")
    p.add_argument("--family", default="triplet", help="results sub-folder of the new run (magquad: the no-ESQ runs, 2026-10-20)")
    a = p.parse_args(argv)

    from spyral_inflector.tracking.deck import particle_rows  # noqa: E402
    from spyral_inflector.tracking.frames import to_deck_frame, rotate_bfield_z  # noqa: E402
    from spyral_inflector.tracking.geometry_point import build_geometry, rmtree_retry  # noqa: E402
    from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402
    from PyPATools.field import Field  # noqa: E402

    RESULTS = os.environ.get("SI_RESULTS") or os.path.join(DECK, "Results")
    PREV = os.path.abspath(a.prev)
    OUT = os.path.join(RESULTS, a.family, a.name)
    GEO, STEPS, BASIS = (os.path.join(OUT, d) for d in ("geometry", "steps", "basis"))
    ROT_BF = os.path.join(OUT, "bfield_rotated.pickle")
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


    s0 = jload(os.path.join(PREV, "geometry", "summary.json"))
    knobs = dict(s0["knobs"])
    t_ent = float(s0["optimizer"]["truncations_deg"][0])
    lv = jload(os.path.join(PREV, "sc_level", "sc_level.json"))
    t_exit, dz = float(lv["t_exit_deg"]), float(lv["dz_m"])
    qrot = [float(x) for x in (a.quads or lv["quad_rotation_deg"])]
    if a.best_knobs:
        bk = jload(a.best_knobs)
        knobs = {k_: v_ for k_, v_ in bk.items() if not k_.startswith("_")}
        t_ent, t_exit = (float(v_) for v_ in bk["_truncations_deg"])
        dz = float(bk["_dz_m"])
    if a.r0 is not None:
        R = float(a.r0)
        r_src = "--r0"
    else:
        sp = jload(os.path.join(PREV, "export_baseline_levelled", "exit_plane_spec.json"))
        if sp.get("rule_name") != "corner" or abs(float(sp["gap_azimuth_deg"]) - a.gap_azimuth) > 1e-9:
            raise SystemExit("{}'s export used rule '{}' at {} deg, not the corner rule at {:g} deg: give --r0".format(
                os.path.basename(PREV), sp.get("rule_name"), sp.get("gap_azimuth_deg"), a.gap_azimuth))
        R = float(sp["rotation_deg"])
        r_src = "{}'s corner export".format(os.path.basename(PREV))
    ENERGY = float(particle_rows(a.particles, core=True)[:, 8].mean())
    stamp("CORNER PASS {} from {}: truncations {:.6f} / {:.6f} deg, dz {:+.4f} mm (the levelled setting, fixed), "
          "R0 {:+.4f} deg ({}); rule corner, {} corner {:g} mm to +x of the {:g} deg ray; quad basis turns {}".format(
              a.name, a.best_knobs or os.path.basename(PREV), t_ent, t_exit, 1e3 * dz, R, r_src, a.corner, a.corner_offset_mm,
              a.gap_azimuth, qrot))
    if a.best_knobs:
        stamp("  knobs from the FullOpt winner: " + ", ".join("{} {}".format(k_, knobs.get(k_)) for k_ in
                                                              ("tilt", "angling_profile", "gamma", "sigma_profile", "aspect", "gap")))

    # 1. rotation fixed point: field rotated by R -> geometry in it -> R from the corner rule
    done_fp = os.path.join(OUT, "fixed_point.json")
    if not os.path.exists(done_fp):
        history = [R]
        fd = to_deck_frame(Field.from_file(a.bfield), log=stamp)
        for it in range(a.max_iter):
            stamp("fixed point {}: the field a system rotated by {:+.4f} deg sees".format(it + 1, R))
            rotate_bfield_z(fd, R, log=stamp).save(ROT_BF)
            rmtree_retry(GEO, log=stamp)
            if os.path.isdir(STEPS):
                shutil.rmtree(STEPS)
            build_geometry(GEO, STEPS, ROT_BF, ENERGY, knobs=knobs, fix_truncations=(t_ent, t_exit), fix_dz=dz, maxiter=0,
                           res=a.geo_res, h=a.h, log=stamp)
            spec = exit_plane_spec(STEPS, os.path.join(GEO, "state.pickle"), a.gap_azimuth, a.half_gap,
                                   out_json=os.path.join(OUT, "exit_plane.json"), log=stamp, rule="corner",
                                   corner_offset_mm=a.corner_offset_mm, corner_which=a.corner, r_near=R)
            R_new = float(spec["rotation_deg"])
            stamp("  rotation {:+.4f} -> {:+.4f} deg (change {:+.4f}); corner offset {:+.4f} mm, voltage {:.1f} V".format(
                R, R_new, R_new - R, spec["corner"]["offset_mm_achieved"], jload(os.path.join(GEO, "summary.json"))["voltage"]))
            history.append(R_new)
            conv = abs(R_new - R) <= a.r_tol
            R_field = R
            R = R_new
            if conv:
                break
        else:
            raise SystemExit("rotation did not settle within {} iterations: {}".format(a.max_iter, history))
        with open(done_fp, "w") as fh:
            json.dump({"prev": PREV, "R_history": history, "R": R, "R_field": R_field, "t_ent_deg": t_ent, "t_exit_deg": t_exit,
                       "dz_m": dz, "quad_rotation_deg": qrot, "rule": "corner", "gap_azimuth_deg": a.gap_azimuth,
                       "corner": a.corner, "corner_offset_mm": a.corner_offset_mm,
                       "voltage": jload(os.path.join(GEO, "summary.json"))["voltage"]}, fh, indent=2)
    fp = jload(done_fp)
    R = float(fp["R"])
    stamp("  R* = {:+.4f} deg (field rotated by {:+.4f}), history {}".format(R, fp["R_field"], ["{:+.4f}".format(x) for x in fp["R_history"]]))
    if a.no_bases:
        raise SystemExit(0)

    # 2. the seven bases at R, as ABRect.py (a run without electrostatic quads, knob n_quads = 0: the spiral basis alone,
    #    with the beam pipe and every other grounded electrode in it, as MagQuadSetup.py solves it)
    os.makedirs(BASIS, exist_ok=True)
    common = ["--step-dir", STEPS, "--voltages", os.path.join(GEO, "voltages.csv"), "--state", os.path.join(GEO, "state.pickle"),
              "--out-dir", BASIS, "--res", a.res, "--h", a.h, "--bfield", a.bfield, "--rotate-all", R, "--energy-mev", ENERGY]
    if knobs.get("n_quads") is not None and int(knobs["n_quads"]) == 0:
        bases = (("spiral", []),)
    else:
        a1, a2, a3 = qrot
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
    stamp("corner pass {}: geometry, exit_plane.json (R {:+.4f}) and the seven bases are in {}".format(a.name, R, OUT))


if __name__ == "__main__":
    main()
