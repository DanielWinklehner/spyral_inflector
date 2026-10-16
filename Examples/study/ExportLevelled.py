"""Baseline-frame export of a space-charge-levelled geometry (follow-up 4 of the vfocus campaign).

SCLevel.py leaves the levelled geometry as the STEP files of <run>/sc_level/t<trunc>/steps (rebuilt at
the levelled exit truncation, at the run's dz) plus an axial shift applied at tracking time. The
plate rule that fixes the system rotation R depends only on the x-y position of the beam's crossing
of the exit plate's outer face, so the shift leaves R untouched, but the truncation moves the exit
face crossing slightly: R is re-solved here on the levelled STEP set and its design orbit, then
every electrode is rotated by that R (quads by their extra angles), shifted by the levelled dz in
the deck frame, mirrored into the machine frame and written in millimetres. The levelled hand-off
files (already tracked with the shift) are copied alongside.

    python ExportLevelled.py --run <deck>\\Results\\final\\vfocus1_final
Writes <run>/export_baseline_levelled/ (steps_baseline_mm/, exit_plane_spec.json, handoff_*_level.h5, README.txt).
"""
import argparse
import json
import os
import shutil
import sys
import time

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(SCRIPTS)))

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run", required=True)
p.add_argument("--level", default="sc_level", help="levelling folder inside the run (with sc_level.json)")
p.add_argument("--quads", type=float, nargs="+", default=[16.0, 24.0], help="quad rotations on top of R [deg], one per quad (2 for a doublet, 3 for a triplet)")
p.add_argument("--rotation-rule", choices=("offset", "parallel", "corner"), default="offset",
               help="which plate rule re-solves R on the levelled STEP set. It MUST match the rule the run was built with -- ExportLevelled had no such flag until 2026-09-19 and silently used offset, which on the parallel-rule e016b run returned -87.0 deg instead of -61.4 and left the plate 25.6 deg off parallel, with every STEP file and the E-field rotated by it.")
p.add_argument("--corner-offset-mm", type=float, default=3.25,
               help="rule=corner: how far to the +x side of the gap ray the chosen corner sits")
p.add_argument("--corner", choices=("outer", "inner"), default="outer",
               help="rule=corner: which end of the exit plate is pinned")
p.add_argument("--r-near", type=float, default=None,
               help="rule=corner: pick the branch nearest this rotation (default: the run's own R)")
p.add_argument("--gap-azimuth", type=float, default=34.0)
p.add_argument("--half-gap", type=float, default=0.005)
p.add_argument("--out", default=None, help="default: <run>/export_baseline_levelled")
a = p.parse_args()

import numpy as np  # noqa: E402

from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402
from spyral_inflector.tracking.export import (HOUSING_OUTLINE_NAME, steps_to_baseline_mm,  # noqa: E402
                                              field_to_machine_frame, housing_outline_to_baseline_mm)
from spyral_inflector.tracking.deck import superpose_basis  # noqa: E402
from PyPATools.field import Field  # noqa: E402

RUN = os.path.abspath(a.run)
LVL = os.path.join(RUN, a.level)
OUT = a.out or os.path.join(RUN, "export_baseline_levelled")
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log_export.txt"), "a", encoding="utf-8")


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


with open(os.path.join(LVL, "sc_level.json")) as fh:
    lv = json.load(fh)
if not lv.get("converged", False):
    stamp("WARNING: the levelling of {} did not converge; exporting its last setting anyway".format(os.path.basename(RUN)))
t_exit, dz, dz0, shift = lv["t_exit_deg"], lv["dz_m"], lv["dz0_m"], lv["shift_z_m"]
steps = os.path.join(LVL, "t{:.3f}".format(t_exit), "steps")
state = os.path.join(LVL, "t{:.3f}".format(t_exit), "geometry", "state.pickle")
if not os.path.isdir(steps) or not os.path.exists(state):
    # SCLevel only builds a t<trunc>/ folder when the levelling MOVED the truncation; when it
    # converges at the run's own value -- which is what a fully converged rotation/levelling
    # fixed point looks like -- the levelled geometry simply IS the run's own.
    run_steps = os.path.join(RUN, "steps")
    run_state = os.path.join(RUN, "geometry", "state.pickle")
    t_run = float(json.load(open(os.path.join(RUN, "geometry", "summary.json")))["optimizer"]["truncations_deg"][1])
    if abs(t_exit - t_run) < 1e-6 and os.path.isdir(run_steps) and os.path.exists(run_state):
        stamp("levelling did not move the truncation ({:.4f} deg): the run's own geometry IS the "
              "levelled one".format(t_exit))
        steps, state = run_steps, run_state
    else:
        raise SystemExit("levelled geometry not found: {} / {} (run truncation {:.4f} deg)".format(
            steps, state, t_run))
with open(os.path.join(RUN, "exit_plane.json")) as fh:
    R_run = json.load(fh)["rotation_deg"]
stamp("EXPORT LEVELLED {}: exit truncation {:.3f} deg, dz {:+.2f} mm (run {:+.2f}, shift {:+.2f} mm), run R {:+.3f} deg".format(
    os.path.basename(RUN), t_exit, 1e3 * dz, 1e3 * dz0, 1e3 * shift, R_run))

# 1. the plate rule on the levelled STEP set (z-shift invariant; the truncation moves the crossing a little)
if a.r_near is not None:
    _r_near = a.r_near
else:                                  # the run's own rotation, so the branch does not flip
    with open(os.path.join(RUN, "exit_plane.json")) as _fh:
        _r_near = float(json.load(_fh)["rotation_deg"])
spec = exit_plane_spec(steps, state, a.gap_azimuth, a.half_gap, out_json=os.path.join(OUT, "exit_plane_spec.json"),
                       log=stamp, rule=a.rotation_rule, corner_offset_mm=a.corner_offset_mm,
                       corner_which=a.corner, r_near=_r_near)
R = spec["rotation_deg"]
stamp("  R re-solved on the levelled geometry: {:+.3f} deg (run {:+.3f}, change {:+.3f})".format(R, R_run, R - R_run))
if abs(R - R_run) > 0.2:
    stamp("  NOTE: |dR| > 0.2 deg -- the levelled geometry was rebuilt in the field rotated by the run's R; a re-optimization "
          "at the new R (RunFinalPush phase 5 loop) would be self-consistent")

# 2. STEP export: rotate by R (quads on top), shift by the levelled dz in the deck frame, mirror, millimetres
files = steps_to_baseline_mm(steps, os.path.join(OUT, "steps_baseline_mm"), rotation_deg=R, quad_rotation=tuple(a.quads),
                             combined="{}_inflector_{}_levelled_baseline_mm.stp".format(os.environ.get("SI_MACHINE", "HCHC60"),
                                                                                         os.path.basename(RUN)),
                             shift_z=shift, log=stamp)

# 2b. the housing OUTLINE through the same transform: the polygon the wall is extruded from and the four
#     corners of the exit opening. The central region builds its own housing from these (no top and bottom
#     plates, open radial ends) instead of importing 002_Housing.stp.
ol_src = os.path.join(steps, HOUSING_OUTLINE_NAME)
has_outline = os.path.exists(ol_src)
if has_outline:
    housing_outline_to_baseline_mm(ol_src, os.path.join(OUT, "steps_baseline_mm", HOUSING_OUTLINE_NAME),
                                   rotation_deg=R, shift_z=shift,
                                   comments=("levelled export of {}".format(os.path.basename(RUN)),))
    stamp("  housing outline -> steps_baseline_mm/{}".format(HOUSING_OUTLINE_NAME))
else:
    stamp("  NOTE: no {} in the levelled STEP folder (built before the outline was exported); "
          "the central region has only 002_Housing.stp".format(HOUSING_OUTLINE_NAME))

# 3. the levelled E-field (Daniel, 2026-09-15 20:30: "add the field and step export"): the levelled spiral basis
#    (re-solved by SCLevel at the levelled truncation) plus the quads at the levelled voltages and rotations, times the
#    operating spiral voltage scale, translated by the levelling shift exactly as bunch.py tracks it, then mirrored into
#    the machine frame. Same file name as RunFinalPush phase 10, so the central-region scripts take it via
#    HCHC60_EFIELD with HCHC60_EFIELD_ZSHIFT_MM = 0 (the shift is already in the field).
basis_dir = lv.get("basis") or ""
if not os.path.isdir(basis_dir):                        # a levelling done on another machine records that machine's path
    basis_dir = os.path.join(LVL, "t{:.3f}".format(t_exit), "basis")
if not os.path.isdir(basis_dir):
    basis_dir = os.path.join(RUN, "basis")
    stamp("  NOTE: no levelled basis folder found; the field is superposed from the run's own bases (unlevelled truncation)")
quads_lv = [float(v) for v in lv["quads"]]
vscale = float(lv.get("vscale", 1.0))
e_lvl = superpose_basis(basis_dir, *quads_lv, vscale=vscale)
g = e_lvl.grid
e_lvl = Field.from_arrays(grid={"x": np.asarray(g["x"], dtype=float), "y": np.asarray(g["y"], dtype=float),
                                "z": np.asarray(g["z"], dtype=float) + shift},
                          values={c: np.asarray(e_lvl.grid_values[c]) for c in "xyz"}, dim=3, units="m",
                          label="{} shifted {:+.2f} mm (levelled)".format(getattr(e_lvl, "label", "E-field"), 1e3 * shift))
field_to_machine_frame(e_lvl, os.path.join(OUT, "efield_fine_best_baseline.pickle"),
                       label="levelled E-field of {} (machine frame)".format(os.path.basename(RUN)))
stamp("  E-field: {} bases from {} at quads {} x vscale {:.4f}, shifted {:+.3f} mm, grid {} x {} x {} -> efield_fine_best_baseline.pickle".format(
    1 + len(quads_lv) // 2, basis_dir, [round(q) for q in quads_lv[0::2]], vscale, 1e3 * shift, len(g["x"]), len(g["y"]), len(g["z"])))
with open(os.path.join(OUT, "levelled_export.json"), "w") as fh:
    json.dump({"run": RUN, "level": LVL, "basis": basis_dir, "quads": quads_lv, "quad_rotation_deg": list(a.quads), "vscale": vscale,
               "t_exit_deg": t_exit, "dz_m": dz, "dz0_m": dz0, "shift_z_m": shift, "rotation_deg": R, "rotation_run_deg": R_run,
               "efield": "efield_fine_best_baseline.pickle", "efield_zshift_mm": 0.0, "housing": "steps_baseline_mm/002_Housing.stp",
               "housing_outline": ("steps_baseline_mm/" + HOUSING_OUTLINE_NAME) if has_outline else None}, fh, indent=2)
del e_lvl

# 4. hand-offs of the levelled full runs (tracked with the shift, machine frame at the electrode exit)
copied = []
for fn in ("handoff_vac_all_level.h5", "handoff_sc_all_level.h5"):
    src = os.path.join(LVL, fn)
    if os.path.exists(src):
        shutil.copyfile(src, os.path.join(OUT, fn))
        copied.append(fn)
with open(os.path.join(OUT, "README.txt"), "w", encoding="utf-8") as fh:
    fh.write("Baseline (machine) frame: right-handed, +z up, the beam enters from +z; azimuth counter-clockwise from above.\n"
             "Space-charge-LEVELLED geometry of {}: exit truncation {:.3f} deg, assembly shifted {:+.3f} mm along the deck z "
             "(= {:+.3f} mm along the Baseline z) relative to the run's design orbit (dz {:+.2f} -> {:+.2f} mm).\n"
             "System rotation R = {:+.4f} deg about z (re-solved on the levelled geometry; the run had {:+.4f}); quads {} deg on top.\n"
             "steps_baseline_mm/       one STEP per electrode + the combined assembly, MILLIMETRES (header and coordinates), shift included\n"
             "steps_baseline_mm/SI_Gap_Volume.stp     the NEGATIVE space between the two spiral electrodes (a solid, not an\n"
             "                         electrode; not in the combined assembly), same frame and units\n"
             "steps_baseline_mm/housing_outline.txt   the housing as a SHAPE for the central region: the polygon the wall is extruded\n"
             "                         from (before any extrusion) and the four corners of the exit opening, same frame and units\n"
             "exit_plane_spec.json     the exit plate's outer face at R (unshifted z; add the shift for the machine-frame z)\n"
             "efield_fine_best_baseline.pickle   levelled E-field (spiral x vscale + quads), V/m on a metre grid, shift included (use HCHC60_EFIELD_ZSHIFT_MM = 0)\n"
             "levelled_export.json     quads, vscale, shift, R of this export\n"
             "handoff_*_level.h5       openPMD hand-off at the electrode exit of the levelled full runs (vacuum / {cur_ma:g} mA), machine frame\n"
             "Levelling report: {}\n".format(os.path.basename(RUN), t_exit, 1e3 * shift, -1e3 * shift, 1e3 * dz0, 1e3 * dz, R, R_run,
                                             " / ".join("{:+.1f}".format(q) for q in a.quads), os.path.join(RUN, "sc_level.md"), cur_ma=float(lv.get("current_ma", 8.0))))
stamp("done: {} STEP files + {} hand-offs -> {}".format(len(files), len(copied), OUT))
