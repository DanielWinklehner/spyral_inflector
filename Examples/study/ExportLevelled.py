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
p.add_argument("--quads", type=float, nargs=2, default=[16.0, 24.0], help="quad rotations on top of R [deg]")
p.add_argument("--gap-azimuth", type=float, default=34.0)
p.add_argument("--half-gap", type=float, default=0.005)
p.add_argument("--out", default=None, help="default: <run>/export_baseline_levelled")
a = p.parse_args()

from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402
from spyral_inflector.tracking.export import steps_to_baseline_mm  # noqa: E402

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
    raise SystemExit("levelled geometry not found: {} / {}".format(steps, state))
with open(os.path.join(RUN, "exit_plane.json")) as fh:
    R_run = json.load(fh)["rotation_deg"]
stamp("EXPORT LEVELLED {}: exit truncation {:.3f} deg, dz {:+.2f} mm (run {:+.2f}, shift {:+.2f} mm), run R {:+.3f} deg".format(
    os.path.basename(RUN), t_exit, 1e3 * dz, 1e3 * dz0, 1e3 * shift, R_run))

# 1. the plate rule on the levelled STEP set (z-shift invariant; the truncation moves the crossing a little)
spec = exit_plane_spec(steps, state, a.gap_azimuth, a.half_gap, out_json=os.path.join(OUT, "exit_plane_spec.json"), log=stamp)
R = spec["rotation_deg"]
stamp("  R re-solved on the levelled geometry: {:+.3f} deg (run {:+.3f}, change {:+.3f})".format(R, R_run, R - R_run))
if abs(R - R_run) > 0.2:
    stamp("  NOTE: |dR| > 0.2 deg -- the levelled geometry was rebuilt in the field rotated by the run's R; a re-optimization "
          "at the new R (RunFinalPush phase 5 loop) would be self-consistent")

# 2. STEP export: rotate by R (quads on top), shift by the levelled dz in the deck frame, mirror, millimetres
files = steps_to_baseline_mm(steps, os.path.join(OUT, "steps_baseline_mm"), rotation_deg=R, quad_rotation=tuple(a.quads),
                             combined="HCHC60_inflector_{}_levelled_baseline_mm.step".format(os.path.basename(RUN)),
                             shift_z=shift, log=stamp)

# 3. hand-offs of the levelled full runs (tracked with the shift, machine frame at the electrode exit)
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
             "System rotation R = {:+.4f} deg about z (re-solved on the levelled geometry; the run had {:+.4f}); quads {:+.1f} / {:+.1f} deg on top.\n"
             "steps_baseline_mm/       one STEP per electrode + the combined assembly, MILLIMETRES (header and coordinates), shift included\n"
             "exit_plane_spec.json     the exit plate's outer face at R (unshifted z; add the shift for the machine-frame z)\n"
             "handoff_*_level.h5       openPMD hand-off at the electrode exit of the levelled full runs (vacuum / 8 mA), machine frame\n"
             "Levelling report: {}\n".format(os.path.basename(RUN), t_exit, 1e3 * shift, -1e3 * shift, 1e3 * dz0, 1e3 * dz, R, R_run,
                                             a.quads[0], a.quads[1], os.path.join(RUN, "sc_level.md")))
stamp("done: {} STEP files + {} hand-offs -> {}".format(len(files), len(copied), OUT))
