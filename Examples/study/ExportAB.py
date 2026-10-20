"""Baseline-frame export of an A/B run folder (ABRect.py) -- the ExportLevelled.py outputs for a run that was NOT
levelled itself but tracked at its reference's levelled setting.

Takes the rotation, the levelling shift, the quads and the voltage scale from the A/B's own bunch json (what the
tracking used), the STEP set from <run>/steps, the bases from <run>/basis. Writes <run>/export_baseline/:
steps_baseline_mm/ (one .stp per electrode, the combined assembly, SI_Gap_Volume_negative.stp, housing_outline.txt),
exit_plane_spec.json (the corner rule re-measured on this geometry, for the record), efield_fine_best_baseline.pickle,
export.json, README.txt, and the hand-off if the run saved one.

    python ExportAB.py --run E:\\HCHC1.5\\Results\\triplet\\t1nr_c2_circ --quads 19.4926 21.6445 28.6825
"""
import argparse
import json
import os
import shutil
import time

p = argparse.ArgumentParser()
p.add_argument("--run", required=True)
p.add_argument("--tag", default="sc_all_ab", help="the bunch run whose setting is exported (bunch_<tag>.json)")
p.add_argument("--quads", type=float, nargs="+", required=True, help="TOTAL quad rotations on top of R [deg]")
p.add_argument("--rotation-rule", choices=("offset", "parallel", "corner"), default="corner")
p.add_argument("--gap-azimuth", type=float, default=24.0)
p.add_argument("--half-gap", type=float, default=0.00325)
p.add_argument("--corner", choices=("outer", "inner"), default="outer")
p.add_argument("--corner-offset-mm", type=float, default=3.25)
p.add_argument("--out", default=None)
a = p.parse_args()

import numpy as np  # noqa: E402

from spyral_inflector.tracking.exit_plane import exit_plane_spec  # noqa: E402
from spyral_inflector.tracking.export import (HOUSING_OUTLINE_NAME, steps_to_baseline_mm,  # noqa: E402
                                              field_to_machine_frame, housing_outline_to_baseline_mm)
from spyral_inflector.tracking.deck import superpose_basis  # noqa: E402
from PyPATools.field import Field  # noqa: E402

RUN = os.path.abspath(a.run)
NAME = os.path.basename(RUN)
OUT = a.out or os.path.join(RUN, "export_baseline")
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log_export.txt"), "a", encoding="utf-8")


def stamp(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


bj = json.load(open(os.path.join(RUN, "bunch_{}.json".format(a.tag))))
R, shift, vscale = float(bj["assembly_rotation_deg"]), float(bj.get("shift_z_m") or 0.0), float(bj.get("vscale") or 1.0)
sp = bj["superpose"]
# 'quads' is the full list of (voltage, skew) pairs; the flat q1/alpha1/q2/alpha2 keys stop at quad 2
quads_v = [float(x) for pair in sp["quads"] for x in pair]
steps = os.path.join(RUN, "steps")
state = os.path.join(RUN, "geometry", "state.pickle")
stamp("EXPORT {} at its tracked setting: R {:+.4f} deg, shift {:+.4f} mm, quads {} x vscale {:.4f}, quad turns {}".format(
    NAME, R, 1e3 * shift, quads_v, vscale, a.quads))

# 1. the plate rule on this geometry, for the record (the A/B tracked at its reference's R)
spec = exit_plane_spec(steps, state, a.gap_azimuth, a.half_gap, out_json=os.path.join(OUT, "exit_plane_spec.json"), log=stamp,
                       rule=a.rotation_rule, corner_offset_mm=a.corner_offset_mm, corner_which=a.corner, r_near=R)
stamp("  {} rule on this geometry: R {:+.4f} deg (tracked at {:+.4f}, difference {:+.4f})".format(
    a.rotation_rule, spec["rotation_deg"], R, spec["rotation_deg"] - R))

# 2. STEP files and the housing outline, placed as tracked
files = steps_to_baseline_mm(steps, os.path.join(OUT, "steps_baseline_mm"), rotation_deg=R, quad_rotation=tuple(a.quads),
                             combined="{}_inflector_{}_baseline_mm.stp".format(os.environ.get("SI_MACHINE", "HCHC60"), NAME),
                             shift_z=shift, log=stamp)
ol_src = os.path.join(steps, HOUSING_OUTLINE_NAME)
if os.path.exists(ol_src):
    housing_outline_to_baseline_mm(ol_src, os.path.join(OUT, "steps_baseline_mm", HOUSING_OUTLINE_NAME), rotation_deg=R,
                                   shift_z=shift, comments=("export of the A/B run {}".format(NAME),))

# 3. the E-field the tracking used: spiral x vscale + the quads, shifted, mirrored
e = superpose_basis(os.path.join(RUN, "basis"), *quads_v, vscale=vscale)
g = e.grid
e = Field.from_arrays(grid={"x": np.asarray(g["x"], dtype=float), "y": np.asarray(g["y"], dtype=float),
                            "z": np.asarray(g["z"], dtype=float) + shift},
                      values={c: np.asarray(e.grid_values[c]) for c in "xyz"}, dim=3, units="m",
                      label="E-field of {} shifted {:+.2f} mm".format(NAME, 1e3 * shift))
field_to_machine_frame(e, os.path.join(OUT, "efield_fine_best_baseline.pickle"), label="E-field of {} (machine frame)".format(NAME))
del e
copied = []
for fn in os.listdir(RUN):
    if fn.startswith("handoff_") and fn.endswith(".h5"):
        shutil.copyfile(os.path.join(RUN, fn), os.path.join(OUT, fn))
        copied.append(fn)
with open(os.path.join(OUT, "export.json"), "w") as fh:
    json.dump({"run": RUN, "tag": a.tag, "rotation_deg": R, "shift_z_m": shift, "quads": quads_v, "quad_rotation_deg": list(a.quads),
               "vscale": vscale, "efield": "efield_fine_best_baseline.pickle", "efield_zshift_mm": 0.0,
               "housing_outline": "steps_baseline_mm/" + HOUSING_OUTLINE_NAME, "rule_rotation_deg": spec["rotation_deg"]}, fh, indent=2)
with open(os.path.join(OUT, "README.txt"), "w", encoding="utf-8") as fh:
    fh.write("Baseline (machine) frame: right-handed, +z up, the beam enters from +z; azimuth counter-clockwise from above.\n"
             "Geometry of the A/B run {name} placed exactly as it was tracked: rotation R = {R:+.4f} deg about z, quads {q} deg on\n"
             "top, assembly shifted {s:+.3f} mm along the deck z (= {ms:+.3f} mm along the Baseline z).\n"
             "steps_baseline_mm/       one STEP per electrode + the combined assembly, MILLIMETRES (header and coordinates)\n"
             "steps_baseline_mm/SI_Gap_Volume_negative.stp   the NEGATIVE space between the two spiral electrodes (a solid, not\n"
             "                         an electrode; not in the combined assembly), same frame and units\n"
             "steps_baseline_mm/housing_outline.txt   the housing as a SHAPE: the wall polygon and the exit opening's points\n"
             "exit_plane_spec.json     the {rule} rule re-measured on this geometry ({rr:+.4f} deg)\n"
             "efield_fine_best_baseline.pickle   the E-field the tracking used (spiral x {v:.4f} + quads), shift included\n"
             "export.json              rotation, shift, quads, vscale of this export\n".format(
                 name=NAME, R=R, q=" / ".join("{:+.4f}".format(x) for x in a.quads), s=1e3 * shift, ms=-1e3 * shift,
                 rule=a.rotation_rule, rr=spec["rotation_deg"], v=vscale))
stamp("done: {} STEP files + {} hand-offs -> {}".format(len(files), len(copied), OUT))
