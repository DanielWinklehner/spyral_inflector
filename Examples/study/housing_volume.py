"""Volume of a run's housing solid (002_Housing.stp), optionally checked against an expected value -- the one test that
catches an exit opening the boolean silently failed to cut (2026-10-20: a one-sided disk tool missed the plate and the
beam ran into a solid wall; the analytic reach check had passed).

    python housing_volume.py <steps_dir> [--expect-cm3 41.309 --tol-cm3 0.02]
"""
import argparse
import os

import gmsh

p = argparse.ArgumentParser()
p.add_argument("steps")
p.add_argument("--expect-cm3", type=float, default=None)
p.add_argument("--tol-cm3", type=float, default=0.02)
a = p.parse_args()
fn = next(os.path.join(a.steps, f) for f in sorted(os.listdir(a.steps)) if f.endswith(("_Housing.stp", "_Housing.step")))
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
ents = gmsh.model.occ.importShapes(fn)
gmsh.model.occ.synchronize()
vol = 1e6 * sum(gmsh.model.occ.getMass(3, t) for d, t in ents if d == 3)   # deck STEPs carry metres
gmsh.finalize()
print("housing {}: {:.3f} cm3{}".format(fn, vol, "" if a.expect_cm3 is None else " (expected {:.3f} +- {:.3f})".format(a.expect_cm3, a.tol_cm3)))
if a.expect_cm3 is not None and abs(vol - a.expect_cm3) > a.tol_cm3:
    raise SystemExit("housing volume {:.3f} cm3 is not the expected {:.3f}: the opening is not what was asked for".format(vol, a.expect_cm3))
