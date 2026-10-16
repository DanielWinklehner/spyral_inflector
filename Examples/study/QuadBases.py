"""Solve each quadrupole ONCE in isolation and place it by rotating its field.

The quad line's grounded plates are circular discs with circular holes -- surfaces of revolution
about z -- so a quadrupole turned inside them produces the same field turned by the same angle.
Solving quad i alone with those plates therefore gives a basis that can be rotated to any angle
instead of re-solved, which replaces both the per-rotation re-solve and the cos/sin skew identity
(the skew bases are no longer needed at all). The spiral basis is solved afterwards with the quad
POLES dropped, so nothing is solved twice.

What this trades away, and why it is measured rather than assumed: an isolated solve is missing
the grounded conductors of the rest of the machine (spiral electrodes, housing, the slotted
entrance aperture, the other quads' poles), which terminate the field. --validate differences the
result against the full-assembly bases of an existing run, on the beam corridor, against the
~1.4 % two honest re-solves of the same geometry differ by.

    python QuadBases.py --run E:\\HCHC60\\Results\\triplet\\par24b --validate
    python QuadBases.py --run <run> --out <basis dir> --res 0.0025 --h 0.004
"""
import argparse
import json
import os
import sys
import time

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
sys.path.insert(0, SCRIPTS)

from PyPATools.field import Field                                            # noqa: E402
from spyral_inflector.tracking.bem_reload import solve_step_assembly, DEFAULT_BOX   # noqa: E402
from spyral_inflector.tracking.deck import read_voltages, load_state          # noqa: E402
from spyral_inflector.tracking.frames import rotate_source_z                  # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run", required=True, help="a finished run: <run>/steps, <run>/geometry, <run>/exit_plane.json")
p.add_argument("--out", default=None, help="basis directory to write (default <run>/basis_iso)")
p.add_argument("--quads", type=float, nargs="+", default=[16.0, 24.0, 24.0], help="quad rotations on top of R [deg]")
p.add_argument("--rotation", type=float, default=None, help="R [deg] (default: the run's exit_plane.json)")
p.add_argument("--bfield", default=os.path.join(DECK, "Fields", "nocut_baseline_1mm.pickle"))
p.add_argument("--res", type=float, default=0.0025)
p.add_argument("--h", type=float, default=0.004)
p.add_argument("--box", type=float, nargs=6, default=list(DEFAULT_BOX))
p.add_argument("--energy-mev", type=float, default=None, help="default: the run's geometry summary")
p.add_argument("--unit", type=float, default=3500.0, help="pole voltage each basis is solved at")
p.add_argument("--no-spiral", action="store_true", help="quad bases only")
p.add_argument("--place", choices=("solve", "rotate"), default="solve",
               help="how a quad is put at its angle in the machine. solve (default): solved there "
                    "directly -- an isolated solve is only ~5 s, and measured 2026-09-19 it sits "
                    "0.4-0.55 %% of peak from the full-assembly field, under the 1.4 %% re-solve "
                    "floor. rotate: solved at zero and turned by R, which saves ~3 s per quad and "
                    "costs 1.6-1.9 %% of peak -- worth it only when many angles are needed")
p.add_argument("--validate", action="store_true",
               help="difference every rotated basis against <run>/basis, the full-assembly solve")
a = p.parse_args()

RUN = a.run
STEPS = os.path.join(RUN, "steps")
GEO = os.path.join(RUN, "geometry")
OUT = a.out or os.path.join(RUN, "basis_iso")
os.makedirs(OUT, exist_ok=True)
NQ = len(a.quads)


def stamp(msg):
    print("[{}] {}".format(time.strftime("%H:%M:%S"), msg), flush=True)


R = a.rotation
if R is None:
    with open(os.path.join(RUN, "exit_plane.json")) as fh:
        R = float(json.load(fh)["rotation_deg"])
ENERGY = a.energy_mev
if ENERGY is None:
    with open(os.path.join(GEO, "summary.json")) as fh:
        ENERGY = float(json.load(fh)["energy_mev"])

VOLTS = read_voltages(os.path.join(GEO, "voltages.csv"))
STATE = load_state(os.path.join(GEO, "state.pickle"))
common = dict(res=a.res, h=a.h, box=tuple(a.box), bfield=a.bfield, energy_mev=ENERGY, log=stamp)

stamp("R = {:+.3f} deg, quad rotations {} deg; placing each quad by {}".format(
    R, a.quads, "solving it there" if a.place == "solve" else "rotating a zero-angle solve by R"))
stamp("grid {:.1f} mm, mesh {:.1f} mm, box {}".format(1e3 * a.res, 1e3 * a.h, a.box))

# ---------------------------------------------------------------- the isolated quad bases
timings = {}
for i in range(NQ):
    tag = "q{}".format(i + 1)
    final = os.path.join(OUT, "ef_itp_{}.pickle".format(tag))
    if os.path.exists(final):
        stamp("{} already built, skipping".format(tag))
        continue
    t0 = time.time()
    qv = [0.0] * NQ
    qv[i] = a.unit
    rot = [0.0] * NQ
    rot[i] = a.quads[i]
    # The quad already carries its own angle a_i, and the sub-assembly sits at zero. Placing it in
    # the machine is then a RIGID rotation by R of poles and plates together, which is exact
    # whatever their shape -- unlike turning the quad alone inside the slotted entrance aperture.
    angle = R + a.quads[i]          # where the poles end up, either way
    if a.place == "solve":
        s = solve_step_assembly(STEPS, VOLTS, STATE, OUT, tag, spiral_voltage=0.0,
                                quad_voltages=qv, rotate=rot, rotate_all=R, isolate_quad=i,
                                test_particle=False, **common)
    else:
        s = solve_step_assembly(STEPS, VOLTS, STATE, OUT, tag + "_iso", spiral_voltage=0.0,
                                quad_voltages=qv, rotate=rot, rotate_all=0.0, isolate_quad=i,
                                test_particle=False, **common)
        f = Field.from_file(os.path.join(OUT, "ef_itp_{}_iso.pickle".format(tag)))
        stamp("  rotating {} rigidly into the machine by R = {:+.3f} deg".format(tag, R))
        rotate_source_z(f, R, label="{} isolated at {:+.1f} deg, rotated by R".format(
            tag, a.quads[i]), log=stamp).save(final)
    timings[tag] = {"wall_s": time.time() - t0, "n_triangles": s["n_triangles"],
                    "solve_s": s["solve_s"], "potential_s": s["potential_s"], "angle_deg": angle}
    stamp("  {} ({:.0f} s, {} triangles)".format(tag, timings[tag]["wall_s"], s["n_triangles"]))

# ---------------------------------------------------------------- the spiral basis, quads dropped
if not a.no_spiral and not os.path.exists(os.path.join(OUT, "ef_itp_spiral.pickle")):
    t0 = time.time()
    stamp("spiral basis at R = {:+.3f} deg with the quad POLES dropped (the plates stay: they are "
          "grounded conductors on the beam path)".format(R))
    s = solve_step_assembly(STEPS, VOLTS, STATE, OUT, "spiral", quad_voltages=[0.0] * NQ,
                            rotate=None, rotate_all=R, drop=["D*"], test_particle=True, **common)
    timings["spiral"] = {"wall_s": time.time() - t0, "n_triangles": s["n_triangles"],
                         "solve_s": s["solve_s"], "potential_s": s["potential_s"], "angle_deg": R}
    stamp("  spiral ({:.0f} s, {} triangles)".format(timings["spiral"]["wall_s"], s["n_triangles"]))

with open(os.path.join(OUT, "basis_meta.json"), "w") as fh:
    json.dump({"rotatable": True, "isolated_quads": True, "rotation_deg": R,
               "quad_rotations_deg": list(a.quads),
               "quad_angles_deg": [R + q for q in a.quads],
               "solved_with_own_angle": True, "rigid_rotation_deg": R, "place": a.place,
               "unit_V": a.unit, "res_m": a.res, "h_m": a.h, "box": list(a.box),
               "spiral_drops_quad_poles": not a.no_spiral, "timings": timings}, fh, indent=1)
stamp("wrote {}".format(os.path.join(OUT, "basis_meta.json")))
if timings:
    stamp("total {:.0f} s: ".format(sum(t["wall_s"] for t in timings.values())) +
          ", ".join("{} {:.0f} s ({} tri)".format(k, v["wall_s"], v["n_triangles"])
                    for k, v in timings.items()))

# ---------------------------------------------------------------- validation
if a.validate:
    ref_dir = os.path.join(RUN, "basis")
    print()
    print("=" * 100)
    print("ISOLATED + ROTATED vs FULL-ASSEMBLY SOLVE   (reference {})".format(ref_dir))
    print("=" * 100)
    print("{:>6}{:>12}{:>14}{:>14}{:>16}{:>16}".format(
        "basis", "peak |E|", "rms diff", "max diff", "rms, r<25mm", "of local field"))
    print("-" * 100)
    tags = ["q{}".format(i + 1) for i in range(NQ)] + ([] if a.no_spiral else ["spiral"])
    for tag in tags:
        ref_fn = os.path.join(ref_dir, "ef_itp_{}.pickle".format(tag))
        if not os.path.exists(ref_fn):
            print("{:>6}   no reference at {}".format(tag, ref_fn))
            continue
        ref = Field.from_file(ref_fn)
        new = Field.from_file(os.path.join(OUT, "ef_itp_{}.pickle".format(tag)))
        ra = [np.asarray(ref.grid_values[c], dtype=float) for c in "xyz"]
        na = [np.asarray(new.grid_values[c], dtype=float) for c in "xyz"]
        if ra[0].shape != na[0].shape:
            print("{:>6}   grid mismatch {} vs {}".format(tag, ra[0].shape, na[0].shape))
            continue
        m = np.sqrt(sum(v ** 2 for v in ra))
        diff = np.sqrt(sum((u - v) ** 2 for u, v in zip(ra, na)))
        peak = float(m.max())
        live = m > 0.01 * peak
        x = np.asarray(ref.grid["x"], dtype=float)
        y = np.asarray(ref.grid["y"], dtype=float)
        X, Y = np.meshgrid(x, y, indexing="ij")
        corr = (np.hypot(X, Y) < 0.025)[:, :, None] & np.ones(m.shape, dtype=bool)
        loc = corr & (m > 1e-3 * peak)
        print("{:>6}{:>12.4g}{:>13.3f}%{:>13.3f}%{:>15.3f}%{:>15.2f}%".format(
            tag, peak,
            100.0 * float(np.sqrt((diff[live] ** 2).mean())) / peak,
            100.0 * float(diff[live].max()) / peak,
            100.0 * float(np.sqrt((diff[corr] ** 2).mean())) / peak,
            100.0 * float(np.mean(diff[loc] / m[loc]))))
    print("-" * 100)
    print("for scale: two honest re-solves of the same geometry differ by ~1.4 % of peak, and "
          "dropping the quad line from the spiral basis moved it 0.28 % rms")
