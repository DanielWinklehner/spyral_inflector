"""Exports for the central-region model in the MACHINE frame.

The deck simulates the mirror image of the machine through the median plane (the beam is
tracked upward in +z through the magnet's lower half; see Docs/coordinate_convention.md):
x, y, azimuths and the sense of rotation agree with the machine, only z is reversed. So
every export applies z -> -z: field maps (E_z changes sign, E_x/E_y do not; the grid is
re-ordered so z still ascends), STEP solids (a mirror through the median plane) and
particle states (z, v_z / p_z negated).

    python -m spyral_inflector.tracking.export --field ef_itp_final.pickle --steps STEPS_DIR --out OUT_DIR
"""
import argparse
import os

import numpy as np

from PyPATools.field import Field

MACHINE_FRAME_NOTE = ("machine frame: origin at the machine centre in the median plane, +z up (towards the RFQ), "
                      "+x = azimuth 0 on a magnet hill, azimuth counter-clockwise from above; the beam comes down "
                      "the axis (-z) and circulates counter-clockwise; B_z < 0 (down). Mirror image (z -> -z) of "
                      "the spyral_inflector deck frame.")


def field_to_machine_frame(field_in, field_out=None, label=None):
    """Mirror a deck-frame Field pickle (or Field) through the median plane and save it."""
    f = Field.from_file(field_in) if isinstance(field_in, str) else field_in
    g, v = f.grid, f.grid_values
    z = np.asarray(g["z"], dtype=float)
    order = np.argsort(-z)                                    # -z ascending
    grid = {"x": np.asarray(g["x"], dtype=float), "y": np.asarray(g["y"], dtype=float), "z": -z[order]}
    values = {"x": np.asarray(v["x"])[:, :, order], "y": np.asarray(v["y"])[:, :, order], "z": -np.asarray(v["z"])[:, :, order]}
    out = Field.from_arrays(grid=grid, values=values, dim=3, units="m", label=label or "{} ({})".format(getattr(f, "label", "E-field"), "machine frame"))
    if field_out:
        out.save(field_out)
        with open(os.path.splitext(field_out)[0] + "_FRAME.txt", "w", encoding="utf-8") as fh:
            fh.write(MACHINE_FRAME_NOTE + "\nSource: {}\nTransform: z -> -z, E_z -> -E_z, E_x/E_y unchanged, grid re-ordered.\n".format(field_in))
    return out


def mirror_step_z(step_in, step_out):
    """Write the mirror image (through the median plane z = 0) of a STEP solid with OCC."""
    from OCC.Core.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
    from OCC.Core.IFSelect import IFSelect_RetDone
    from OCC.Core.gp import gp_Trsf, gp_Ax2, gp_Pnt, gp_Dir
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    reader = STEPControl_Reader()
    if reader.ReadFile(step_in) != IFSelect_RetDone:
        raise RuntimeError("could not read " + step_in)
    reader.TransferRoots()
    shape = reader.OneShape()
    trsf = gp_Trsf()
    trsf.SetMirror(gp_Ax2(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0)))     # mirror through the plane z = 0
    mirrored = BRepBuilderAPI_Transform(shape, trsf, True).Shape()
    writer = STEPControl_Writer()
    writer.Transfer(mirrored, STEPControl_AsIs)
    if writer.Write(step_out) != IFSelect_RetDone:
        raise RuntimeError("could not write " + step_out)
    return step_out


def steps_to_machine_frame(steps_dir, out_dir, only=None):
    """Mirror every NNN_<Name>.step of a folder (optionally only names in `only`)."""
    os.makedirs(out_dir, exist_ok=True)
    done = []
    for fn in sorted(os.listdir(steps_dir)):
        if not fn.lower().endswith(".step") or "assembly" in fn.lower():
            continue
        name = os.path.splitext(fn)[0].split("_", 1)[1]
        if only and name not in only:
            continue
        mirror_step_z(os.path.join(steps_dir, fn), os.path.join(out_dir, fn))
        done.append(fn)
    with open(os.path.join(out_dir, "FRAME.txt"), "w", encoding="utf-8") as fh:
        fh.write(MACHINE_FRAME_NOTE + "\nSource: {}\nTransform: mirror through z = 0 (the spiral's handedness flips with it).\nFiles: {}\n".format(steps_dir, ", ".join(done)))
    return done


def state_to_machine_frame(r, v=None):
    """Mirror particle positions (N,3) and velocities/momenta (N,3): z and the z-component negated."""
    r2 = np.array(r, dtype=float)
    r2[:, 2] *= -1.0
    if v is None:
        return r2
    v2 = np.array(v, dtype=float)
    v2[:, 2] *= -1.0
    return r2, v2


def main(argv=None):
    p = argparse.ArgumentParser(description="machine-frame exports (z mirrored) for the central region")
    p.add_argument("--field", default=None, help="deck-frame E-field pickle (ef_itp_<tag>.pickle)")
    p.add_argument("--steps", default=None, help="deck-frame STEP folder")
    p.add_argument("--only", nargs="*", default=None, help="electrode names to export from the STEP folder (default: all)")
    p.add_argument("--out", required=True, help="output folder")
    a = p.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    if a.field:
        out = os.path.join(a.out, os.path.splitext(os.path.basename(a.field))[0] + "_machine.pickle")
        field_to_machine_frame(a.field, out)
        print("wrote", out)
    if a.steps:
        done = steps_to_machine_frame(a.steps, os.path.join(a.out, "steps_machine"), only=a.only)
        print("mirrored {} STEP files into {}".format(len(done), os.path.join(a.out, "steps_machine")))


if __name__ == "__main__":
    main()
