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


def steps_to_baseline_mm(steps_dir, out_dir, rotation_deg=0.0, quad_rotation=(0.0, 0.0), source_unit="m",
                         combined="assembly_baseline_mm.step", log=print):
    """Every NNN_<Name>.step of steps_dir rotated about z by rotation_deg (quad 1 / quad 2 by
    their extra angles on top), mirrored into the machine (Baseline) frame and written in
    MILLIMETRES with a matching unit declaration: one file per electrode plus one combined
    compound. The deck's own STEP exports carry metre-valued coordinates under a
    millimetre unit header (CAD imports them 1000x too small); source_unit="m" scales
    them by 1000 here so the files import at true size. Returns the written file names."""
    from OCC.Core.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
    from OCC.Core.IFSelect import IFSelect_RetDone
    from OCC.Core.Interface import Interface_Static
    from OCC.Core.gp import gp_Trsf, gp_Ax1, gp_Ax2, gp_Pnt, gp_Dir
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCC.Core.BRep import BRep_Builder
    from OCC.Core.TopoDS import TopoDS_Compound
    from .deck import QUAD1, QUAD2

    scale = {"m": 1000.0, "mm": 1.0}[source_unit]
    os.makedirs(out_dir, exist_ok=True)
    zax = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    builder = BRep_Builder()
    compound = TopoDS_Compound()
    builder.MakeCompound(compound)
    done = []
    for fn in sorted(os.listdir(steps_dir)):
        if not fn.lower().endswith(".step") or "assembly" in fn.lower():
            continue
        name = os.path.splitext(fn)[0].split("_", 1)[1]
        ang = rotation_deg + (quad_rotation[0] if name in QUAD1 else quad_rotation[1] if name in QUAD2 else 0.0)
        rd = STEPControl_Reader()
        if rd.ReadFile(os.path.join(steps_dir, fn)) != IFSelect_RetDone:
            raise RuntimeError("cannot read " + fn)
        rd.TransferRoots()
        t_rot = gp_Trsf()
        t_rot.SetRotation(zax, np.radians(ang))
        t_mir = gp_Trsf()
        t_mir.SetMirror(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)))
        t_scl = gp_Trsf()
        t_scl.SetScale(gp_Pnt(0, 0, 0), scale)
        shp = BRepBuilderAPI_Transform(rd.OneShape(), t_scl.Multiplied(t_mir.Multiplied(t_rot)), True).Shape()
        builder.Add(compound, shp)
        Interface_Static.SetCVal("write.step.unit", "MM")
        w = STEPControl_Writer()
        w.Transfer(shp, STEPControl_AsIs)
        if w.Write(os.path.join(out_dir, fn)) != IFSelect_RetDone:
            raise RuntimeError("cannot write " + fn)
        done.append(fn)
        log("  {:22s} rotated {:+8.3f} deg -> {}".format(name, ang, fn))
    if combined:
        Interface_Static.SetCVal("write.step.unit", "MM")
        w = STEPControl_Writer()
        w.Transfer(compound, STEPControl_AsIs)
        if w.Write(os.path.join(out_dir, combined)) != IFSelect_RetDone:
            raise RuntimeError("cannot write the combined assembly")
        done.append(combined)
    with open(os.path.join(out_dir, "FRAME.txt"), "w", encoding="utf-8") as fh:
        fh.write(MACHINE_FRAME_NOTE + "\nUnits: MILLIMETRES (coordinates and STEP header).\nSource: {}\n"
                 "Transform: rotation about z by {:+.4f} deg (quads {:+.2f} / {:+.2f} deg on top), mirror through "
                 "z = 0, scale x{:g}.\nFiles: {}\n".format(steps_dir, rotation_deg, quad_rotation[0], quad_rotation[1],
                                                           scale, ", ".join(done)))
    return done


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
