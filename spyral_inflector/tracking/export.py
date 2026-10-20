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
    """Mirror every electrode file NNN_<Name>.stp/.step of a folder (optionally only names in `only`); written as .stp."""
    from .deck import electrode_cad_files, cad_name
    os.makedirs(out_dir, exist_ok=True)
    done = []
    for fn in electrode_cad_files(steps_dir):
        name = cad_name(fn)
        if only and name not in only:
            continue
        out_fn = os.path.splitext(fn)[0] + ".stp"
        mirror_step_z(os.path.join(steps_dir, fn), os.path.join(out_dir, out_fn))
        done.append(out_fn)
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
                         combined="assembly_baseline_mm.stp", shift_z=0.0, log=print):
    """Every NNN_<Name>.stp (or older .step) of steps_dir rotated about z by rotation_deg (each quadrupole by
    its own extra angle from quad_rotation, one entry per quad, on top), shifted by shift_z along z in the DECK frame (metres of the
    source geometry; a levelled assembly's axial shift), mirrored into the machine (Baseline)
    frame and written in MILLIMETRES with a matching unit declaration: one file per electrode
    plus one combined compound. The deck's own STEP exports carry metre-valued coordinates
    under a millimetre unit header (CAD imports them 1000x too small); source_unit="m" scales
    them by 1000 here so the files import at true size. Returns the written file names."""
    from OCC.Core.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
    from OCC.Core.IFSelect import IFSelect_RetDone
    from OCC.Core.Interface import Interface_Static
    from OCC.Core.gp import gp_Trsf, gp_Ax1, gp_Ax2, gp_Pnt, gp_Dir, gp_Vec
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCC.Core.BRep import BRep_Builder
    from OCC.Core.TopoDS import TopoDS_Compound
    from .deck import quad_index, electrode_cad_files, cad_name, CAD_EXTENSIONS, NEGATIVE_SUFFIX

    scale = {"m": 1000.0, "mm": 1.0}[source_unit]
    os.makedirs(out_dir, exist_ok=True)
    zax = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    builder = BRep_Builder()
    compound = TopoDS_Compound()
    builder.MakeCompound(compound)
    done = []
    for fn in electrode_cad_files(steps_dir):
        name = cad_name(fn)
        qi = quad_index(name)          # None for anything that is not a quadrupole pole
        qr = list(quad_rotation or ())
        ang = rotation_deg + (qr[qi] if qi is not None and qi < len(qr) else 0.0)
        rd = STEPControl_Reader()
        if rd.ReadFile(os.path.join(steps_dir, fn)) != IFSelect_RetDone:
            raise RuntimeError("cannot read " + fn)
        rd.TransferRoots()
        t_rot = gp_Trsf()
        t_rot.SetRotation(zax, np.radians(ang))
        t_shf = gp_Trsf()
        t_shf.SetTranslation(gp_Vec(0.0, 0.0, float(shift_z)))      # deck frame, before the mirror
        t_mir = gp_Trsf()
        t_mir.SetMirror(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)))
        t_scl = gp_Trsf()
        t_scl.SetScale(gp_Pnt(0, 0, 0), scale)
        shp = BRepBuilderAPI_Transform(rd.OneShape(), t_scl.Multiplied(t_mir.Multiplied(t_shf.Multiplied(t_rot))), True).Shape()
        builder.Add(compound, shp)
        Interface_Static.SetCVal("write.step.unit", "MM")
        w = STEPControl_Writer()
        w.Transfer(shp, STEPControl_AsIs)
        out_fn = os.path.splitext(fn)[0] + ".stp"
        if w.Write(os.path.join(out_dir, out_fn)) != IFSelect_RetDone:
            raise RuntimeError("cannot write " + out_fn)
        done.append(out_fn)
        log("  {:22s} rotated {:+8.3f} deg -> {}".format(name, ang, out_fn))
    # the negative space between the spiral electrodes (an extra CAD object, not an electrode): the same
    # transform as the spiral electrodes (R only, no quad angle), and kept OUT of the combined assembly
    for neg in sorted(f for f in os.listdir(steps_dir)
                      if os.path.splitext(f)[1].lower() in CAD_EXTENSIONS
                      and os.path.splitext(f)[0].lower().endswith(NEGATIVE_SUFFIX)):
        neg_out = os.path.splitext(neg)[0] + ".stp"
        rd = STEPControl_Reader()
        if rd.ReadFile(os.path.join(steps_dir, neg)) != IFSelect_RetDone:
            raise RuntimeError("cannot read " + neg)
        rd.TransferRoots()
        t_rot = gp_Trsf()
        t_rot.SetRotation(zax, np.radians(rotation_deg))
        t_shf = gp_Trsf()
        t_shf.SetTranslation(gp_Vec(0.0, 0.0, float(shift_z)))
        t_mir = gp_Trsf()
        t_mir.SetMirror(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)))
        t_scl = gp_Trsf()
        t_scl.SetScale(gp_Pnt(0, 0, 0), scale)
        shp = BRepBuilderAPI_Transform(rd.OneShape(), t_scl.Multiplied(t_mir.Multiplied(t_shf.Multiplied(t_rot))), True).Shape()
        Interface_Static.SetCVal("write.step.unit", "MM")
        w = STEPControl_Writer()
        w.Transfer(shp, STEPControl_AsIs)
        if w.Write(os.path.join(out_dir, neg_out)) != IFSelect_RetDone:
            raise RuntimeError("cannot write " + neg_out)
        done.append(neg_out)
        log("  {:22s} rotated {:+8.3f} deg -> {} (not an electrode)".format(os.path.splitext(neg)[0], rotation_deg, neg_out))
    if combined:
        Interface_Static.SetCVal("write.step.unit", "MM")
        w = STEPControl_Writer()
        w.Transfer(compound, STEPControl_AsIs)
        if w.Write(os.path.join(out_dir, combined)) != IFSelect_RetDone:
            raise RuntimeError("cannot write the combined assembly")
        done.append(combined)
    with open(os.path.join(out_dir, "FRAME.txt"), "w", encoding="utf-8") as fh:
        fh.write(MACHINE_FRAME_NOTE + "\nUnits: MILLIMETRES (coordinates and STEP header).\nSource: {}\n"
                 "Transform: rotation about z by {:+.4f} deg (quads {:+.2f} / {:+.2f} deg on top), shift {:+.3f} mm along "
                 "the deck z (= {:+.3f} mm along the Baseline z), mirror through z = 0, scale x{:g}.\nFiles: {}\n".format(
                     steps_dir, rotation_deg, quad_rotation[0], quad_rotation[1], 1e3 * shift_z, -1e3 * shift_z, scale,
                     ", ".join(done)))
    return done


HOUSING_OUTLINE_NAME = "housing_outline.txt"
# The negative space between the two spiral electrodes (geometry.SIGapVolume), exported next to the electrode
# files. All CAD exports are .stp (Daniel, 2026-10-16); the stem suffix "_negative" is what tells every loader
# (deck.is_electrode_cad) that this solid is NOT an electrode -- it fills the beam channel and must never enter a
# field solve or a collision check.
GAP_VOLUME_NAME = "SI_Gap_Volume_negative.stp"


def write_housing_outline(path, outer, aperture, curves=None, frame="deck", units="m",
                          comments=()):
    """The housing as a SHAPE: the outer outline and the four aperture corners, as text.

    The central region does not want the housing solid - it wants what the solid was made
    from, so it can build its own wall (no top and bottom plates, open radial ends) and cut
    its own opening. `outer` is the outline the wall is extruded from, (N, 2), taken before
    any extrusion; `aperture` is the four corners of the exit opening, (4, 3), going around
    the rectangle. Both in `frame` / `units`.

    `curves` says how to join the outer points: a list of (kind, i0, i1) with kind 'line' or
    'spline' and indices into `outer`, the last wrapping back to 0. The spiral inflector now
    simplifies the outline itself - one straight plate at the exit and one spline round the
    rest - so this carries the shape it actually built rather than a polygon the reader would
    have to simplify again, and differently. Without it the points are a closed polygon, which
    is what older files hold.

    Note the outline is NOT convex any more: the exit plate is a chord across the hull, so
    there is a re-entrant kink of a few degrees where it ends.

    Format: '#' comments, then 'key value' lines, then a section header 'name count'
    followed by that many rows. Deliberately trivial to parse without this module.
    """
    outer = np.asarray(outer, dtype=float).reshape(-1, 2)
    aperture = np.asarray(aperture, dtype=float).reshape(-1, 3)
    machine = os.environ.get("SI_MACHINE", "HCHC60")             # the file-name tag of the deck that runs
    machine = {"HCHC60": "HCHC-60", "HCHC15": "HCHC-1.5"}.get(machine, machine)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# " + machine + " spiral inflector -> central region: the housing outline.\n"
                 "# outer     the outline the housing wall is extruded from, BEFORE any\n"
                 "#           extrusion: one 'x y' per point, counter-clockwise, not closed.\n"
                 "# curves    how to join them: 'kind i0 i1', kind line|spline, indices into\n"
                 "#           outer, each curve starting where the last one ended and the last\n"
                 "#           wrapping back to 0. A spline interpolates every point from i0 to\n"
                 "#           i1. Absent means a plain closed polygon.\n"
                 "# aperture  the four corners of the housing's exit opening, 'x y z', around the\n"
                 "#           rectangle; they span the opening's mid-plane through the wall.\n")
        for line in comments:
            fh.write("# {}\n".format(line))
        fh.write("frame {}\nunits {}\n".format(frame, units))
        fh.write("outer {}\n".format(len(outer)))
        for x, y in outer:
            fh.write("  {:+.12g} {:+.12g}\n".format(x, y))
        if curves:
            fh.write("curves {}\n".format(len(curves)))
            for kind, i0, i1 in curves:
                fh.write("  {} {:d} {:d}\n".format(kind, int(i0), int(i1)))
        fh.write("aperture {}\n".format(len(aperture)))
        for x, y, z in aperture:
            fh.write("  {:+.12g} {:+.12g} {:+.12g}\n".format(x, y, z))
    return path


def read_housing_outline(path):
    """Inverse of write_housing_outline: {'frame', 'units', 'outer' (N,2), 'aperture' (4,3)}
    and, when the file has them, 'curves' as a list of (kind, i0, i1)."""
    out = {"frame": "deck", "units": "m"}
    rows, name, left = [], None, 0
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            if left > 0:
                if name == "curves":
                    kind, i0, i1 = line.split()
                    rows.append((kind, int(i0), int(i1)))
                else:
                    rows.append([float(v) for v in line.split()])
                left -= 1
                if left == 0:
                    out[name] = rows if name == "curves" else np.array(rows, dtype=float)
                continue
            key, _, rest = line.partition(" ")
            rest = rest.strip()
            if key in ("outer", "aperture", "curves"):
                name, left, rows = key, int(rest), []
                if left == 0:
                    out[key] = [] if key == "curves" else np.zeros(
                        (0, 2 if key == "outer" else 3))
            else:
                out[key] = rest
    for key in ("outer", "aperture"):
        if key not in out:
            raise ValueError("{}: no '{}' section".format(path, key))
    return out


def housing_outline_to_baseline_mm(src, dst, rotation_deg=0.0, shift_z=0.0, comments=()):
    """A deck-frame outline (metres) through the same transform steps_to_baseline_mm applies
    to the STEP solids: rotate about z by rotation_deg, shift by shift_z along the deck z,
    mirror through z = 0, scale to millimetres. The housing is not a quadrupole, so there is
    no extra per-quad angle here."""
    o = read_housing_outline(src)
    # carry the source's description of the exit opening (a CIRCLE's rim points look like any other polygon)
    with open(src, encoding="utf-8") as _fh:
        comments = tuple(comments) + tuple(ln[1:].strip() for ln in _fh if ln.startswith("# the exit opening"))
    if o.get("frame") != "deck" or o.get("units") != "m":
        raise ValueError("{}: expected a deck-frame outline in metres, got {} / {}".format(
            src, o.get("frame"), o.get("units")))
    c, s = np.cos(np.radians(rotation_deg)), np.sin(np.radians(rotation_deg))
    xy = lambda p: np.column_stack([1e3 * (c * p[:, 0] - s * p[:, 1]),                    # noqa: E731
                                    1e3 * (s * p[:, 0] + c * p[:, 1])])
    ap = np.asarray(o["aperture"], dtype=float)
    aperture = np.column_stack([xy(ap), -1e3 * (ap[:, 2] + float(shift_z))])
    # A rotation about z, a scaling and the z mirror all leave the point ORDER alone, so the
    # curve indices carry over untouched.
    return write_housing_outline(dst, xy(np.asarray(o["outer"], dtype=float)), aperture,
                                 curves=o.get("curves"), frame="baseline", units="mm",
                                 comments=tuple(comments) + (
                                     MACHINE_FRAME_NOTE,
                                     "from {} by rotation {:+.4f} deg about z, shift {:+.3f} mm along the deck z, "
                                     "mirror through z = 0, scale x1000.".format(src, rotation_deg, 1e3 * shift_z)))


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
