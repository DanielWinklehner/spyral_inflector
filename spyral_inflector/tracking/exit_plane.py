"""The housing's exit plate, and the rotation of the whole inflector that places it at the
first accelerating gap.

    python -m spyral_inflector.tracking.exit_plane --steps STEPS --state state.pickle --out spec.json

The rule (2026-09-10): the point where the design orbit crosses the OUTER face of the
exit plate, pushed half a gap along the plate's downstream normal, must lie on the radial
plane of the first gap. Since that plane contains the axis, the condition is simply that
the pushed point's azimuth equals the gap azimuth, which gives the rotation in closed form.
"""
import argparse
import json
import os

import numpy as np

from .deck import load_step_assembly, mesh_assembly, load_state


def rotz(deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def outer_face(steps_dir, state, corridor=0.045, ahead=(-0.005, 0.070), align=0.90,
               plane_tol=0.2, face_tol=3.0e-4):
    """Fit the outer face of the housing's exit plate where the design orbit leaves.

    Housing triangles in a corridor around the beam ray past the exit whose normal is
    within acos(align) of the beam are candidates; the largest group of them sharing a normal
    to within `plane_tol` degrees IS the plate, and that group splits into two parallel faces
    by their offset along the normal, the far one being the outer face. Returns the beam
    crossing point on it, the downstream unit normal, the beam direction, the plate
    thickness and the crossing's r / azimuth (deck frame, metres / degrees)."""
    st = load_state(state) if not isinstance(state, dict) else state
    trj, vd = np.asarray(st["trj_design"]), np.asarray(st["v_design"])
    p0, d = trj[-1], vd[-1] / np.linalg.norm(vd[-1])

    a = load_step_assembly(steps_dir)
    for u in [u for u, e in a.electrodes.items() if e.name != "Housing"]:
        a.electrodes.pop(u)
    if not a.electrodes:
        raise RuntimeError("no electrode named 'Housing' in {}".format(steps_dir))
    mesh_assembly(a)
    e = next(iter(a.electrodes.values()))
    V = np.asarray(e._gmsh_msh["vertices"], dtype=float)
    T = np.asarray(e._gmsh_msh["elements"], dtype=int)
    A, B, C = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    cen = (A + B + C) / 3.0
    nr = np.cross(B - A, C - A)
    area = 0.5 * np.linalg.norm(nr, axis=1)
    nr = nr / np.linalg.norm(nr, axis=1)[:, None]

    rel = cen - p0
    along = rel @ d
    perp = np.linalg.norm(rel - along[:, None] * d, axis=1)
    sel = (along > ahead[0]) & (along < ahead[1]) & (perp < corridor) & (np.abs(nr @ d) > align)
    if sel.sum() < 20:
        raise RuntimeError("only {} plate triangles found along the beam exit".format(int(sel.sum())))
    fc, fn, fa = cen[sel], nr[sel], area[sel]
    fs = np.sign(fn @ d)[:, None]

    # The plate is ONE plane, so take the largest group of triangles that share a normal, not the
    # average of everything pointing roughly the right way. Since the wall's exit plate became
    # exactly perpendicular to the beam, curved wall beside it passes the `align` test too, and
    # averaging it in swung the normal 0.975 deg and put the measured thickness at 3.94 mm
    # against a 2.00 mm wall. A degree matters here: the rotation is solved to make this normal
    # parallel to the gap, so a degree of error is a degree of non-parallel plate.
    # plane_tol 0.5 -> 0.2 on 2026-09-19: on the LEVELLED e016 wall (exit truncation 0.609 deg)
    # 0.5 still merged a few stray triangles and reported the thickness as 3.75 mm against the
    # 2.00 mm wall. The NORMAL was unaffected -- 174.833 vs 174.832 deg over tol 0.5 down to
    # 0.05, so that geometry's 1.26 deg plate drift was real -- but the thickness is the signal
    # that catches this class of bug at all, and it has to stay trustworthy.
    groups, order = [], np.argsort(-fa)
    for i in order:
        u = fn[i] * fs[i]
        for g in groups:
            if float(u @ g["n"]) > np.cos(np.radians(plane_tol)):
                g["i"].append(i)
                break
        else:
            groups.append({"n": u, "i": [i]})
    for g in groups:
        k = np.array(g["i"])
        g["area"] = float(fa[k].sum())
        g["n"] = (fa[k][:, None] * fn[k] * fs[k]).sum(0)
        g["n"] /= np.linalg.norm(g["n"])
    keep = np.array(max(groups, key=lambda g: g["area"])["i"])
    fc, fn, fa = fc[keep], fn[keep], fa[keep]
    n = (fa[:, None] * fn * np.sign(fn @ d)[:, None]).sum(0)
    n /= np.linalg.norm(n)
    s = fc @ n
    # The plate's two faces are the two OFFSETS that carry the most area (a flat face's triangles
    # all sit at one offset; face_tol separates faces). Splitting at the midpoint of the offset
    # range instead let a single stray parallel sliver decide the outer face: on the HCHC-1.5 t3
    # wall (2026-10-16) an 11 mm2 triangle 4.8 mm beyond the plate's 3,300 mm2 outer face became
    # "the outer face", read the 2.00 mm plate as 5.81 mm, moved the crossing 4.8 mm downstream
    # and turned the rotation fixed point 3.4 deg away from the clean build of the same shape.
    o = np.argsort(s)
    faces = sorted(np.split(o, np.where(np.diff(s[o]) > face_tol)[0] + 1), key=lambda k: -fa[k].sum())
    if len(faces) < 2:
        raise RuntimeError("the exit plate shows only one face along the beam exit")
    inner, outer = sorted(faces[:2], key=lambda k: s[k].mean())
    stray = float(1.0 - (fa[inner].sum() + fa[outer].sum()) / fa.sum())
    s_out = float(np.average(s[outer], weights=fa[outer]))
    s_in = float(np.average(s[inner], weights=fa[inner]))
    t = (s_out - p0 @ n) / (d @ n)
    p = p0 + t * d
    return {"point": p, "normal": n, "beam_dir": d, "thickness": s_out - s_in, "n_triangles": int(sel.sum()),
            "stray_area_fraction": stray, "n_offsets": len(faces),
            "r": float(np.hypot(p[0], p[1])), "azimuth_deg": float(np.degrees(np.arctan2(p[1], p[0]))),
            "orbit_exit": p0, "orbit_exit_azimuth_deg": float(np.degrees(np.arctan2(p0[1], p0[0])))}


def rotation_for_gap(point, normal, gap_azimuth_deg=34.0, half_gap=5.0e-3):
    """Rotation about z [deg] that puts point + half_gap * normal on the radial plane at
    gap_azimuth_deg, i.e. the exit face half a gap before the gap along its own normal."""
    q = np.asarray(point, dtype=float) + half_gap * np.asarray(normal, dtype=float)
    r = gap_azimuth_deg - np.degrees(np.arctan2(q[1], q[0]))
    return float((r + 180.0) % 360.0 - 180.0)


def rotation_parallel_to_gap(normal, gap_azimuth_deg=24.0):
    """Rotation about z [deg] that turns the exit plate PARALLEL to the gap's radial plane.

    That plane contains the axis and the radial direction at `gap_azimuth_deg`, so its normal
    lies at gap_azimuth_deg +- 90 deg; the plate is parallel to it when its own normal does
    too. The two branches differ by half a turn and both are parallel -- the smaller rotation
    is the one that keeps the beam crossing the gap the way it circulates, and it is also the
    one that stays near the rotation the offset rule would have given.

    This cannot also place the plate at the gap: one rotation, one condition. The plate ends up
    parallel to the gap at whatever stand-off the geometry gives (Daniel, 2026-09-19).
    """
    n = np.asarray(normal, dtype=float)
    az = np.degrees(np.arctan2(n[1], n[0]))
    cands = [float((gap_azimuth_deg + s * 90.0 - az + 180.0) % 360.0 - 180.0) for s in (+1, -1)]
    return min(cands, key=abs)


def plate_corners(steps_dir):
    """The two ends of the housing's exit plate, in the deck frame [m], as (outer, inner).

    The simplified wall is a spline plus the straight plate (see housing_outline.py), so the
    plate's two endpoints are its corners. `outer` is the one at the larger radius. Since the
    outline gained straight step edges at the plate ends (sub-mm lines between the plate and
    the spline), the plate is the LONGEST straight segment, not the first: on HCHC-1.5 t1nr
    the first line was the 0.8 mm step edge at r 37 mm, which put R 23 deg off.
    """
    from .export import HOUSING_OUTLINE_NAME, read_housing_outline

    ol = read_housing_outline(os.path.join(steps_dir, HOUSING_OUTLINE_NAME))
    pts = np.asarray(ol["outer"], dtype=float)
    lines = [c for c in (ol.get("curves") or []) if c[0] == "line"]
    if not lines:
        raise RuntimeError("{} has no straight segment: the wall is not the simplified "
                           "spline+plate outline".format(steps_dir))
    i0, i1 = max(((c[1], c[2]) for c in lines), key=lambda ij: np.hypot(*(pts[ij[1]] - pts[ij[0]])))
    a, b = pts[i0], pts[i1]
    return (a, b) if np.hypot(*a) >= np.hypot(*b) else (b, a)


def rotation_for_corner(corner_xy, gap_azimuth_deg=24.0, offset_mm=3.25, r_near=0.0):
    """Rotation about z [deg] putting a housing corner on a line PARALLEL to the gap's radial
    ray and `offset_mm` to the +x side of it.

    The ray from the origin at gap_azimuth_deg has unit vector u and left normal
    n = (-sin, cos); a point's signed perpendicular offset is p.n, positive on the -x side. The
    condition p(R).n = -offset_mm/1000 is one scalar equation in R, exactly like the other two
    rules, and it has a closed form: rotating by R takes p.n to |p| sin(phi + R - gap_azimuth),
    so phi + R - gap_azimuth = asin(target/|p|) or its supplement. Both branches are returned
    by the caller's choice of the one nearest r_near, because they place the housing on opposite
    sides of the machine.

    Unlike the parallel rule this says nothing about the plate's ANGLE -- it pins one corner and
    lets the plate lie where it falls.
    """
    p = np.asarray(corner_xy, dtype=float)[:2]
    rad = float(np.hypot(*p))
    target = -1e-3 * float(offset_mm)                 # +x side is the negative perpendicular side
    if rad < abs(target):
        raise ValueError("corner at r = {:.1f} mm can never be {:.2f} mm from the ray".format(
            1e3 * rad, offset_mm))
    phi = float(np.degrees(np.arctan2(p[1], p[0])))
    base = float(np.degrees(np.arcsin(target / rad)))
    cands = []
    for ang in (base, 180.0 - base):
        cands.append(float((gap_azimuth_deg + ang - phi + 180.0) % 360.0 - 180.0))
    return min(cands, key=lambda R: abs(R - r_near))


def exit_plane_spec(steps_dir, state, gap_azimuth_deg=34.0, half_gap=5.0e-3, out_json=None, log=print,
                    rule="offset", corner_offset_mm=3.25, corner_which="outer", r_near=0.0):
    """Measure the outer face, solve the rotation and describe the exit plane at that
    rotation in the deck frame and in the Baseline (machine) frame; optionally write it.

    rule "offset"   the plate's outer face, pushed half a gap along its normal, lands on the
                    gap's radial plane -- the plate sits at the gap but at an angle to it.
    rule "parallel" the plate is turned parallel to that plane instead, and sits wherever the
                    geometry puts it.
    rule "corner"   one CORNER of the housing -- an end of the exit plate, `corner_which` --
                    is put on a line parallel to the gap's ray, `corner_offset_mm` to the +x
                    side of it. Says nothing about the plate's angle. `r_near` picks the branch:
                    the two solutions put the housing on opposite sides of the machine.
    """
    f = outer_face(steps_dir, state)
    corner = None
    if rule == "parallel":
        R = rotation_parallel_to_gap(f["normal"], gap_azimuth_deg)
    elif rule == "offset":
        R = rotation_for_gap(f["point"], f["normal"], gap_azimuth_deg, half_gap)
    elif rule == "corner":
        outer_c, inner_c = plate_corners(steps_dir)
        corner = outer_c if corner_which == "outer" else inner_c
        R = rotation_for_corner(corner, gap_azimuth_deg, corner_offset_mm, r_near=r_near)
    else:
        raise ValueError("unknown exit-plane rule {!r}".format(rule))
    Rz = rotz(R)
    p, n, d = Rz @ f["point"], Rz @ f["normal"], Rz @ f["beam_dir"]
    q = p + half_gap * n
    mir = np.array([1.0, 1.0, -1.0])

    def frame(pp, nn, dd):
        return {"point_mm": (1e3 * pp).tolist(), "normal": nn.tolist(), "beam_dir": dd.tolist(),
                "r_mm": float(1e3 * np.hypot(pp[0], pp[1])), "azimuth_deg": float(np.degrees(np.arctan2(pp[1], pp[0]))),
                "normal_azimuth_deg": float(np.degrees(np.arctan2(nn[1], nn[0]))),
                "normal_tilt_deg": float(np.degrees(np.arcsin(nn[2])))}

    # how far the plate stands off the gap's radial plane, along the plate's own normal: the
    # signed distance of the face from that plane, positive when the plate is short of it
    n_gap = np.array([np.cos(np.radians(gap_azimuth_deg + 90.0)),
                      np.sin(np.radians(gap_azimuth_deg + 90.0)), 0.0])
    standoff = float(np.dot(n_gap, p))
    spec_angle = float((np.degrees(np.arctan2(n[1], n[0])) - (gap_azimuth_deg + 90.0) + 90.0) % 180.0 - 90.0)
    spec = {"rotation_deg": R, "gap_azimuth_deg": gap_azimuth_deg, "half_gap_mm": 1e3 * half_gap,
            "rule_name": rule,
            "rule": ("the plate is parallel to the gap's radial plane at gap_azimuth_deg; it is NOT placed at "
                     "the gap -- one rotation can only do one of the two" if rule == "parallel" else
                     "beam crossing of the plate's outer face + half_gap along the downstream normal lies on the "
                     "radial plane at gap_azimuth_deg"),
            "pushed_point_azimuth_deg": float(np.degrees(np.arctan2(q[1], q[0]))),
            "plate_to_gap_plane_mm": 1e3 * standoff,
            "plate_vs_gap_deg": spec_angle,
            "plate_thickness_mm": 1e3 * f["thickness"], "n_triangles": f["n_triangles"],
            "stray_area_fraction": f["stray_area_fraction"],
            "unrotated": dict(frame(f["point"], f["normal"], f["beam_dir"]),
                              orbit_exit_mm=(1e3 * f["orbit_exit"]).tolist(),
                              orbit_exit_azimuth_deg=f["orbit_exit_azimuth_deg"]),
            "deck": frame(p, n, d),
            "baseline": dict(frame(p * mir, n * mir, d * mir),
                             note="machine frame: +z up, beam enters from +z; x, y, azimuth as in the deck frame")}
    log("exit plate outer face: r {:.2f} mm, azimuth {:+.3f} deg (orbit exit {:+.3f}), thickness {:.2f} mm, "
        "normal azimuth {:+.2f} deg".format(1e3 * f["r"], f["azimuth_deg"], f["orbit_exit_azimuth_deg"],
                                            1e3 * f["thickness"], spec["unrotated"]["normal_azimuth_deg"]))
    if rule == "corner":
        cr = Rz @ np.array([corner[0], corner[1], 0.0])
        n_ray = np.array([-np.sin(np.radians(gap_azimuth_deg)), np.cos(np.radians(gap_azimuth_deg)), 0.0])
        got = 1e3 * float(cr @ n_ray)
        spec["corner"] = {"which": corner_which, "offset_mm_requested": corner_offset_mm,
                          "offset_mm_achieved": -got,
                          "unrotated_mm": (1e3 * np.asarray(corner[:2])).tolist(),
                          "rotated_mm": (1e3 * cr[:2]).tolist(),
                          "r_mm": float(1e3 * np.hypot(corner[0], corner[1])),
                          "azimuth_deg": float(np.degrees(np.arctan2(cr[1], cr[0])))}
        spec["rule"] = ("the {} corner of the exit plate lies on a line parallel to the gap ray at "
                        "gap_azimuth_deg, {:g} mm to the +x side of it; the plate's ANGLE is not "
                        "constrained".format(corner_which, corner_offset_mm))
        log("rotation for the {} CORNER on a line {:g} mm to +x of the {:g} deg ray: R = {:+.3f} deg "
            "-> corner at r {:.2f} mm, azimuth {:+.3f} deg, offset {:+.3f} mm (asked {:+.3f}); the plate "
            "itself lies {:+.3f} deg from parallel, {:+.1f} mm off the gap plane".format(
                corner_which, corner_offset_mm, gap_azimuth_deg, R, spec["corner"]["r_mm"],
                spec["corner"]["azimuth_deg"], spec["corner"]["offset_mm_achieved"], corner_offset_mm,
                spec_angle, 1e3 * standoff))
    elif rule == "parallel":
        log("rotation PARALLEL to the {:g} deg gap plane: R = {:+.3f} deg -> plate normal at {:+.3f} deg "
            "({:+.4f} deg from parallel), face at azimuth {:+.3f} deg, r {:.2f} mm, standing {:+.1f} mm "
            "off the gap plane".format(gap_azimuth_deg, R, spec["deck"]["normal_azimuth_deg"], spec_angle,
                                       spec["deck"]["azimuth_deg"], spec["deck"]["r_mm"], 1e3 * standoff))
    else:
        log("rotation for the {:g} deg gap, {:.1f} mm along the normal: R = {:+.3f} deg -> face at azimuth {:+.3f} deg, "
            "pushed point at {:+.3f} deg".format(gap_azimuth_deg, 1e3 * half_gap, R, spec["deck"]["azimuth_deg"],
                                                 spec["pushed_point_azimuth_deg"]))
    if out_json:
        os.makedirs(os.path.dirname(os.path.abspath(out_json)), exist_ok=True)
        with open(out_json, "w") as fh:
            json.dump(spec, fh, indent=2)
    return spec


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--steps", required=True)
    p.add_argument("--state", required=True)
    p.add_argument("--gap-azimuth", type=float, default=34.0)
    p.add_argument("--half-gap", type=float, default=5.0e-3, help="[m]")
    p.add_argument("--rule", choices=("offset", "parallel", "corner"), default="offset")
    p.add_argument("--corner-offset-mm", type=float, default=3.25)
    p.add_argument("--corner", choices=("outer", "inner"), default="outer")
    p.add_argument("--r-near", type=float, default=0.0, help="rule corner: the branch nearest this rotation [deg]")
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)
    return exit_plane_spec(a.steps, a.state, a.gap_azimuth, a.half_gap, a.out, rule=a.rule,
                           corner_offset_mm=a.corner_offset_mm, corner_which=a.corner, r_near=a.r_near)


if __name__ == "__main__":
    main()
