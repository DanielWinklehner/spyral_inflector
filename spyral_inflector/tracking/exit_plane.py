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


def outer_face(steps_dir, state, corridor=0.045, ahead=(-0.005, 0.070), align=0.90):
    """Fit the outer face of the housing's exit plate where the design orbit leaves.

    Housing triangles in a corridor around the beam ray past the exit whose normal is
    within acos(align) of the beam form the plate; they split into two parallel faces by
    their offset along the normal, and the far one is the outer face. Returns the beam
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
    n = (fa[:, None] * fn * np.sign(fn @ d)[:, None]).sum(0)
    n /= np.linalg.norm(n)
    s = fc @ n
    mid = 0.5 * (s.min() + s.max())
    outer = s > mid
    s_out = float(np.average(s[outer], weights=fa[outer]))
    s_in = float(np.average(s[~outer], weights=fa[~outer]))
    t = (s_out - p0 @ n) / (d @ n)
    p = p0 + t * d
    return {"point": p, "normal": n, "beam_dir": d, "thickness": s_out - s_in, "n_triangles": int(sel.sum()),
            "r": float(np.hypot(p[0], p[1])), "azimuth_deg": float(np.degrees(np.arctan2(p[1], p[0]))),
            "orbit_exit": p0, "orbit_exit_azimuth_deg": float(np.degrees(np.arctan2(p0[1], p0[0])))}


def rotation_for_gap(point, normal, gap_azimuth_deg=34.0, half_gap=5.0e-3):
    """Rotation about z [deg] that puts point + half_gap * normal on the radial plane at
    gap_azimuth_deg, i.e. the exit face half a gap before the gap along its own normal."""
    q = np.asarray(point, dtype=float) + half_gap * np.asarray(normal, dtype=float)
    r = gap_azimuth_deg - np.degrees(np.arctan2(q[1], q[0]))
    return float((r + 180.0) % 360.0 - 180.0)


def exit_plane_spec(steps_dir, state, gap_azimuth_deg=34.0, half_gap=5.0e-3, out_json=None, log=print):
    """Measure the outer face, solve the rotation and describe the exit plane at that
    rotation in the deck frame and in the Baseline (machine) frame; optionally write it."""
    f = outer_face(steps_dir, state)
    R = rotation_for_gap(f["point"], f["normal"], gap_azimuth_deg, half_gap)
    Rz = rotz(R)
    p, n, d = Rz @ f["point"], Rz @ f["normal"], Rz @ f["beam_dir"]
    q = p + half_gap * n
    mir = np.array([1.0, 1.0, -1.0])

    def frame(pp, nn, dd):
        return {"point_mm": (1e3 * pp).tolist(), "normal": nn.tolist(), "beam_dir": dd.tolist(),
                "r_mm": float(1e3 * np.hypot(pp[0], pp[1])), "azimuth_deg": float(np.degrees(np.arctan2(pp[1], pp[0]))),
                "normal_azimuth_deg": float(np.degrees(np.arctan2(nn[1], nn[0]))),
                "normal_tilt_deg": float(np.degrees(np.arcsin(nn[2])))}

    spec = {"rotation_deg": R, "gap_azimuth_deg": gap_azimuth_deg, "half_gap_mm": 1e3 * half_gap,
            "rule": "beam crossing of the plate's outer face + half_gap along the downstream normal lies on the "
                    "radial plane at gap_azimuth_deg",
            "pushed_point_azimuth_deg": float(np.degrees(np.arctan2(q[1], q[0]))),
            "plate_thickness_mm": 1e3 * f["thickness"], "n_triangles": f["n_triangles"],
            "unrotated": dict(frame(f["point"], f["normal"], f["beam_dir"]),
                              orbit_exit_mm=(1e3 * f["orbit_exit"]).tolist(),
                              orbit_exit_azimuth_deg=f["orbit_exit_azimuth_deg"]),
            "deck": frame(p, n, d),
            "baseline": dict(frame(p * mir, n * mir, d * mir),
                             note="machine frame: +z up, beam enters from +z; x, y, azimuth as in the deck frame")}
    log("exit plate outer face: r {:.2f} mm, azimuth {:+.3f} deg (orbit exit {:+.3f}), thickness {:.2f} mm, "
        "normal azimuth {:+.2f} deg".format(1e3 * f["r"], f["azimuth_deg"], f["orbit_exit_azimuth_deg"],
                                            1e3 * f["thickness"], spec["unrotated"]["normal_azimuth_deg"]))
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
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)
    return exit_plane_spec(a.steps, a.state, a.gap_azimuth, a.half_gap, a.out)


if __name__ == "__main__":
    main()
