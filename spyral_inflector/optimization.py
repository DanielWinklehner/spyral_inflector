# from .global_variables import *
from .vector import Vector
import numpy as np
import time


def optimize_fringe(si, initial_guess=(None, None), maxiter=10, tol=1e-1, res=0.002,
                    apply_shift=False, robust_exit=False):
    """
    This function optimizes the length of the spiral inflector to adjust for fringe fields.
    :param si: the spiral inflector object
    :param initial_guess: tuple, list or ndarray containing angle adjustment for entrance and exit
                          units: degrees
    :param maxiter: Maximum number of iterations to find the optimum entrance and exit adjustment before giving up
    :param tol: maximum deviation tolerated for a successful solution (degrees)
    :param res:
    :return:
    """

    analytic_params = si.analytic_parameters
    analytic_vars = si.analytic_variables
    track_vars = si.track_variables

    assert si.method != "analytical", "You can't optimize using an analytical model!"

    print("Starting the optimization process...")

    # Start with an initial guess for bmin, bmax. If None, use 0.0
    # TODO: Update to use current internal adjustment, if present -DW
    _db_init = np.array([_db if _db is not None else 0.0 for _db in initial_guess])

    # This doesn't work -DW
    # _db_init[np.where(_db_init is None)] = 0.0

    _db = _db_init[:]

    # Half the length of the cube for electric field calculation
    _hcl = (analytic_params["gap"] + 2.0 * analytic_params["dx"])
    # _ion = si._params_analytic["ion"]  # type: IonSpecies
    x_axis = Vector([1.0, 0.0, 0.0])
    y_axis = Vector([0.0, 1.0, 0.0])
    z_axis = Vector([0.0, 0.0, 1.0])

    _r = np.zeros([1, 3])

    # --- Optimize entrance fringe --- #
    deviation_x = 1e20  # Some large number as initialization
    it = 0

    entrance_opt = []  # Store the information from the optimization

    if track_vars["shift"] is None:
        track_vars["shift"] = np.zeros(3)

    while abs(deviation_x) > tol:

        # Apply the angle correction (first time use initial guess)
        si.set_blim(b_min=0.0 + _db[0])
        # (Re-)calculate the new geometry and BEM++ solution
        si.generate_meshed_model()
        si.solve()

        # Trajectory starting point:
        _trj = analytic_vars["trj_design"]  # type: np.ndarray
        _v_des = analytic_vars["v_design"]  # type: np.ndarray

        # Calculate E-Field upwards of starting point
        xs, ys, zs = _trj[0]

        if si.solver == 'bempp':
            si.calculate_potential(limits=((xs - 0.5 * _hcl, xs + 0.5 * _hcl),
                                           (ys - 0.5 * _hcl, ys + 0.5 * _hcl),
                                           (zs - 1.9 * _hcl, zs + 0.1 * _hcl)),
                                   res=res,
                                   domain_decomp=(4, 4, 4),
                                   overlap=0)

        si.calculate_efield()

        _r, _v = si.track(r_start=_trj[0],  # Use starting point of design particle
                          v_start=-_v_des[0],  # Reverse direction of design particle
                          nsteps=5000,
                          dt=1e-11,
                          omit_b=True)

        trj_dir = Vector(_r[-1] - _r[-100])
        deviation_x = 90.0 - np.rad2deg(trj_dir.angle_with(x_axis))
        deviation_y = 90.0 - np.rad2deg(trj_dir.angle_with(y_axis))

        print("Current entrance adjustment: {:.4f}, "
              "deviation from z-axis in x-dir: {:.4f} degrees, "
              "deviation from z-axis in y-dir: {:.4f} degrees".format(_db[0], deviation_x, deviation_y))

        # Assignment, not "-=": this is the lateral offset of the back-tracked
        # particle for the CURRENT geometry. Accumulating it summed one measurement
        # per iteration, so with maxiter=12 the stored shift was ~12x too large.
        track_vars["shift"][0] = -_r[-1][0]
        track_vars["shift"][1] = -_r[-1][1]

        entrance_opt.append([it, _db[0], deviation_x, deviation_y])

        it += 1

        if it == maxiter:
            print("Entrance Fringe: Maximum number of iterations has been reached. Breaking.")
            break

        if it == 1:
            _db[0] += deviation_x
        else:
            _db[0] += 0.5 * deviation_x  # Dampen the oscillations a bit

    entrance_opt = np.array(entrance_opt)
    _db[0] = entrance_opt[np.abs(entrance_opt[:, 2]).argsort()][0, 1]

    # TODO: This is a work in progress
    # # Send the devation from the y-axis for the entrance to the optimziation parameters
    # si._variables_optimization["x_rot"] = -entrance_opt[np.abs(entrance_opt[:, 2]).argsort()][0, 3] * (np.pi / 180.0)
    # if si._params_exp["y_opt"]:
    #     print("Using x rotation angle of {:.4f} deg.".format(np.rad2deg(si._variables_optimization["x_rot"])))

    # Construct the rotation matrix
    # rot = np.array([[1.0, 0.0, 0.0],
    #                 [0.0, np.cos(si._variables_optimization["x_rot"]), -np.sin(si._variables_optimization["x_rot"])],
    #                 [0.0, np.sin(si._variables_optimization["x_rot"]), np.cos(si._variables_optimization["x_rot"])]])

    # --- Optimize exit fringe --- #
    deviation = 1e20  # Some large number as initialization
    it = 0

    exit_opt = []  # Store the information from the optimization
    _hist = []  # (b_max adjustment, deviation) samples, for the robust_exit solver

    while abs(deviation) > tol:

        # Apply the angle correction (first time use initial guess)
        si.set_blim(b_max=90.0 + _db[1])
        # (Re-)calculate the new geometry and BEM++ solution
        si.generate_meshed_model()
        si.solve()

        # Trajectory starting point:
        _trj = analytic_vars["trj_design"]  # type: np.ndarray
        _v_des = analytic_vars["v_design"]  # type: np.ndarray

        # Rotate the initial trajectory point and velocity
        # _trj = np.matmul(rot, _trj.T).T
        # _v_des = np.matmul(rot, _v_des.T).T

        _ns = analytic_params["ns"]  # type: int
        start_idx = int(0.9 * _ns)  # start the tracking "10 percent" into the spiral inflector exit
        xs, ys, zs = rs = _trj[start_idx]
        vs = _v_des[start_idx]

        if si.solver == 'bempp':
            # Calculate E-Field
            # TODO: Better way to determine the fringe field region
            si.calculate_potential(limits=((xs - 2.0 * _hcl, xs + 2.0 * _hcl),
                                           (ys - 2.0 * _hcl, ys + 2.0 * _hcl),
                                           (zs - _hcl, zs + _hcl)),
                                   res=res,
                                   domain_decomp=(4, 4, 4),
                                   overlap=0)

        si.calculate_efield()

        _r, _v = si.track(r_start=rs,  # Use point close to exit along design particle
                          v_start=vs,  # Regular direction of design particle
                          nsteps=5000,
                          dt=1e-11,
                          omit_b=False)

        trj_dir = Vector(_r[-1] - _r[-100])
        deviation = 90.0 - np.rad2deg(trj_dir.angle_with(z_axis))

        print("Current exit adjustment: {:.4f}, "
              "deviation from xy-plane: {:.4f} degrees".format(_db[1], deviation))

        # See the entrance loop: assignment rather than accumulation.
        track_vars["shift"][2] = -_r[-1, 2]

        exit_opt.append([it, _db[1], deviation])

        it += 1

        if it == maxiter:
            print("Exit Fringe: Maximum number of iterations has been reached. Breaking.")
            break

        if robust_exit:
            # Secant on deviation(b_max), with bisection once zero is bracketed.
            # The damped fixed-point step below overshoots badly: the response is
            # strongly nonlinear, and a -2.15 deg step was measured to take the
            # deviation from -3.31 to +0.33 deg.
            _hist.append((_db[1], deviation))

            _bracket = [(x, f) for x, f in _hist if f > 0.0], [(x, f) for x, f in _hist if f < 0.0]

            if _bracket[0] and _bracket[1]:
                # Bracketed: bisect between the closest straddling pair.
                _pos = min(_bracket[0], key=lambda p: abs(p[1]))
                _neg = min(_bracket[1], key=lambda p: abs(p[1]))
                _db[1] = 0.5 * (_pos[0] + _neg[0])

            elif len(_hist) >= 2 and _hist[-1][1] != _hist[-2][1]:
                # Secant step through the last two samples.
                (_x0, _f0), (_x1, _f1) = _hist[-2], _hist[-1]
                _db[1] = _x1 - _f1 * (_x1 - _x0) / (_f1 - _f0)

            else:
                _db[1] += deviation

        elif it == 1:
            _db[1] += deviation
        else:
            _db[1] += 0.65 * deviation  # Dampen the oscillations a bit

    exit_opt = np.array(exit_opt)
    _db[1] = exit_opt[np.abs(exit_opt[:, 2]).argsort()][0, 1]

    print("Applied entrance adjustment: {:.4f}, and exit adjustment: {:.4f}.".format(_db[0], _db[1]))

    if si.debug:
        print("Entrance Optimization:")
        print(entrance_opt)
        print("Exit Optimization:")
        print(exit_opt)

    if apply_shift:
        # The entrance loop measured the lateral offset of a back-tracked particle
        # and the exit loop the vertical one. Store it where generate_solid_assembly
        # picks it up, so the electrodes are translated to centre the incoming beam
        # on the machine axis; doing it here rather than on the assembly directly
        # means it survives the re-meshing below and reaches the exported geometry.
        # Lateral AND axial. The inflector assembly is aligned to the beam in all
        # three axes; the axial component brings the exit back onto the median plane.
        _applied = np.array(track_vars["shift"], dtype=float)

        si.track_variables["shift_applied"] = _applied

        # The quadrupoles and their apertures stay centred on the incoming beam
        # laterally, but follow the inflector axially so the drift between quad 2 and
        # the inflector entrance -- 38 mm, where the matching optics live -- is
        # preserved. The source-to-quad-1 gap is only ~5 mm, over which the beam
        # barely evolves, so changing that instead is optically irrelevant.
        si.track_variables["shift_applied_lab"] = np.array([0.0, 0.0, _applied[2]])

        print("Applying shift: inflector dx={:.4f} dy={:.4f} dz={:.4f} mm, "
              "quadrupoles dz={:.4f} mm".format(
                  1000.0 * _applied[0], 1000.0 * _applied[1], 1000.0 * _applied[2],
                  1000.0 * _applied[2]))

    # Both loops pick their best-scoring adjustment above, but nothing ever put it
    # back: the limits still in effect were whatever the final iteration happened
    # to set, which is not the same thing once a loop overshoots and comes back.
    si.set_blim(b_min=0.0 + _db[0], b_max=90.0 + _db[1])

    # Recalculate the new geometry and BEM++ solution one last time
    si.initialize()
    si.generate_meshed_model()
    si.solve()

    print("Done optimizing!")

    si.track_variables = track_vars

    return track_vars["shift"]


# --------------------------------------------------------------------------------------- #
# Full-trajectory fringe optimization
# --------------------------------------------------------------------------------------- #

def _safe_ratio(num, den, mask):
    """num / den where mask is set and den is non-zero, else 0. Degenerate
    (zero-area) triangles make some of the edge denominators vanish."""
    ok = mask & (den != 0.0)
    return np.where(ok, num / np.where(ok, den, 1.0), 0.0)


def _point_triangle_distance(p, a, b, c):
    """Exact distance from each point p[i] to the triangle (a[i], b[i], c[i]).

    Vectorised form of the closest-point-on-triangle test from Ericson, Real-Time
    Collision Detection, ch. 5.1.5. The seven Voronoi regions are applied as masks
    in reverse order of the original early returns, so the first test in the
    original (vertex A) is assigned last and wins where regions touch.
    """
    ab, ac, ap = b - a, c - a, p - a
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)

    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2

    with np.errstate(divide="ignore", invalid="ignore"):
        denom = va + vb + vc
        ok = denom != 0.0
        v = np.where(ok, vb / np.where(ok, denom, 1.0), 0.0)
        w = np.where(ok, vc / np.where(ok, denom, 1.0), 0.0)
        q = a + ab * v[:, None] + ac * w[:, None]  # face interior

        m = (va <= 0.0) & (d4 - d3 >= 0.0) & (d5 - d6 >= 0.0)  # edge BC
        t = _safe_ratio(d4 - d3, (d4 - d3) + (d5 - d6), m)
        q[m] = b[m] + (c[m] - b[m]) * t[m][:, None]

        m = (vb <= 0.0) & (d2 >= 0.0) & (d6 <= 0.0)  # edge AC
        t = _safe_ratio(d2, d2 - d6, m)
        q[m] = a[m] + ac[m] * t[m][:, None]

        m = (d6 >= 0.0) & (d5 <= d6)  # vertex C
        q[m] = c[m]

        m = (vc <= 0.0) & (d1 >= 0.0) & (d3 <= 0.0)  # edge AB
        t = _safe_ratio(d1, d1 - d3, m)
        q[m] = a[m] + ab[m] * t[m][:, None]

        m = (d3 >= 0.0) & (d4 <= d3)  # vertex B
        q[m] = b[m]

        m = (d1 <= 0.0) & (d2 <= 0.0)  # vertex A
        q[m] = a[m]

    return np.linalg.norm(p - q, axis=1)


def _surface_distance(points, verts, tris, k=32):
    """Distance from each point to a triangulated surface.

    The k triangles with the nearest centroids are found with a KD-tree and the exact
    point-triangle distance is taken over those. At the mesh sizes used here (h of a
    few mm, distances of ~10 mm) the closest triangle is always among them.
    """
    from scipy.spatial import cKDTree

    points = np.asarray(points, dtype=float)
    tri_pts = verts[tris]  # (M, 3, 3)
    k = int(min(k, len(tris)))
    _, idx = cKDTree(tri_pts.mean(axis=1)).query(points, k=k)
    idx = np.asarray(idx).reshape(len(points), k)

    cand = tri_pts[idx.ravel()]
    p = np.repeat(points, k, axis=0)
    d = _point_triangle_distance(p, cand[:, 0], cand[:, 1], cand[:, 2])

    return d.reshape(len(points), k).min(axis=1)


def _electrode_surface(si, name_fragment):
    """Vertices (N, 3) and triangles (M, 3) of one electrode, from the assembled
    BEM mesh -- i.e. the surfaces the field was actually solved on, translation
    included."""
    numerical_vars = si.numerical_variables
    mesh = numerical_vars["full mesh"]

    domain = None
    for _electrode in numerical_vars["objects"].electrodes.values():
        if name_fragment.lower() in _electrode.name.lower():
            domain = _electrode.bempp_domain
            break

    if domain is None:
        raise RuntimeError("No electrode with '{}' in its name in the assembly".format(name_fragment))

    verts = np.asarray(mesh["verts"], dtype=float).T
    elems = np.asarray(mesh["elems"], dtype=int).T
    domns = np.asarray(mesh["domns"], dtype=int)

    return verts, elems[domns == domain]


def _optical_frame(v, b, kp):
    """Width (h) and gap (u) unit vectors along an orbit, as generate_numerical_geometry
    builds them: the optical frame of the velocity, rotated about it by the tilt
    angle atan(kp sin b). +u points towards the upper electrode (geo[0..4])."""
    vx, vy, vz = v[:, 0], v[:, 1], v[:, 2]
    v2 = np.hypot(vx, vy)
    v3 = np.linalg.norm(v, axis=1)

    h = np.tile([0.0, -1.0, 0.0], (len(v), 1))
    u = np.tile([-1.0, 0.0, 0.0], (len(v), 1))
    ok = v2 > 0.0
    h[ok] = np.column_stack([vy[ok] / v2[ok], -vx[ok] / v2[ok], np.zeros(ok.sum())])
    u[ok] = np.column_stack([-(vx * vz)[ok] / (v3 * v2)[ok],
                             -(vy * vz)[ok] / (v3 * v2)[ok],
                             (v3 ** 2 - vz ** 2)[ok] / (v3 * v2)[ok]])

    nemo = np.arctan(kp * np.sin(b))
    cn, sn = np.cos(nemo)[:, None], np.sin(nemo)[:, None]

    return cn * h - sn * u, cn * u + sn * h


def _bounded_newton_step(jac, f, x, lo, hi, max_step):
    """Newton step -J^-1 f, with variables that would leave [lo, hi] clamped to the
    bound and the remaining ones re-solved (least squares), then scaled down
    uniformly so no component exceeds max_step."""
    free = np.ones(len(x), dtype=bool)
    dx = np.zeros(len(x))

    for _ in range(len(x)):
        rhs = -(f + jac[:, ~free] @ dx[~free])
        try:
            dx[free] = np.linalg.lstsq(jac[:, free], rhs, rcond=None)[0]
        except np.linalg.LinAlgError:
            dx[free] = 0.0

        x_try = x + dx
        viol = ((x_try < lo) | (x_try > hi)) & free
        if not viol.any():
            break
        dx[viol] = np.clip(x_try[viol], lo[viol], hi[viol]) - x[viol]
        free[viol] = False

    ratio = np.max(np.abs(dx) / np.asarray(max_step, dtype=float))
    if ratio > 1.0:
        dx /= ratio

    return dx


def _measure_trajectory(si, r, v, ctx, shift):
    """Objectives and diagnostics of one tracked test particle."""
    analytic_vars = si.analytic_variables
    numerical_vars = si.numerical_variables

    out = {"failure": None}

    # Electrode ends: the truncated design trajectory is what the electrodes were
    # swept along, so its end points (plus the applied shift) are the end planes.
    trj = analytic_vars["trj_design"] + shift
    vdes = analytic_vars["v_design"]
    p_in, p_out = trj[0], trj[-1]
    n_in = vdes[0] / np.linalg.norm(vdes[0])
    n_out = vdes[-1] / np.linalg.norm(vdes[-1])

    seg = np.linalg.norm(np.diff(r, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])

    d_in = (r - p_in) @ n_in
    d_out = (r - p_out) @ n_out

    _i = np.where(d_in >= 0.0)[0]
    if _i.size == 0:
        out["failure"] = "never reached the entrance plane"
        return out
    i_in = int(_i[0])

    _o = np.where(d_out[i_in:] >= 0.0)[0]
    if _o.size == 0:
        out["failure"] = "never reached the exit plane"
        return out
    i_out = int(_o[0]) + i_in

    def _cross(i, d):
        # path length at which d changes sign between i-1 and i
        if i == 0 or d[i] == d[i - 1]:
            return s[i]
        return s[i - 1] + (s[i] - s[i - 1]) * (-d[i - 1]) / (d[i] - d[i - 1])

    s_in, s_out = _cross(i_in, d_in), _cross(i_out, d_out)
    out["s_in"], out["s_out"] = s_in, s_out
    out["exit_point"] = r[i_out].copy()
    out["exit_angle_deg"] = float(np.degrees(np.arcsin(v[i_out, 2] / np.linalg.norm(v[i_out]))))

    # --- outside: flatness and vertical offset ---
    w0 = s_out + ctx["outside_start"]
    w1 = w0 + ctx["outside_window"]
    win = (s >= w0) & (s <= w1)
    if win.sum() < 5:
        out["failure"] = "track too short for the outside window (raise nsteps)"
        return out

    vnorm = np.linalg.norm(v[win], axis=1)
    out["angle_deg"] = float(np.mean(np.degrees(np.arcsin(np.clip(v[win, 2] / vnorm, -1.0, 1.0)))))
    out["z_offset"] = float(np.mean(r[win, 2]))
    out["outside_window_s"] = (w0, w1)

    e_win = np.linalg.norm(numerical_vars["ef_itp"](r[win]), axis=1)
    out["outside_efield_rel"] = float(e_win.max() / analytic_vars["ef_design"])
    out["outside_ez_rel"] = float(np.abs(numerical_vars["ef_itp"](r[win])[:, 2]).max()
                                 / analytic_vars["ef_design"])
    lim = numerical_vars["limits"]
    out["outside_in_box"] = bool(np.all((r[win] >= lim[:, 0]) & (r[win] <= lim[:, 1])))

    # --- inside: clearance to the electrode surfaces ---
    # Contiguous stretch between the end planes: after the exit the orbit curls
    # round and can come back behind the exit plane, which must not count as inside.
    inside = np.arange(i_in, i_out + 1)
    frac = (s[inside] - s_in) / max(s_out - s_in, 1e-12)

    verts_a, tris_a = _electrode_surface(si, "anode")
    verts_c, tris_c = _electrode_surface(si, "cathode")

    d_a = _surface_distance(r[inside], verts_a, tris_a)
    d_c = _surface_distance(r[inside], verts_c, tris_c)

    out["clearance_anode"] = d_a
    out["clearance_cathode"] = d_c
    out["inside_fraction"] = frac
    out["min_clearance"] = float(min(d_a.min(), d_c.min()))

    cw = ctx["centering_window"]
    sel = (frac >= cw[0]) & (frac <= cw[1])
    if sel.sum() < 3:
        out["failure"] = "too few points in the centering window"
        return out

    # +: closer to the cathode than to the anode
    out["centering"] = float(np.mean(0.5 * (d_a[sel] - d_c[sel])))
    out["clearance_anode_exit"] = float(np.mean(d_a[sel]))
    out["clearance_cathode_exit"] = float(np.mean(d_c[sel]))

    if out["min_clearance"] < ctx["min_clearance"]:
        out["failure"] = "particle within {:.1f} mm of an electrode (min clearance {:.2f} mm)".format(
            1e3 * ctx["min_clearance"], 1e3 * out["min_clearance"])
        return out

    # --- diagnostic: offset from the untruncated design orbit, split into the gap
    # direction (should agree with the clearance difference) and the width direction
    # (unmeasured by the objectives; reported so a drift along the electrodes is seen).
    from scipy.spatial import cKDTree

    trj_full = analytic_vars["trj_design_full"] + shift
    h_dir, u_dir = _optical_frame(analytic_vars["v_design_full"], analytic_vars["b_full"],
                                  analytic_vars["kp"])
    _, j = cKDTree(trj_full).query(r[inside])
    off = r[inside] - trj_full[j]
    out["design_offset_gap"] = np.einsum("ij,ij->i", off, u_dir[j])
    out["design_offset_width"] = np.einsum("ij,ij->i", off, h_dir[j])

    # Same window as the clearance-based value. The sign is flipped so that both
    # centering measures are positive when the particle sits closer to the cathode.
    out["centering_clearance"] = out["centering"]
    out["centering_design"] = -float(np.mean(out["design_offset_gap"][sel]))
    out["width_offset"] = float(np.mean(out["design_offset_width"][sel]))
    if ctx["centering"] == "design":
        out["centering"] = out["centering_design"]

    return out


def _solve_voltage(track, ctx):
    """Voltage scale that zeroes the gap centering, from tracking alone.

    Used when the field is exactly proportional to the spiral-electrode voltage
    (every other conductor at 0 V): track(scale) re-tracks the particle in the
    scaled field. Secant steps, regula falsi once the zero is bracketed. Returns
    (scale, measurements, number of tracks).
    """
    lo, hi = ctx["volt_bounds"]
    tol = ctx["volt_tol_centering"]

    hist = []

    def _try(s):
        s = float(np.clip(s, lo, hi))
        m = track(s)
        hist.append((s, m))
        return s, m

    s0, m0 = _try(ctx["volt_scale0"])
    if m0["failure"] is None and abs(m0["centering"]) <= tol:
        return s0, m0, 1

    _try(s0 + (0.01 if s0 + 0.01 <= hi else -0.01))

    for _ in range(ctx["volt_maxiter"]):
        pts = [(s, m["centering"]) for s, m in hist if m["failure"] is None]
        if len(pts) < 2:
            break

        pos = [p for p in pts if p[1] > 0.0]
        neg = [p for p in pts if p[1] < 0.0]
        if pos and neg:
            sp = min(pos, key=lambda p: p[1])
            sn = max(neg, key=lambda p: p[1])
            s_new = sp[0] - sp[1] * (sn[0] - sp[0]) / (sn[1] - sp[1])
        else:
            (sa, fa), (sb, fb) = pts[-2], pts[-1]
            if fb == fa:
                break
            s_new = sb - fb * (sb - sa) / (fb - fa)

        s_new = float(np.clip(s_new, lo, hi))
        if any(abs(s_new - s) < 1e-5 for s, _ in hist):
            break

        _, m_new = _try(s_new)
        if m_new["failure"] is None and abs(m_new["centering"]) <= tol:
            break

    good = [(s, m) for s, m in hist if m["failure"] is None]
    if not good:
        return hist[0][0], hist[0][1], len(hist)
    s_best, m_best = min(good, key=lambda p: abs(p[1]["centering"]))
    return s_best, m_best, len(hist)


def _evaluate_trajectory(si, x, ctx):
    """Build the geometry for knobs x = (db_entrance [deg], db_exit [deg], dz [m]
    [, volt_scale]), solve, track the test particle and return (residuals,
    measurements)."""
    analytic_vars = si.analytic_variables
    analytic_params = si.analytic_parameters

    t0 = time.time()
    db_ent, db_exit, dz = x[:3]

    # Knob 1 and 2: electrode truncation. The exit is referenced to the b_max the
    # untruncated orbit reaches, so db_exit = 0 means "no truncation".
    si.set_blim(b_min=db_ent, b_max=np.rad2deg(ctx["b_max_achieved"]) - db_exit)

    # Knob 3: vertical position of the inflector assembly. generate_solid_assembly
    # applies shift_applied to the inflector electrodes and shift_applied_lab to the
    # quadrupoles and their apertures.
    shift = np.array([ctx["dx0"], ctx["dy0"], dz], dtype=float)
    si.track_variables["shift_applied"] = shift
    si.track_variables["shift_applied_lab"] = (np.array([0.0, 0.0, dz]) if ctx["move_quadrupoles"]
                                               else np.zeros(3))

    # Knob 4: operating voltage relative to the design voltage. The geometry stays
    # that of the design voltage; only the boundary values change. Modes:
    #   outer  the voltage is x[3], solved by the outer solver together with the
    #          geometry (needed when other conductors are at non-zero voltage)
    #   inner  every other conductor is at 0 V, so the field is exactly proportional
    #          to the voltage: it is computed once at the design voltage, and the
    #          voltage that centres the particle is found by re-tracking alone
    #   fixed  ctx["volt_scale0"]
    mode = ctx["voltage_mode"]
    volt_scale = float(x[3]) if mode == "outer" else ctx["volt_scale0"]
    si.numerical_parameters["volt_scale"] = 1.0 if mode == "inner" else volt_scale

    # The on-axis test particle sees no quadrupole field, so the quadrupoles and
    # their apertures can be left out of the mesh during the optimization.
    quads_saved = si.numerical_parameters["make_quadrupoles"]
    if ctx["exclude_quadrupoles"]:
        si.numerical_parameters["make_quadrupoles"] = False

    try:
        si.generate_meshed_model()
        t_mesh = time.time()
        si.solve()
        t_solve = time.time()

        # Potential on a box around the whole (untruncated) orbit plus a fringe
        # margin, from the launch point up past the exit.
        trj_full = analytic_vars["trj_design_full"] + shift
        m = ctx["fringe_margin"] * (analytic_params["gap"] + 2.0 * analytic_params["dx"])
        lo = np.minimum(trj_full.min(axis=0), 0.0) - m
        hi = np.maximum(trj_full.max(axis=0), 0.0) + m
        limits = ((lo[0], hi[0]), (lo[1], hi[1]), (ctx["z_start"], hi[2]))

        si.calculate_potential(limits=limits, res=ctx["res"], domain_decomp=ctx["domain_decomp"],
                               overlap=0)
        si.calculate_efield()
        t_pot = time.time()

        efield = si.numerical_variables["ef_itp"]
        r_start = np.array([0.0, 0.0, ctx["z_start"]])
        v_start = np.array([0.0, 0.0, ctx["v0"]])

        def _track(scale):
            efield.scaling = scale
            r, v = si.fast_track(r_start=r_start, v_start=v_start,
                                 nsteps=ctx["nsteps"], dt=ctx["dt"])
            meas = _measure_trajectory(si, r, v, ctx, shift)
            meas["r"], meas["v"] = r, v
            return meas

        if mode == "inner":
            volt_scale, meas, n_tracks = _solve_voltage(_track, ctx)
            efield.scaling = volt_scale
            si.numerical_parameters["volt_scale"] = volt_scale
            ctx["volt_scale0"] = volt_scale  # warm start for the next geometry
        else:
            meas = _track(1.0)
            n_tracks = 1
    finally:
        si.numerical_parameters["make_quadrupoles"] = quads_saved

    # Each solve holds a dense N x N matrix (5 GB at 26k triangles). Reference cycles
    # inside the bempp objects can keep the previous one alive until the cyclic
    # collector runs, which it decides by object count, not size.
    import gc
    gc.collect()
    meas["x"] = np.array(x, dtype=float)
    meas["volt_scale"] = volt_scale
    meas["voltage"] = analytic_params["volt"] * volt_scale
    meas["n_tracks"] = n_tracks
    meas["quadrupoles_included"] = bool(quads_saved and not ctx["exclude_quadrupoles"])
    meas["b_lim_deg"] = np.rad2deg(analytic_params["b_lim"]).copy()
    meas["n_triangles"] = int(si.numerical_variables["full mesh"]["elems"].shape[1])
    meas["timing"] = {"mesh": t_mesh - t0, "solve": t_solve - t_mesh,
                      "potential": t_pot - t_solve, "track": time.time() - t_pot,
                      "total": time.time() - t0}

    if meas["failure"] is not None:
        return np.full(4, np.nan), meas

    return np.array([meas["angle_deg"], meas["centering"], meas["z_offset"],
                     meas["width_offset"]]), meas


def optimize_trajectory(si, initial_guess=None, maxiter=15, solver="auto", vary_voltage=True,
                        exclude_quadrupoles=True,
                        tol_angle=0.05, tol_offset=0.2e-3, tol_width=0.5e-3,
                        res=0.005, domain_decomp=(3, 3, 3), z_start=None, dt=1e-10, nsteps=1600,
                        fd_steps=(0.5, 0.5, 1.0e-3, 0.01),
                        bounds=((0.0, 15.0), (0.0, 15.0), (-15.0e-3, 15.0e-3), (0.9, 1.1)),
                        max_step=(3.0, 3.0, 6.0e-3, 0.03),
                        outside_start=None, outside_window=0.025,
                        centering_window=(0.75, 1.0), centering="clearance",
                        fringe_margin=1.5, min_clearance=1.0e-3,
                        move_quadrupoles=True, verbose=True):
    """Set the fringe-field corrections of the spiral inflector from one test particle
    tracked through the whole system.

    Four knobs are solved against four residuals (three knobs, least squares, with
    vary_voltage=False). Each evaluation is a remesh, BEM++ solve, potential grid and
    single-particle track, so expect half a minute to a minute per evaluation; a
    solve from zero takes about 10-15 evaluations.

    Knobs (the solution vector x)
      db_entrance [deg]  entrance truncation, b_lim[0]
      db_exit     [deg]  exit truncation, measured back from the b_max the untruncated
                         design orbit actually reaches (generate_numerical_trajectory
                         stops when vz < 0, typically short of 90 deg), so 0 means
                         "no truncation" and the number reported is the real shortening
      dz          [m]    vertical position of the inflector assembly (the electrodes,
                         housing and entrance aperture); the quadrupoles and their
                         apertures follow axially when move_quadrupoles is True so the
                         quad-to-inflector drift is preserved. The inflector height is
                         a design parameter and is not varied.
      volt_scale  [1]    operating voltage of the spiral electrodes relative to the
                         design voltage (numerical parameter "volt_scale"). The
                         geometry is built for the design voltage; this only changes
                         the boundary values. It is the knob that sets the curvature
                         inside the electrodes and with it the gap-direction
                         centering, which the truncations hardly touch. When every
                         other conductor is at 0 V (the default, with the
                         quadrupoles excluded) the field is exactly proportional to
                         it, so it is solved per geometry by re-tracking alone and the
                         outer solver sees three knobs; otherwise it is a fourth
                         outer knob. Omitted when vary_voltage is False.

    Objectives (residual vector, all driven to zero)
      angle    [deg]  angle between the outgoing velocity and the median plane, averaged
                      over the outside window
      centering [m]   0.5 * (clearance to anode - clearance to cathode), i.e. how far the
                      particle sits from the middle of the gap, averaged over the
                      centering_window fraction of the inside path (default: last quarter).
                      Measured against the actual anode and cathode surface meshes.
      z_offset [m]    mean height of the outside part of the trajectory
      width    [m]    offset from the untruncated design orbit along the electrode
                      width, over the same window. On the HCHC-60 deck this reaches
                      8 mm without entrance truncation and is what the entrance
                      truncation actually controls; the gap-direction centering above
                      hardly responds to any of the three knobs, so without this
                      residual the split between entrance and exit truncation is
                      poorly determined.

    Solvers
      "dfols"    DFO-LS (Cartis et al.), a model-based derivative-free least-squares
                 solver with a trust region and noise handling, which suits expensive
                 residuals that carry remeshing noise. Needs the dfols package.
      "broyden"  built-in: Newton step on a finite-difference Jacobian, then Broyden
                 updates, with bounds, a per-iteration step limit and a fallback to
                 the best point when a step makes things much worse. No dependencies.
      "auto"     dfols if it can be imported, else broyden.
    Residuals are scaled by the tolerances for both, so the objective is balanced.

    The test particle starts at x = y = 0, z = z_start (default: two gaps below the
    entrance), purely along +z, at the mean beam energy. On axis the quadrupole
    field vanishes by symmetry, so it reaches the inflector undeflected. "Outside" is
    the path-length window [outside_start, outside_start + outside_window] past the
    electrode end plane; outside_start defaults to 2 gaps, which is past the housing
    wall and where the fringe field has decayed, and the window ends long before the
    orbit curls back into the housing. The maximum field in the window, relative to
    the design field, is reported as a check.

    No rotation of the inflector is adjusted: the azimuth at which the beam is placed
    in the central region is a design parameter. The entrance angle in the plane of
    the electrode width is not a knob either, since the three exit objectives absorb
    it; the offset along the width is reported as a diagnostic.

    :param si: the spiral inflector object (numerical method, bempp solver)
    :param initial_guess: (db_entrance [deg], db_exit [deg], dz [m], volt_scale); None
                          takes the state the object is in (current b_lim, applied
                          shift and volt_scale)
    :param maxiter: broyden: Newton/Broyden iterations after the Jacobian; dfols:
                    evaluations after its initial sample of n_knobs + 1. 0 evaluates
                    the starting point only and returns.
    :param solver: "auto", "dfols" or "broyden"
    :param vary_voltage: include volt_scale as a knob
    :param exclude_quadrupoles: leave the quadrupoles and their apertures out of the
                                mesh while optimizing. The on-axis test particle sees
                                no quadrupole field, the mesh is a third smaller, and
                                the voltage knob becomes free (see volt_scale). The
                                final rebuild puts them back and reports the residuals
                                of the full assembly.
    :param tol_angle: convergence tolerance on the exit angle (deg)
    :param tol_offset: convergence tolerance on centering and z offset (m)
    :param tol_width: convergence tolerance on the width offset (m). The tolerances
                      also scale the residuals for the least-squares solve.
    :param res: resolution of the potential grid (m); this also sets most of the
                cost per evaluation
    :param domain_decomp: passed to calculate_potential
    :param z_start: launch height of the test particle (m)
    :param dt, nsteps: tracking time step and number of steps
    :param fd_steps: broyden: finite-difference steps for the initial Jacobian, in
                     knob units; dfols: initial trust-region radius per knob
    :param bounds: (lo, hi) per knob; truncations cannot be negative
    :param max_step: broyden: largest change of each knob per iteration
    :param outside_start, outside_window: outside window along the path (m)
    :param centering_window: fraction of the inside path over which centering is
                             averaged
    :param centering: "clearance" measures centering against the anode and cathode
                      surface meshes (equal clearance); "design" measures the offset
                      from the untruncated design orbit along the gap direction, i.e.
                      the centreline the electrodes were built around. Both are
                      reported every evaluation; this picks the one that is solved.
    :param fringe_margin: potential box margin around the orbit, in units of
                          gap + 2 * electrode thickness
    :param min_clearance: an evaluation whose particle comes closer than this to an
                          electrode is treated as failed
    :param move_quadrupoles: give the quadrupoles and their apertures the same dz
    :param verbose: print one line per evaluation
    :return: dict with the solution ("db_entrance", "db_exit", "dz", "b_lim_deg",
             "shift"), the final residuals, the full evaluation history and the
             measurements of the best evaluation. Also stored in
             si.track_variables["optimize_trajectory"].
    """
    assert si.method != "analytical", "You can't optimize using an analytical model!"
    assert si.solver == "bempp", "optimize_trajectory needs the bempp solver"

    analytic_params = si.analytic_parameters
    analytic_vars = si.analytic_variables
    track_vars = si.track_variables

    if solver == "auto":
        try:
            import dfols  # noqa: F401
            solver = "dfols"
        except ImportError:
            solver = "broyden"
    assert solver in ("dfols", "broyden"), "solver must be 'auto', 'dfols' or 'broyden'"

    if analytic_vars.get("trj_design_full") is None or analytic_vars.get("b_max_achieved") is None:
        si.generate_design_trajectory()

    gap = analytic_params["gap"]
    height = analytic_vars["height"]
    ion = analytic_params["ion"]

    # Existing lateral shift is kept; only the axial component is a knob here.
    _prev = track_vars.get("shift_applied")
    dx0, dy0, dz0 = (0.0, 0.0, 0.0) if _prev is None else np.asarray(_prev, dtype=float)
    numerical_params = si.numerical_parameters
    volt_scale0 = float(numerical_params.get("volt_scale", 1.0))
    if initial_guess is not None and len(initial_guess) > 3 and initial_guess[3] is not None:
        volt_scale0 = float(initial_guess[3])

    quads_present = bool(numerical_params.get("make_quadrupoles", False))
    exclude_quadrupoles = bool(exclude_quadrupoles and quads_present)

    # The field is proportional to the spiral-electrode voltage only when every
    # other conductor is at 0 V; then the voltage is solved per geometry by tracking.
    _others_grounded = (
        (not quads_present or exclude_quadrupoles)
        and (not numerical_params.get("make_aperture", False)
             or numerical_params["aperture_params"].get("voltage", 0.0) == 0.0)
        and (not numerical_params.get("make_housing", False)
             or numerical_params["housing_params"].get("voltage", 0.0) == 0.0)
        and (not numerical_params.get("make_cylinder", False)
             or numerical_params["cylinder_params"].get("voltage", 0.0) == 0.0))
    if not vary_voltage:
        voltage_mode = "fixed"
    elif _others_grounded:
        voltage_mode = "inner"
    else:
        voltage_mode = "outer"
    n_x = 4 if voltage_mode == "outer" else 3

    ctx = {"b_max_achieved": analytic_vars["b_max_achieved"],
           "dx0": dx0, "dy0": dy0, "volt_scale0": volt_scale0, "vary_voltage": bool(vary_voltage),
           "voltage_mode": voltage_mode, "exclude_quadrupoles": exclude_quadrupoles,
           "volt_bounds": tuple(bounds[3]) if len(bounds) > 3 else (0.9, 1.1),
           "volt_tol_centering": 0.5 * tol_offset, "volt_maxiter": 6,
           "z_start": -(height + 2.0 * gap) if z_start is None else float(z_start),
           "v0": ion.v_mean_m_per_s,
           "dt": dt, "nsteps": int(nsteps), "res": res, "domain_decomp": domain_decomp,
           "outside_start": 2.0 * gap if outside_start is None else float(outside_start),
           "outside_window": float(outside_window),
           "centering_window": centering_window,
           "centering": centering,
           "fringe_margin": float(fringe_margin),
           "min_clearance": float(min_clearance),
           "move_quadrupoles": bool(move_quadrupoles),
           "solver": solver}

    _b_lim = np.rad2deg(analytic_params["b_lim"])
    x_state = [_b_lim[0], max(0.0, np.rad2deg(ctx["b_max_achieved"]) - _b_lim[1]), dz0,
               volt_scale0]
    if initial_guess is None:
        x = np.array(x_state[:n_x], dtype=float)
    else:
        x = np.array([_g if _g is not None else x_state[k]
                      for k, _g in enumerate(list(initial_guess)[:n_x])], dtype=float)
        if len(x) < n_x:
            x = np.concatenate([x, x_state[len(x):n_x]])

    lo = np.array([b[0] for b in bounds][:n_x], dtype=float)
    hi = np.array([b[1] for b in bounds][:n_x], dtype=float)
    fd_steps = np.asarray(fd_steps, dtype=float)[:n_x]
    max_step = np.asarray(max_step, dtype=float)[:n_x]
    x = np.clip(x, lo, hi)
    tol = np.array([tol_angle, tol_offset, tol_offset, tol_width])
    n_res = len(tol)
    knob_names = ["db_entrance [deg]", "db_exit [deg]", "dz [mm]", "volt_scale"][:n_x]

    print("Starting the full-trajectory optimization ({}, voltage {}, quadrupoles {}): b_max "
          "achieved {:.3f} deg, test particle from z = {:+.1f} mm at {:.4e} m/s, outside window "
          "{:.0f}..{:.0f} mm past the exit.".format(
              solver, voltage_mode, "excluded" if exclude_quadrupoles else "included",
              np.rad2deg(ctx["b_max_achieved"]), 1e3 * ctx["z_start"], ctx["v0"],
              1e3 * ctx["outside_start"], 1e3 * (ctx["outside_start"] + ctx["outside_window"])),
          flush=True)

    history = []

    def _norm(f):
        return float(np.sqrt(np.sum((f / tol) ** 2)))

    def _converged(f):
        return bool(np.all(np.abs(f) <= tol))

    def _eval(xx, tag):
        f, meas = _evaluate_trajectory(si, xx, ctx)
        meas["tag"] = tag
        meas["index"] = len(history)
        meas["residual"] = f
        meas["norm"] = _norm(f) if np.all(np.isfinite(f)) else np.inf
        history.append(meas)

        if verbose:
            _t = meas["timing"]
            _knobs = "ent {:6.3f}  exit {:6.3f} deg  dz {:+7.3f} mm  V x{:.4f}".format(
                xx[0], xx[1], 1e3 * xx[2], meas["volt_scale"])
            if meas["failure"] is None:
                print("[{:>6s} {:2d}] {} | "
                      "angle {:+7.3f} deg  centre {:+6.3f} mm  z {:+6.3f} mm  "
                      "width {:+6.2f} mm | "
                      "clr A/C {:.2f}/{:.2f} mm  design gap {:+.2f} mm  "
                      "E_out {:.1e} (Ez {:.1e}) | "
                      "{:.0f} s (mesh {:.0f}, solve {:.0f}, pot {:.0f}, {} tracks {:.0f})".format(
                          tag, meas["index"], _knobs,
                          f[0], 1e3 * f[1], 1e3 * f[2], 1e3 * f[3],
                          1e3 * meas["clearance_anode_exit"], 1e3 * meas["clearance_cathode_exit"],
                          1e3 * meas["centering_design"],
                          meas["outside_efield_rel"], meas["outside_ez_rel"],
                          _t["total"], _t["mesh"], _t["solve"], _t["potential"],
                          meas["n_tracks"], _t["track"]), flush=True)
            else:
                print("[{:>6s} {:2d}] {} | FAILED: {} | {:.0f} s".format(
                    tag, meas["index"], _knobs, meas["failure"], _t["total"]), flush=True)
            if meas["failure"] is None and not meas["outside_in_box"]:
                print("    warning: part of the outside window lies beyond the potential box "
                      "(field taken as zero there); raise fringe_margin to check.", flush=True)

        return f, meas

    f, meas = _eval(x, "start")
    if not np.all(np.isfinite(f)):
        raise RuntimeError("The starting point could not be evaluated: {}".format(meas["failure"]))

    best = {"x": x.copy(), "f": f.copy(), "norm": _norm(f), "index": meas["index"]}

    def _update_best(xx, ff, idx):
        if np.all(np.isfinite(ff)) and _norm(ff) < best["norm"]:
            best.update({"x": xx.copy(), "f": ff.copy(), "norm": _norm(ff), "index": idx})

    status = "max iterations reached"

    if _converged(f):
        status = "converged at the starting point"

    elif maxiter <= 0:
        status = "starting point evaluated only (maxiter = 0)"

    elif solver == "dfols":
        import dfols

        x_start, f_start = x.copy(), f.copy()

        def _objfun(xx):
            if np.allclose(xx, x_start):
                return f_start / tol  # DFO-LS evaluates x0 itself; reuse the start
            ff, mm = _eval(np.asarray(xx, dtype=float), "dfols")
            _update_best(np.asarray(xx, dtype=float), ff, mm["index"])
            if not np.all(np.isfinite(ff)):
                # A failed evaluation (particle lost, track too short): hand the
                # trust region a large residual so it retreats from that point.
                worst = max([h["norm"] for h in history if np.isfinite(h["norm"])] + [1.0])
                return np.full(n_res, 10.0 * worst)
            return ff / tol

        # scaling_within_bounds maps each knob to [0, 1] over its bounds, so the
        # trust-region radius is set in those units: fd_steps as the initial radius,
        # the useful resolution of each knob as the final one.
        span = hi - lo
        rhobeg = float(np.min(fd_steps / span))
        rhoend = float(np.min(np.array([0.02, 0.02, 0.05e-3, 1e-3])[:n_x] / span))

        soln = dfols.solve(_objfun, x, bounds=(lo, hi),
                           rhobeg=rhobeg, rhoend=rhoend,
                           maxfun=int(maxiter) + n_x + 1,
                           objfun_has_noise=True, scaling_within_bounds=True,
                           user_params={"model.abs_tol": 1.0},  # all residuals within tol
                           print_progress=False)
        status = "dfols: {} (flag {})".format(str(soln.msg).strip(), soln.flag)

    else:
        # Finite-difference Jacobian, one evaluation per knob. Forward steps, unless
        # that would leave the bounds. Four residuals for three knobs, so the Newton
        # step below is a Gauss-Newton (least-squares) step.
        jac = np.zeros((n_res, n_x))
        for k in range(n_x):
            step = fd_steps[k] if x[k] + fd_steps[k] <= hi[k] else -fd_steps[k]
            xk = x.copy()
            xk[k] += step
            fk, mk = _eval(xk, "jac{}".format(k))
            if not np.all(np.isfinite(fk)):
                raise RuntimeError("Jacobian evaluation {} failed: {}".format(k, mk["failure"]))
            jac[:, k] = (fk - f) / step
            _update_best(xk, fk, mk["index"])

        if verbose:
            with np.printoptions(precision=4, suppress=True):
                print("Jacobian (rows: angle [deg], centering [mm], z [mm], width [mm]; "
                      "columns: {}):".format(", ".join(knob_names)))
                print(jac * np.array([[1.0], [1e3], [1e3], [1e3]])
                      * np.array([[1.0, 1.0, 1e-3, 1.0][:n_x]]))

        for it in range(int(maxiter)):
            dx = _bounded_newton_step(jac, f, x, lo, hi, max_step)

            if np.linalg.norm(dx / max_step) < 1e-9:
                status = "no step possible within the bounds"
                break

            f_new, m_new = None, None
            for _try in range(3):
                x_new = np.clip(x + dx, lo, hi)
                f_new, m_new = _eval(x_new, "iter{}".format(it))
                if np.all(np.isfinite(f_new)):
                    break
                dx *= 0.5  # failed evaluation: halve the step and try again

            if not np.all(np.isfinite(f_new)):
                status = "evaluation kept failing after halving the step"
                break

            s_vec = x_new - x
            y_vec = f_new - f
            if s_vec @ s_vec > 0.0:
                jac += np.outer(y_vec - jac @ s_vec, s_vec) / (s_vec @ s_vec)

            _update_best(x_new, f_new, m_new["index"])

            if _norm(f_new) > 2.0 * _norm(f):
                # Much worse: continue from the best point with the updated Jacobian.
                x, f = best["x"].copy(), best["f"].copy()
            else:
                x, f = x_new, f_new

            if _converged(f):
                status = "converged"
                break

    # Leave the object in the best state found: the quadrupoles back in (if the deck
    # has them), the voltage fixed at the value found and the BEM solution at that
    # voltage. With the quadrupoles restored this is also the verification of the
    # solution in the full assembly.
    best_meas_opt = history[best["index"]]
    f_final, final_meas = best["f"], best_meas_opt
    if (voltage_mode == "inner" or exclude_quadrupoles
            or not np.allclose(history[-1]["x"], best["x"])):
        if verbose:
            print("Rebuilding at the best evaluation ({}) with the full assembly, voltage "
                  "x{:.5f}.".format(best["index"], best_meas_opt["volt_scale"]), flush=True)
        ctx["voltage_mode"] = "fixed"
        ctx["volt_scale0"] = float(best_meas_opt["volt_scale"])
        ctx["exclude_quadrupoles"] = False
        f_final, final_meas = _eval(best["x"][:3], "final")
        if not np.all(np.isfinite(f_final)):
            final_meas = best_meas_opt

    x = best["x"]
    b_lim_deg = np.rad2deg(analytic_params["b_lim"]).copy()
    shift = np.array(track_vars["shift_applied"], dtype=float)
    track_vars["shift"] = shift.copy()

    result = {"status": status,
              "db_entrance": float(x[0]), "db_exit": float(x[1]), "dz": float(x[2]),
              "volt_scale": float(final_meas["volt_scale"]), "voltage": float(final_meas["voltage"]),
              "b_lim_deg": b_lim_deg, "b_max_achieved_deg": float(np.rad2deg(ctx["b_max_achieved"])),
              "shift": shift, "shift_lab": np.array(track_vars["shift_applied_lab"], dtype=float),
              "residual": best["f"], "residual_final": f_final,
              "converged": _converged(best["f"]), "voltage_mode": voltage_mode,
              "n_evaluations": len(history),
              "history": history, "measurements": final_meas,
              "measurements_opt": best_meas_opt, "context": ctx}
    track_vars["optimize_trajectory"] = result
    si.track_variables = track_vars

    print("Done optimizing ({}), {} evaluations.".format(status, len(history)))
    print("  entrance truncation {:.4f} deg, exit truncation {:.4f} deg (b_lim = {:.4f}..{:.4f} deg), "
          "dz = {:+.4f} mm, voltage x{:.5f} = {:+.1f} V".format(
              x[0], x[1], b_lim_deg[0], b_lim_deg[1], 1e3 * x[2],
              final_meas["volt_scale"], final_meas["voltage"]))
    print("  optimized:  exit angle {:+.4f} deg, centering {:+.4f} mm, z offset {:+.4f} mm, width "
          "offset {:+.4f} mm".format(best["f"][0], 1e3 * best["f"][1], 1e3 * best["f"][2],
                                    1e3 * best["f"][3]))
    print("  final state ({}): exit angle {:+.4f} deg, centering {:+.4f} mm, z offset {:+.4f} mm, "
          "width offset {:+.4f} mm; clearance at the exit anode {:.3f} / cathode {:.3f} mm, "
          "minimum inside {:.3f} mm".format(
              "quadrupoles included" if final_meas["quadrupoles_included"] else "as optimized",
              f_final[0], 1e3 * f_final[1], 1e3 * f_final[2], 1e3 * f_final[3],
              1e3 * final_meas["clearance_anode_exit"], 1e3 * final_meas["clearance_cathode_exit"],
              1e3 * final_meas["min_clearance"]), flush=True)

    return result
