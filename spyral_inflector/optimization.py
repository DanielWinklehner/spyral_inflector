# from .global_variables import *
from .vector import Vector
import numpy as np


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
        # Lateral only. The request is a centred ENTRANCE trajectory and a LEVEL exit;
        # levelness is the b_max angle handled above, and translating the inflector
        # vertically would just move it out of the median plane. shift[2] is still
        # measured and returned, it is simply not applied here.
        _applied = np.array(track_vars["shift"], dtype=float)
        _applied[2] = 0.0

        si.track_variables["shift_applied"] = _applied

        print("Applying centering shift dx={:.4f} mm, dy={:.4f} mm "
              "(dz={:.4f} mm measured, not applied)".format(
                  1000.0 * _applied[0], 1000.0 * _applied[1],
                  1000.0 * track_vars["shift"][2]))

    # Recalculate the new geometry and BEM++ solution one last time
    si.initialize()
    si.generate_meshed_model()
    si.solve()

    print("Done optimizing!")

    si.track_variables = track_vars

    return track_vars["shift"]
