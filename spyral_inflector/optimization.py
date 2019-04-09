from dans_pymodules import *


def optimize_fringe(si, initial_guess=(None, None), maxiter=10, tol=1e-1, res=0.002):
    """
    This function optimizes the length of the spiral inflector to adjust for fringe fields.
    :param initial_guess: tuple, list or ndarray containing angle adjustment for entrance and exit
                          units: degrees
    :param maxiter: Maximum number of iterations to find the optimum entrance and exit adjustment before giving up
    :param tol: maximum deviation tolerated for a successful solution (degrees)
    :param res:
    :return:
    """

    assert si._method is not "analytical", "You can't optimize using an analytical model!"

    print("Starting the optimization process...")

    # Start with an initial guess for bmin, bmax. If None, use 0.0
    # TODO: Update to use current internal adjustment, if present -DW
    # TODO: Fix the angle in the y direction!
    _db_init = np.array(initial_guess)
    _db_init[np.where(_db_init == None)] = 0.0
    _db = _db_init[:]

    # Half the length of the cube for electric field calculation
    _hcl = (si._params_analytic["gap"] + 2.0 * si._params_analytic["dx"])
    # _ion = si._params_analytic["ion"]  # type: IonSpecies
    x_axis = Vector([1.0, 0.0, 0.0])
    y_axis = Vector([0.0, 1.0, 0.0])
    z_axis = Vector([0.0, 0.0, 1.0])

    _r = np.zeros([1, 3])

    # --- Optimize entrance fringe --- #
    deviation_x = 1e20  # Some large number as initialization
    it = 0

    if si._variables_track["shift"] is None:
        si._variables_track["shift"] = np.zeros(3)

    while abs(deviation_x) > tol:

        # Apply the angle correction (first time use initial guess)
        si._set_blim(b_min=0.0 + _db[0])
        # (Re-)calculate the new geometry and BEM++ solution
        si.generate_meshed_model()
        si.solve_bempp()

        # Trajectory starting point:
        _trj = si._variables_analytic["trj_design"]  # type: np.ndarray
        _v_des = si._variables_analytic["v_design"]  # type: np.ndarray

        # Calculate E-Field upwards of starting point
        xs, ys, zs = _trj[0]

        si.calculate_efield(limits=((xs - 0.5 * _hcl, xs + 0.5 * _hcl),
                                    (ys - 0.5 * _hcl, ys + 0.5 * _hcl),
                                    (zs - 1.9 * _hcl, zs + 0.1 * _hcl)),
                            res=res,
                            domain_decomp=(4, 4, 4),
                            overlap=0)

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

        si._variables_track["shift"][0] -= _r[-1][0]
        si._variables_track["shift"][1] -= _r[-1][1]

        it += 1

        if it == maxiter:
            print("Entrance Fringe: Maximum number of iterations has been reached. Breaking.")
            break

        if it == 1:
            _db[0] += deviation_x
        else:
            _db[0] += 0.5 * deviation_x  # Dampen the oscillations a bit

    # --- Optimize exit fringe --- #
    deviation = 1e20  # Some large number as initialization
    it = 0
    while abs(deviation) > tol:

        # Apply the angle correction (first time use initial guess)
        si._set_blim(b_max=90.0 + _db[1])
        # (Re-)calculate the new geometry and BEM++ solution
        si.generate_meshed_model()
        si.solve_bempp()

        # Trajectory starting point:
        _trj = si._variables_analytic["trj_design"]  # type: np.ndarray
        _v_des = si._variables_analytic["v_design"]  # type: np.ndarray

        _ns = si._params_analytic["ns"]  # type: int
        start_idx = int(0.9 * _ns)  # start the tracking "10 percent" into the spiral inflector exit
        xs, ys, zs = rs = _trj[start_idx]
        vs = _v_des[start_idx]

        # Calculate E-Field
        # TODO: Better way to determine the fringe field region
        si.calculate_efield(limits=((xs - 2.0 * _hcl, xs + 2.0 * _hcl),
                                    (ys - 2.0 * _hcl, ys + 2.0 * _hcl),
                                    (zs - _hcl, zs + _hcl)),
                            res=res,
                            domain_decomp=(4, 4, 4),
                            overlap=0)

        _r, _v = si.track(r_start=rs,  # Use point close to exit along design particle
                          v_start=vs,  # Regular direction of design particle
                          nsteps=5000,
                          dt=1e-11,
                          omit_b=False)

        trj_dir = Vector(_r[-1] - _r[-100])
        deviation = 90.0 - np.rad2deg(trj_dir.angle_with(z_axis))

        print("Current exit adjustment: {:.4f}, "
              "deviation from xy-plane: {:.4f} degrees".format(_db[1], deviation))

        si._variables_track["shift"][2] -= _r[-1, 2]

        it += 1

        if it == maxiter:
            print("Exit Fringe: Maximum number of iterations has been reached. Breaking.")
            break

        if it == 1:
            _db[1] += deviation
        else:
            _db[1] += 0.65 * deviation  # Dampen the oscillations a bit

    # Recalculate the new geometry and BEM++ solution one last time
    si.initialize()
    si.generate_meshed_model()
    si.solve_bempp()

    print("Done optimizing!")

    return si._variables_track["shift"]