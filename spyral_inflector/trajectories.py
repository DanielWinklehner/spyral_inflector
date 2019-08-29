from dans_pymodules import *
# import multiprocessing as mp
# import time


def z_rotate(v, angle):
    rot = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                    [np.sin(angle), np.cos(angle), 0.0],
                    [0.0, 0.0, 1.0]])
    if len(v) == 2:
        _v = np.array([v[0], v[1], 1.0])
        return np.matmul(rot, _v[np.newaxis].T)[:2, 0]
    else:
        return np.matmul(rot, v[np.newaxis].T)[:, 0]


def track(si, r_start=None, v_start=None, nsteps=10000, dt=1e-12, omit_b=False, omit_e=False):
    assert (r_start is not None and v_start is not None), "Have to specify r_start and v_start for now!"

    analytic_params = si.analytic_parameters
    numerical_vars = si.numerical_variables
    track_vars = si.track_variables

    if not omit_e:
        if numerical_vars["ef_itp"] is None:
            print("No E-Field has been generated. Cannot track!")
            return 1

    pusher = ParticlePusher(analytic_params["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    if omit_e:
        efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        efield1 = numerical_vars["ef_itp"]  # type: Field

    if omit_b:
        bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        bfield1 = analytic_params["bf_itp"]  # type: Field

    r = np.zeros([nsteps + 1, 3])
    v = np.zeros([nsteps + 1, 3])
    r[0, :] = r_start[:]
    v[0, :] = v_start[:]

    # initialize the velocity half a step back:
    ef = efield1(r[0])
    bf = bfield1(r[0])
    _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

    # Track for n steps
    for i in range(nsteps):
        # print(r[i])

        ef = efield1(r[i])
        bf = bfield1(r[i])

        r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)

    track_vars["trj_tracker"] = r
    si.track_variables = track_vars

    return r, v


def fast_track(si, r_start=None, v_start=None, nsteps=10000, dt=1e-12, omit_b=False, omit_e=False):
    assert (r_start is not None and v_start is not None), "Have to specify r_start and v_start for now!"

    analytic_params = si.analytic_parameters
    numerical_vars = si.numerical_variables
    track_vars = si.track_variables

    if numerical_vars["ef_itp"] is None:
        print("No E-Field has been generated. Cannot track!")
        return 1

    pusher = ParticlePusher(analytic_params["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    if omit_e:
        efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        efield1 = numerical_vars["ef_itp"]  # type: Field

    if omit_b:
        bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        bfield1 = analytic_params["bf_itp"]  # type: Field

    pusher.set_efield(efield1)
    pusher.set_bfield(bfield1)

    r, v = pusher.track(r_start, v_start, nsteps, dt)

    track_vars["trj_tracker"] = r
    si.track_variables = track_vars

    return r, v


def fast_track_with_termination(si, r_start=None, v_start=None,
                                nsteps=10000, dt=1e-12,
                                omit_b=False, omit_e=False):
    assert (r_start is not None and v_start is not None), "Have to specify r_start and v_start for now!"

    analytic_params = si.analytic_parameters
    numerical_vars = si.numerical_variables
    track_vars = si.track_variables

    if numerical_vars["ef_itp"] is None:
        print("No E-Field has been generated. Cannot track!")
        return 1

    pusher = ParticlePusher(analytic_params["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    if omit_e:
        efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        efield1 = numerical_vars["ef_itp"]  # type: Field

    if omit_b:
        bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        bfield1 = analytic_params["bf_itp"]  # type: Field

    pusher.set_efield(efield1)
    pusher.set_bfield(bfield1)
    pusher.set_bds(si.bempp_variables["objects"])  # 'objects' is now a PyElectrodeAssembly

    r, v = pusher.track(r_start, v_start, nsteps, dt)

    track_vars["trj_tracker"] = r
    si.track_variables = track_vars

    return r, v


def modulated_track(cr,
                    r_start=None,
                    v_start=None,
                    nsteps=10000,
                    dt=1e-12,
                    freq=0.0,
                    phase=0.0,
                    omit_b=False,
                    omit_e=False):
    """
    Method for tracking in the central region, applies a cosine modulation to the efield.
    :param cr: CentralRegion object
    :param r_start: Initial position
    :param v_start: Initial velocity
    :param nsteps: Number of steps to track
    :param dt: Time step (s)
    :param freq: Frequency of the efield oscillations (Hz)
    :param phase: Phase offset for the modulation
    :param omit_b: If true, use a zero bfield
    :param omit_e: If true, use a zero efield
    :return: Tuple of trajectory r and v
    """
    assert (r_start is not None and v_start is not None), "Have to specify r_start and v_start for now!"

    analytic_params = cr.analytic_parameters
    numerical_vars = cr.numerical_variables
    track_vars = cr.track_variables

    pusher = ParticlePusher(analytic_params["ion"], "boris")

    if omit_e:
        efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        efield1 = numerical_vars["ef_itp"]  # type: Field

    if omit_b:
        bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        bfield1 = analytic_params["bf_itp"]  # type: Field

    r = np.zeros([nsteps + 1, 3])
    v = np.zeros([nsteps + 1, 3])
    r[0, :] = r_start[:]
    v[0, :] = v_start[:]

    # initialize the velocity half a step back:
    ef = efield1(r[0]) * np.cos(phase)
    bf = bfield1(r[0])
    _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

    # Track for n steps
    for i in range(nsteps):
        ef = efield1(r[i]) * np.cos(2.0 * np.pi * freq * i * dt + phase)
        bf = bfield1(r[i])

        r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)
        # v[i + 1] = cr.dee_crossing(r[i], r[i + 1], v[i + 1], dt * i)

    track_vars["trj_tracker"] = r
    cr.track_variables = track_vars

    return r, v


def term_track(cr, r_init, v_init, starting_angle=0.0,
               maxsteps=20000, dt=1e-11, omit_e=True, omit_b=False):
    """
    Method for tracking in the central region, not meant to use the modulated electric field.
    :param cr: CentralRegion object
    :param r_start: Initial position
    :param v_start: Initial velocity
    :param starting_angle:
    :param maxsteps:
    :param dt:
    :param omit_e:
    :param omit_b:
    :return:
    """
    from .central_region import CRSegment

    ion = cr.analytic_parameters["ion"]

    r_start = r_init
    v_start = v_init

    ra, rb = np.array([0.0, 0.0]), np.array([99.0, 0.0])
    segment = CRSegment(ra=ra, rb=rb)  # Make a segment to check for intersections

    pusher = ParticlePusher(ion, "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    if omit_e:
        efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        efield1 = cr.numerical_variables["ef_itp"]
    if omit_b:
        bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        bfield1 = cr.analytic_parameters["bf_itp"]

    r = np.zeros([maxsteps + 1, 3])
    v = np.zeros([maxsteps + 1, 3])
    r[0, :] = r_start[:]
    v[0, :] = v_start[:]

    # initialize the velocity half a step back:
    ef = efield1(r[0])
    bf = bfield1(r[0])
    _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

    # Track until termination
    intersections = 0
    i = 0
    while intersections == 0 and i < maxsteps:
        ef = efield1(r[i])
        bf = bfield1(r[i])

        r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)
        intersections = segment.check_intersection(r[i, :2], r[i + 1, :2])  # 1 if there's an intersection, 0 otherwise
        i += 1

    return r[:i], v[:i]


def orbit_finder(cr, energy_mev, verbose=False, radial_limit=0.3, fatol=1e-2, xatol=1e-2):
    """
    Use the scipy Nelder-Mead algorithm and particle tracking methods to find static equilibrium orbits.
    :param cr: CentralRegion object
    :param energy_mev: Energy of the ion in MeV
    :param verbose: Print results at the end
    :return: Starting point on x-axis, angle of the velocity w.r.t. the y-axis
    """
    import scipy.optimize
    import scipy.interpolate

    freq_orbit = cr.rf_freq / cr.harmonic  # Ideal frequency of an orbit
    omega_orbit = 2.0 * np.pi * freq_orbit  # Angular orbit frequency

    # Make a new ion with the same type as the CR, with energy_mev
    ion = IonSpecies(cr.analytic_parameters["ion"].name(), energy_mev)
    errors = []  # Running list of the errors from the optimization

    def optimization_function(x):

        if verbose:
            print("Initial R: {:.9f}".format(x[0]))
            print("Initial v angle: {:.9f}".format(x[1]))

        if x[0] > radial_limit:
            return 100 + 100 * x[0]  # TODO: Not sure how well this check works -PW

        r_init = np.array([x[0], 1e-12, 0.0])  # 1e-12 to put it just above the x-axis
        v_angle = x[1]

        v_init = np.array([ion.v_m_per_s() * np.sin(v_angle),
                           ion.v_m_per_s() * np.cos(v_angle),
                           0.0])

        r_final, v_final = term_track(cr,
                                      r_init=r_init,
                                      v_init=v_init,
                                      dt=1e-10)

        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # ax.plot(r_final[:, 0], r_final[:, 1])
        # ax.set_xlim([-0.3, 0.3])
        # ax.set_ylim([-0.3, 0.3])
        # plt.show()
        #
        # theta = np.arctan2(r_final[:, 1], r_final[:, 0])  # Calculate the angle of the particle w.r.t. x-axis
        # theta[theta < 0] += 2 * np.pi  # Change the angles from (-pi, pi) to (0, 2*pi)
        # r_itp = scipy.interpolate.interp1d(theta,  # Interpolate the radius as a function of angle
        #                                    np.sqrt(r_final[:, 0] ** 2 + r_final[:, 1] ** 2),
        #                                    fill_value="extrapolate")
        #
        # # Interpolate the three components of velocity as a function of angle
        # vx_itp = scipy.interpolate.interp1d(theta,
        #                                     v_final[:, 0],
        #                                     fill_value="extrapolate")
        # vy_itp = scipy.interpolate.interp1d(theta,
        #                                     v_final[:, 1],
        #                                     fill_value="extrapolate")
        # vz_itp = scipy.interpolate.interp1d(theta,
        #                                     v_final[:, 2],
        #                                     fill_value="extrapolate")
        #
        # # Velocity at the end of tracking
        # v_end = np.array([vx_itp(2 * np.pi), vy_itp(2 * np.pi), vz_itp(2 * np.pi)])
        #
        # # Calculate the angular difference between start and end velocities
        # v_diff = np.arccos(np.dot(v_init[:2], v_end[:2]) / (np.linalg.norm(v_init[:2]) * np.linalg.norm(v_end[:2])))
        # error = (1e2 * r_itp(2.0 * np.pi) - 1e2 * x[0]) ** 2 + v_diff ** 2  # Calculate error

        v_diff = np.arccos((v_init[0]*v_final[-1, 0] + v_init[1]*v_final[-1, 1]) / (np.linalg.norm(v_init[:2])*np.linalg.norm(v_final[-1, :2])))
        error = (np.linalg.norm(r_final[-1, :2] - r_init[:2]) / np.linalg.norm(r_init[:2]))**2 + \
                (v_diff)**2
        print("Error: {:.5f}".format(error))
        # TODO: Look into how the interpolation affects the results, rather than just getting v[-1]
        errors.append(error)

        return error

    sol = scipy.optimize.minimize(optimization_function,
                                  method="Nelder-Mead",
                                  x0=np.array([1.15 * ion.v_m_per_s() / omega_orbit, 0.0]),  # 1.15 is a fudge factor
                                  options={'fatol': fatol, 'xatol': xatol})

    if verbose:
        print("Calculated orbit for {:.3f} MeV <<<".format(energy_mev))
        print(sol.x)

    return sol.x


def generate_analytical_trajectory(si):
    if not si.initialized:
        si.initialize()

    analytic_params = si.analytic_parameters
    analytic_vars = si.analytic_variables

    print("Calculating Design Trajectory... ", end="")

    h = analytic_vars["height"]
    tilt = analytic_params["tilt"]  # type: float
    ion_vel = analytic_params["ion"].v_m_per_s()

    analytic_vars["kp"] = np.tan(np.deg2rad(tilt))  # Tilt parameter
    analytic_vars["k"] = ((h / analytic_vars["r_cyc"]) + analytic_vars["kp"]) / 2.0

    cp = analytic_vars["c+"] = (2.0 * analytic_vars["k"] + 1.0)
    cm = analytic_vars["c-"] = -(2.0 * analytic_vars["k"] - 1.0)

    # --- Trajectory coordinates --- #
    _x = +0.5 * h * ((2.0 / (1.0 - (4.0 * (analytic_vars["k"] ** 2.0)))) -
                     (np.cos(cp * analytic_vars["b"]) / cp) - np.cos(
                -cm * analytic_vars["b"]) / cm)
    _vx = 0.5 * ion_vel * (np.sin((2 * analytic_vars["k"] + 1) * analytic_vars["b"]) +
                           np.sin(-(2 * analytic_vars["k"] - 1) * analytic_vars["b"]))

    _y = -0.5 * h * (np.sin(cp * analytic_vars["b"]) / cp +
                     np.sin(-cm * analytic_vars["b"]) / cm)
    _vy = 0.5 * ion_vel * (np.cos((2 * analytic_vars["k"] + 1) * analytic_vars["b"]) -
                           np.cos(-(2 * analytic_vars["k"] - 1) * analytic_vars["b"]))

    _z = - h * (1.0 - np.sin(analytic_vars["b"]))
    _vz = ion_vel * np.sqrt(1 - np.sin(analytic_vars["b"]) ** 2)

    analytic_vars["trj_design"] = np.array([_x, _y, _z]).T
    analytic_vars["trj_vel"] = np.array([_vx, _vy, _vz]).T

    # Rotation/flip
    if not ((analytic_vars["bf_design"] < 0.0) ^ (analytic_params["ion"].q() < 0.0)):
        if si.debug:
            print("Flipping direction of cyclotron motion...", end="")
        analytic_vars["trj_design"][:, 1] = -analytic_vars["trj_design"][:, 1]

    # Orbit center calculation
    xc, yc = calculate_orbit_center(analytic_vars["k"], analytic_vars["kp"], analytic_vars["height"])

    analytic_vars["orbit_center"] = (xc, yc)

    print("Done!")

    if si.debug:
        print("Design Trajectory:")
        print(analytic_vars["trj_design"])
        print("")

    si.analytic_parameters = analytic_params
    si.analytic_variables = analytic_vars

    return analytic_vars["trj_design"]


def generate_numerical_trajectory(si, bf=None, nsteps=15000, dt=1e-11):
    # TODO: Make sure the nsteps and dt are being consistent throughout the code

    analytic_params = si.analytic_parameters
    analytic_vars = si.analytic_variables
    track_params = si.track_parameters

    if "nsteps" in track_params:
        nsteps = track_params["nsteps"]
    if "dt" in track_params:
        dt = track_params["dt"]

    pusher = ParticlePusher(analytic_params["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    tilt = analytic_params["tilt"]  # type: float

    analytic_vars["kp"] = np.tan(np.deg2rad(tilt))  # Tilt parameter
    analytic_vars["k"] = ((analytic_vars["height"] / analytic_vars["r_cyc"]) + analytic_vars["kp"]) / 2.0

    r_start = np.array([0.0, 0.0, -analytic_vars["height"]])
    v_start = np.array([0.0, 0.0, analytic_params["ion"].v_m_per_s()])

    _r = np.zeros([nsteps + 1, 3])
    _v = np.zeros([nsteps + 1, 3])
    _b = np.zeros([nsteps + 1])  # Store the "b" angle for geometry generation
    _r[0, :] = r_start[:]
    _v[0, :] = v_start[:]

    # Create a new electric field, which will be repeatedly re-defined
    field_val = analytic_vars["ef_design"] * np.sign(analytic_params["ion"].q())
    efield1 = Field(dim=0, field={"x": field_val, "y": 0.0, "z": 0.0})

    if bf is not None:
        bfield1 = bf
    else:
        bfield1 = analytic_params["bf_itp"]

    # initialize the velocity half a step back:
    ef = efield1(_r[0])
    bf = bfield1(_r[0])
    _, _v[0] = pusher.push(_r[0], _v[0], ef, bf, -0.5 * dt)

    i = 0
    while i < nsteps:
        ef = efield1(_r[i])
        bf = bfield1(_r[i])

        _r[i + 1], _v[i + 1] = pusher.push(_r[i], _v[i], ef, bf, dt)

        vx, vy, vz = _v[i + 1]
        vo = np.sqrt(vx ** 2.0 + vy ** 2.0 + vz ** 2.0)
        _b[i + 1] = i * dt * vo / analytic_vars["height"]

        # Toprek theory with surgery
        Eh = field_val * analytic_vars["kp"] * np.sin(_b[i + 1])
        Ehx = -Eh * vy / (np.sqrt(vo ** 2.0 - vz ** 2.0))
        Ehy = Eh * vx / (np.sqrt(vo ** 2.0 - vz ** 2.0))

        ex = field_val * vx * np.abs(vz) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0)) + Ehx
        ey = field_val * vy * np.abs(vz) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0)) + Ehy
        ez = -field_val * (vo ** 2.0 - vz ** 2.0) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0))
        efield1 = Field(dim=0, field={"x": ex, "y": ey, "z": ez})
        if vz < 0:  # Stop when the z-component of the velocity is zero
            if si.debug:
                print(_r[i + 1, :])  # Print the final position
            break
        i += 1
    ns = i

    try:
        i_init = np.where(_b >= analytic_params["b_lim"][0])[0][0]
        if i_init == 0:
            i_init = 1
    except IndexError:
        i_init = 1

    try:
        i_final = np.where(_b >= analytic_params["b_lim"][1])[0][0]
    except IndexError:
        i_final = i

    if si.debug:
        print("Design Trajectory: Initial index: %i, Final index: %i.".format(i_init, i_final))

    # The arrays cannot be as large as they're made initially.
    # This would cause the BEMPP routines to perform the computations
    # with incredibly high resolution and would never finish. -PW

    from scipy.interpolate import interp1d
    _r, _v, _b = _r[i_init:i_final, :], _v[i_init:i_final, :], _b[i_init:i_final]

    _b_itp = interp1d(_v[:, 2], _b, fill_value='extrapolate')
    _b_end = _b_itp(0.0)

    _z_wrt_v = interp1d(_v[:, 2], _r[:, 2], fill_value='extrapolate')
    shift = np.array([0.0, 0.0, - _z_wrt_v(0.0)])
    height = np.abs(_r[0, 2] - _z_wrt_v(0.0))

    _vz = interp1d(_b, _v[:, 2], fill_value='extrapolate')

    _x, _y, _z = interp1d(_b, _r[:, 0], fill_value='extrapolate'), \
                 interp1d(_b, _r[:, 1], fill_value='extrapolate'), \
                 interp1d(_b, _r[:, 2], fill_value='extrapolate')

    _vx, _vy, _vz = interp1d(_b, _v[:, 0], fill_value='extrapolate'), \
                    interp1d(_b, _v[:, 1], fill_value='extrapolate'), \
                    interp1d(_b, _v[:, 2], fill_value='extrapolate')

    b = np.linspace(0.0, _b_end, analytic_params["ns"])

    r = np.array([_x(b), _y(b), _z(b)]).T
    v = np.array([_vx(b), _vy(b), _vz(b)]).T

    r += shift  # TODO: Make sure this doesn't mess up shifts anywhere else (i.e. py_electrodes)

    analytic_vars["trj_design"] = r
    analytic_vars["trj_vel"] = v

    si.analytic_parameters = analytic_params
    si.analytic_variables = analytic_vars

    return r, v


def calculate_orbit_center(k, kp, height):
    xc = height * ((1 - 2 * k * np.sin(k * np.pi)) / (1 - 4 * k ** 2) - np.sin(k * np.pi) / (2 * k - kp))
    yc = height * (2 * k / (1 - 4 * k ** 2) + 1 / (2 * k - kp)) * np.cos(k * np.pi)

    return xc, yc


# not_so_simple_tracker
def simple_tracker(cr, r_start=None, v_start=None, dt=1e-11, nturns=1, phase=0.0):
    from .central_region import TwoDeeField, Sectors

    gaps = cr._abstract_dees
    sectors = Sectors(gaps)

    df = TwoDeeField(left_voltage=0.0, right_voltage=70.0e3)
    dee_field = df._efield

    def get_dee_field(pos):
        prev_seg, next_seg = sectors.get_sector(pos)  # Get the segments before and after the particle pos

        # Get the vector between rb and ra
        dr1 = prev_seg.rb - prev_seg.ra
        dr2 = next_seg.rb - next_seg.ra

        # Calculate the distance from pos to this edge
        d1 = calculate_distance_to_edge(pos - prev_seg.ra, dr1)
        d2 = calculate_distance_to_edge(pos - next_seg.ra, dr2)

        # These fields are x,y: x = direction normal to dee gap, y = z
        # Includes the inherent phase shifts from dee-->dummy dee and vice versa
        # ef2 will always be -distance, since it's approached from the left in those coordinates
        ef1 = prev_seg.phase_shift * dee_field(np.array([d1, pos[2], 0.0]))
        ef2 = next_seg.phase_shift * dee_field(np.array([-d2, pos[2], 0.0]))

        # Calculate the angle that these fields need to be rotated
        rot1 = np.arctan2(dr1[1], dr1[0]) + np.pi / 2.0
        rot2 = np.arctan2(dr2[1], dr2[0]) + np.pi / 2.0

        # Calculate the rotated fields
        rot_ef1 = np.array([ef1[0] * np.cos(rot1), ef1[0] * np.sin(rot1), ef1[1]])
        rot_ef2 = np.array([ef2[0] * np.cos(rot2), ef2[0] * np.sin(rot2), ef2[1]])

        return rot_ef1, rot_ef2

    def calculate_distance_to_edge(pos, edge):
        n_vec = np.array([-edge[1], edge[0]])
        d = np.abs(np.dot(pos[:2], n_vec)) / np.linalg.norm(n_vec)
        return d

    analytic_params = cr.analytic_parameters
    track_vars = cr.track_variables

    pusher = ParticlePusher(analytic_params["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    bfield1 = analytic_params["bf_itp"]  # type: Field

    omega_rf = 2.0 * np.pi * cr.rf_freq
    rf_phase = phase
    maxsteps = int(4.0 * nturns * 2.0 * np.pi / (dt * omega_rf))

    r = np.zeros([maxsteps + 1, 3])
    v = np.zeros([maxsteps + 1, 3])

    r[0, :] = r_start[:]
    v[0, :] = v_start[:]

    # initialize the velocity half a step back:
    ef1, ef2 = get_dee_field(r[0])
    ef = ef1 + ef2
    ef *= np.cos(rf_phase)
    bf = bfield1(r[0])
    _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

    _ef1, _ef2 = [], []

    i = 0
    while i < maxsteps:
        ef1, ef2 = get_dee_field(r[i, :])

        ef = ef1 + ef2
        ef *= np.cos(omega_rf * i * dt + rf_phase)

        bf = bfield1(r[i])

        r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)

        i += 1

    track_vars["trj_tracker"] = r
    cr.track_variables = track_vars

    return r, v
