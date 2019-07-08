from dans_pymodules import *
# import multiprocessing as mp
import time


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
    # TODO: For now break if r_start or v_start are not given, later get from class properties?
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
    # TODO: For now break if r_start or v_start are not given, later get from class properties?
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
    # TODO: For now break if r_start or v_start are not given, later get from class properties?
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


def central_region_track(cr,
                         r_start=None,
                         v_start=None,
                         nsteps=10000,
                         dt=1e-12,
                         freq=0.0,
                         phase=0.0,
                         omit_b=False,
                         omit_e=False):
    assert (r_start is not None and v_start is not None), "Have to specify r_start and v_start for now!"

    analytic_params = cr.analytic_parameters
    numerical_vars = cr.numerical_variables
    track_vars = cr.track_variables

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


def track_segment(cr,
                  initial_segment=None,
                  end_segment=None,
                  r_start=None,
                  v_start=None,
                  dt=1e-11,
                  maxsteps=10000,
                  omit_b=False,
                  omit_e=True):
    # Calculate gap field from theta_seg, theta_rf
    # Track through gap field
    # Initialize r0, v0
    # Calculate r(i), v(i)          <--|
    # Intersection with segment?  No --|
    # --> Yes: Is current RF within tolerance?
    # If yes: Stop
    # If no: Check if early or late
    # If early: Decrease theta_seg
    # If late: increase theta_seg

    # seg1_field = initial_segment.get_field()

    assert (r_start is not None and v_start is not None), "Have to specify r_start and v_start for now!"
    assert end_segment is not None, "Needs a segment to terminate on!"

    analytic_params = cr.analytic_parameters
    numerical_vars = cr.numerical_variables
    track_vars = cr.track_variables

    pusher = ParticlePusher(analytic_params["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    if omit_e:
        efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        efield1 = numerical_vars["ef_itp"]  # type: Field

    if omit_b:
        bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        bfield1 = analytic_params["bf_itp"]  # type: Field

    r = np.zeros([maxsteps + 1, 3])
    v = np.zeros([maxsteps + 1, 3])
    r[0, :] = r_start[:]
    v[0, :] = v_start[:]

    # initialize the velocity half a step back:
    ef = efield1(r[0])
    bf = bfield1(r[0])
    _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

    # Track for n steps
    intersected = False
    i = 0
    while not intersected and i < maxsteps:
        ef = efield1(r[i])
        bf = bfield1(r[i])

        r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)
        intersected = end_segment.check_intersection(r[i, :2], r[i + 1, :2])
        # v[i + 1] = cr.dee_crossing(r[i], r[i + 1], v[i + 1], dt * i)
        i += 1

    track_vars["trj_tracker"] = r[:i, :]
    cr.track_variables = track_vars

    return r[:i], v[:i]


def cr_track(ion, r_init, v_init, symmetry="full", starting_angle=0.0,
             end_type="termination", nterms=1, input_bfield=None,
             maxsteps=20000, dt=1e-11, omit_e=True, omit_b=False):
    from .central_region import CRSegment

    r_start = r_init
    v_start = v_init

    ra, rb = np.array([0.0, 0.0]), np.array([99.0, 0.0])

    if symmetry == "quarter":
        segment = CRSegment(ra=ra, rb=z_rotate(rb, starting_angle + np.pi / 2.0))
    elif symmetry == "half":
        segment = CRSegment(ra=ra, rb=z_rotate(rb, starting_angle + np.pi))
    else:
        segment = CRSegment(ra=ra, rb=z_rotate(rb, starting_angle))
        # segment = CRSegment(ra=np.array([0.0, 0.0]), rb=np.array([0.0, 5.0]))

    pusher = ParticlePusher(ion, "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    if omit_e:
        efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        efield1 = None

    if omit_b:
        bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        bfield1 = input_bfield

    r = np.zeros([maxsteps + 1, 3])
    v = np.zeros([maxsteps + 1, 3])
    r[0, :] = r_start[:]
    v[0, :] = v_start[:]

    # initialize the velocity half a step back:
    ef = efield1(r[0])
    bf = bfield1(r[0])
    _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

    # Track for n steps
    if end_type == "steps":

        pusher.set_efield(efield1)
        pusher.set_bfield(bfield1)

        r, v = pusher.track(r_start, v_start, maxsteps, dt)

        return r, v

    # Track until termination
    elif end_type == "termination":
        intersections = 0
        i = 0
        while intersections < nterms and i < maxsteps:
            ef = efield1(r[i])
            bf = bfield1(r[i])

            r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)
            intersections += segment.check_intersection(r[i, :2], r[i + 1, :2])
            # v[i + 1] = cr.dee_crossing(r[i], r[i + 1], v[i + 1], dt * i)
            i += 1

        # print("Tracked through {} steps.".format(i))
        # track_vars["trj_tracker"] = r[:i, :]
        # cr.track_variables = track_vars

        return r[:i], v[:i]

    return 0


def orbit_finder(cr, energy_mev, verbose=False):
    import scipy.optimize
    import scipy.interpolate

    freq_orbit = cr.rf_freq / cr.harmonic
    omega_orbit = 2.0 * np.pi * freq_orbit
    ion = IonSpecies("H2_1+", energy_mev)  # TODO: Generalize -PW

    errors = []

    def optimization_function(x):

        if verbose:
            print("Initial R: {:.9f}".format(x[0]))
            print("Initial v angle: {:.9f}".format(x[1]))

        if x[0] > 0.3:
            return 100 + 100*x[0]

        r_init = np.array([x[0], 1e-12, 0.0])
        v_angle = x[1]

        v_init = np.array([ion.v_m_per_s() * np.sin(v_angle),
                           ion.v_m_per_s() * np.cos(v_angle),
                           0.0])

        r_final, v_final = cr_track(ion,
                                    r_init=r_init,
                                    v_init=v_init,
                                    symmetry="full",
                                    dt=1e-10,
                                    nterms=1,
                                    input_bfield=cr.analytic_parameters["bf_itp"])

        theta = np.arctan2(r_final[:, 1], r_final[:, 0])
        theta[theta < 0] += 2 * np.pi
        r_itp = scipy.interpolate.interp1d(theta,
                                           np.sqrt(r_final[:, 0] ** 2 + r_final[:, 1] ** 2),
                                           fill_value="extrapolate")

        vx_itp = scipy.interpolate.interp1d(theta,
                                            v_final[:, 0],
                                            fill_value="extrapolate")
        vy_itp = scipy.interpolate.interp1d(theta,
                                            v_final[:, 1],
                                            fill_value="extrapolate")
        vz_itp = scipy.interpolate.interp1d(theta,
                                            v_final[:, 2],
                                            fill_value="extrapolate")

        v_end = np.array([vx_itp(2*np.pi), vy_itp(2*np.pi), vz_itp(2*np.pi)])
        v_diff = np.arccos(np.dot(v_init[:2], v_end[:2]) / (np.linalg.norm(v_init[:2]) * np.linalg.norm(v_end[:2])))
        error = (1e2 * r_itp(2.0 * np.pi) - 1e2 * x[0]) ** 2 + (v_diff) ** 2

        errors.append(error)

        return error

    sol = scipy.optimize.minimize(optimization_function,
                                  method="Nelder-Mead",
                                  x0=np.array([1.15*ion.v_m_per_s() / omega_orbit, 0.0]),
                                  options={'fatol': 1e-5, 'xatol': 1e-5})

    if verbose:
        print(">>> Finished {:.3f} MeV <<<".format(energy_mev))
        print(sol.x)
    return sol.x


def gordon_algorithm(cr,
                     energy_mev=1.0,
                     symmetry_mode="full",
                     starting_angle=0.0,
                     verbose=False,
                     tune_calc=False,
                     retrack=False):
    import scipy.interpolate

    if symmetry_mode == "full":
        theta_range = 2 * np.pi
    elif symmetry_mode == "half":
        theta_range = np.pi
    elif symmetry_mode == "quarter":
        theta_range = 0.5 * np.pi
    else:
        print("Symmetry mode not recognized.")
        return 1

    freq_orbit = cr.rf_freq / cr.harmonic
    omega_orbit = 2.0 * np.pi * freq_orbit
    ion = IonSpecies("H2_1+", energy_mev)

    a = (clight / omega_orbit) * 1000  # mm
    # b = omega_orbit / ion.q_over_m()

    total_momentum = ion.gamma() * ion.mass_kg() * ion.v_m_per_s()

    cx = 1000.0  # Conversion from m to mm
    cpx = a / (ion.mass_kg() * clight)  # Conversion from SI to Gordon's units

    dri, dpri = 0.0, 0.0
    initial_r = ion.beta() * a / cx
    initial_pr = 0.0
    err = 1
    it = 0

    while err > 1e-6:
        print("Iteration {}:".format(it))
        initial_r += dri / cx
        # print("* Initial r = {}".format(initial_r))
        initial_pr += dpri / cpx
        # print("* Intial pr = {}".format(initial_pr))

        r_init = z_rotate(np.array([initial_r, 1E-9, 0.0]), starting_angle)
        pr_init = z_rotate(np.array([initial_pr, np.sqrt(total_momentum ** 2 - initial_pr ** 2), 0.0]), starting_angle)

        # r_init = np.array([-1E-9, initial_r, 0.0])
        # pr_init = np.array([-np.sqrt(total_momentum ** 2 - initial_pr ** 2), initial_pr, 0.0])

        # print("pr_init", pr_init)
        if verbose:
            print("* Trial Orbit *")
            print("-> Initial r: {}".format(r_init))
            print("-> Initial pr: {}".format(pr_init))
        vr_init = pr_init / (ion.gamma() * ion.mass_kg())

        # Trial orbit
        r, vr = cr_track(ion=ion,
                         r_init=r_init,
                         v_init=vr_init,
                         symmetry=symmetry_mode,
                         starting_angle=starting_angle,
                         end_type="termination",
                         input_bfield=cr._params_analytic["bf_itp"],
                         maxsteps=100000,
                         dt=5e-11,
                         omit_e=True,
                         omit_b=False)

        r_mag = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2)

        pr = ion.gamma() * ion.mass_kg() * vr
        # pr_mag = np.sum(r[:, :2] * pr[:, :2] / np.linalg.norm(r[:, :2]), axis=1)

        pr_mag = []
        for i in range(len(r_mag)):
            rhat = r[i, :2] / np.linalg.norm(r[i, :2])
            _p = np.dot(pr[i, :2], rhat)
            pr_mag.append(_p)
        pr_mag = np.array(pr_mag)

        r1_init = r_init + z_rotate(np.array([0.001, 0.0, 0.0]), starting_angle)  # Offset by 1 mm
        # r1_init = r_init + np.array([0.0, 0.001, 0.0])  # Offset by 1 mm
        initial_pr1 = initial_pr + 0.0
        # pr1_init = np.array([initial_pr1, np.sqrt(total_momentum ** 2 - initial_pr1 ** 2), 0.0])
        pr1_init = pr_init

        if verbose:
            print("* X-Offset Orbit *")
            print("-> Initial r: {}".format(r1_init))
            print("-> Initial pr: {}".format(pr1_init))
        vr1_init = pr1_init / (ion.gamma() * ion.mass_kg())

        r1, vr1 = cr_track(ion=ion,
                           r_init=r1_init,
                           v_init=vr1_init,
                           symmetry=symmetry_mode,
                           starting_angle=starting_angle,
                           end_type="termination",
                           input_bfield=cr._params_analytic["bf_itp"],
                           maxsteps=100000,
                           dt=5e-11,
                           omit_e=True,
                           omit_b=False)

        r1_mag = np.sqrt(r1[:, 0] ** 2 + r1[:, 1] ** 2)
        pr1 = ion.gamma() * ion.mass_kg() * vr1
        # pr1_mag = np.sum(r1[:, :2] * pr1[:, :2] / np.linalg.norm(r1[:, :2]), axis=1)

        pr1_mag = []

        for i in range(len(r1_mag)):
            rhat = r1[i, :2] / np.linalg.norm(r1[i, :2])
            _p = np.dot(pr1[i, :2], rhat)
            pr1_mag.append(_p)
        pr1_mag = np.array(pr1_mag)

        r2_init = r_init
        # initial_pr2 = initial_pr + total_momentum * 1E-4
        initial_pr2 = initial_pr + ion.mass_kg() * clight / a  # Offset by one unit of momentum in Gordon's units
        pr2_init = z_rotate(np.array([initial_pr2, np.sqrt(total_momentum ** 2 - initial_pr2 ** 2), 0.0]),
                            starting_angle)
        # pr2_init = np.array([-np.sqrt(total_momentum ** 2 - initial_pr2 ** 2), initial_pr2, 0.0])

        if verbose:
            print("* PX-Offset Orbit *")
            print("-> Initial r: {}".format(r2_init))
            print("-> Initial pr: {}".format(pr2_init))

        vr2_init = pr2_init / (ion.gamma() * ion.mass_kg())

        r2, vr2 = cr_track(ion=ion,
                           r_init=r2_init,
                           v_init=vr2_init,
                           symmetry=symmetry_mode,
                           starting_angle=starting_angle,
                           end_type="termination",
                           input_bfield=cr._params_analytic["bf_itp"],
                           maxsteps=100000,
                           dt=5e-11,
                           omit_e=True,
                           omit_b=False)

        r2_mag = np.sqrt(r2[:, 0] ** 2 + r2[:, 1] ** 2)
        pr2 = ion.gamma() * ion.mass_kg() * vr2
        # pr2_mag = np.sum(r2[:, :2] * pr2[:, :2] / np.linalg.norm(r2[:, :2]), axis=1)

        pr2_mag = []
        for i in range(len(r2_mag)):
            rhat = r2[i, :2] / np.linalg.norm(r2[i, :2])
            _p = np.dot(pr2[i, :2], rhat)
            pr2_mag.append(_p)
        pr2_mag = np.array(pr2_mag)

        theta = np.arctan2(r[:, 1], r[:, 0])
        theta1 = np.arctan2(r1[:, 1], r1[:, 0])
        theta2 = np.arctan2(r2[:, 1], r2[:, 0])

        theta[theta < starting_angle] += 2 * np.pi
        theta1[theta1 < starting_angle] += 2 * np.pi
        theta2[theta2 < starting_angle] += 2 * np.pi

        r_itp = scipy.interpolate.interp1d(theta, r_mag, kind='cubic', fill_value="extrapolate")
        r1_itp = scipy.interpolate.interp1d(theta1, r1_mag, kind='cubic', fill_value="extrapolate")
        r2_itp = scipy.interpolate.interp1d(theta2, r2_mag, kind='cubic', fill_value="extrapolate")

        pr_itp = scipy.interpolate.interp1d(theta, pr_mag, kind='cubic', fill_value="extrapolate")
        pr1_itp = scipy.interpolate.interp1d(theta1, pr1_mag, kind='cubic', fill_value="extrapolate")
        pr2_itp = scipy.interpolate.interp1d(theta2, pr2_mag, kind='cubic', fill_value="extrapolate")

        sampling_angles = np.linspace(starting_angle, theta_range + starting_angle - 1E-9, 1000)
        # sampling_angles = np.linspace(np.pi/4, 2*np.pi + np.pi/4 - 1E-9, 1000)

        r_sampled = r_itp(sampling_angles)
        r1_sampled = r1_itp(sampling_angles)
        r2_sampled = r2_itp(sampling_angles)

        pr_sampled = pr_itp(sampling_angles)
        pr1_sampled = pr1_itp(sampling_angles)
        pr2_sampled = pr2_itp(sampling_angles)

        x1 = r1_sampled - r_sampled
        x2 = r2_sampled - r_sampled
        px1 = pr1_sampled - pr_sampled
        px2 = pr2_sampled - pr_sampled

        # plt.plot(sampling_angles, x1)
        # plt.plot(sampling_angles, x2)
        # plt.show()

        tmr = np.array([[cx * x1[-1], cx * x2[-1]],
                        [cpx * px1[-1], cpx * px2[-1]]])

        _tmr = np.array([[cx * x1, cx * x2],
                         [cpx * px1, cpx * px2]])

        tm_det = tmr[0, 0] * tmr[1, 1] - tmr[0, 1] * tmr[1, 0]

        ep1 = (r_sampled[-1] - np.sqrt(r_init[0] ** 2 + r_init[1] ** 2)) * cx
        ep2 = (pr_sampled[-1] - initial_pr) * cpx
        # Positive ep1: Final radius > Initial radius
        # Negative ep1: Final radius < Initial radius
        # Positive ep2: Final radial momentum > Initial radial momentum
        # Negative ep2: Final radial momentum < Initial radial momentum

        dri = (((tmr[1, 1] - 1) * ep1 - tmr[0, 1] * ep2) / (tmr[0, 0] + tmr[1, 1] - 1 - tm_det))
        dpri = ((-ep2 - dri * tmr[1, 0]) / (tmr[1, 1] - 1))

        err = (1.0 / (initial_r * cx)) * np.sqrt(dri ** 2 + (dpri ** 2))

        if tune_calc:
            sign_r = np.sign(x2)
            nr = np.sum(((np.roll(sign_r, 1) - sign_r) != 0).astype(int))

        print(" Total Error: {}".format(err))
        print("*** Epsilon 1: {}".format(ep1))
        print("*** Epsilon 2: {}".format(ep2))

        it += 1

    print("### Final SEO Parameters ###")
    print("Energy (MeV/n): {}".format(energy_mev))
    print("# of iterations: {}".format(it))
    print("Final error: {:.6f}".format(err))
    print("Intial r: {:.6f} m".format(initial_r))
    # print("Initial pr: {:.6f} kg*m/s\n".format(initial_pr))  # TODO: Fix this output

    if tune_calc:
        z_init = np.array([initial_r, 1e-6, 0.0])
        pz_init = np.array([initial_pr, np.sqrt(total_momentum ** 2 - initial_pr ** 2), 0.0])
        vz_init = pz_init / (ion.gamma() * ion.mass_kg())

        z, vz = cr_track(ion=ion,
                         r_init=z_init,
                         v_init=vz_init,
                         symmetry="full",
                         end_type="termination",
                         input_bfield=cr._params_analytic["bf_itp"],
                         maxsteps=5000,
                         dt=1e-10,
                         omit_e=True,
                         omit_b=False)

        z_mag = np.sqrt(z[:, 0] ** 2 + z[:, 1] ** 2)
        pz = ion.gamma() * ion.mass_kg() * vz

        pz_mag = []
        for i in range(len(z_mag)):
            rhat = z[i, :2] / np.linalg.norm(z[i, :2])
        _p = np.dot(pz[i, :2], rhat)
        pz_mag.append(_p)
        pz_mag = np.array(pz_mag)

        z1_init = np.array([initial_r, 1e-6, 0.001])
        pz1_init = np.array([initial_pr, np.sqrt(total_momentum ** 2 - initial_pr ** 2), 0.0])
        vz1_init = pz1_init / (ion.gamma() * ion.mass_kg())

        z1, vz1 = cr_track(ion=ion,
                           r_init=z1_init,
                           v_init=vz1_init,
                           symmetry="full",
                           end_type="termination",
                           input_bfield=cr._params_analytic["bf_itp"],
                           maxsteps=5000,
                           dt=1e-10,
                           omit_e=True,
                           omit_b=False)

        z1_mag = np.sqrt(z1[:, 0] ** 2 + z1[:, 1] ** 2)
        pz1 = ion.gamma() * ion.mass_kg() * vz1

        pz1_mag = []
        for i in range(len(z1_mag)):
            rhat = z1[i, :2] / np.linalg.norm(z1[i, :2])
        _p = np.dot(pz1[i, :2], rhat)
        pz1_mag.append(_p)
        pz1_mag = np.array(pz1_mag)

        z2_init = np.array([initial_r, 1e-6, 0.0])
        pz_offset = initial_pr + ion.mass_kg() * clight / a
        pz2_init = np.array([initial_pr,
                             np.sqrt(total_momentum ** 2 - initial_pr ** 2 - pz_offset ** 2),
                             pz_offset])
        vz2_init = pz2_init / (ion.gamma() * ion.mass_kg())

        z2, vz2 = cr_track(ion=ion,
                           r_init=z2_init,
                           v_init=vz2_init,
                           symmetry="full",
                           end_type="termination",
                           input_bfield=cr._params_analytic["bf_itp"],
                           maxsteps=5000,
                           dt=1e-10,
                           omit_e=True,
                           omit_b=False)

        z2_mag = np.sqrt(z2[:, 0] ** 2 + z2[:, 1] ** 2)
        pz2 = ion.gamma() * ion.mass_kg() * vz2

        pz2_mag = []
        for i in range(len(z2_mag)):
            rhat = z2[i, :2] / np.linalg.norm(z2[i, :2])
        _p = np.dot(pz2[i, :2], rhat)
        pz2_mag.append(_p)
        pz2_mag = np.array(pz2_mag)

        theta = np.arctan2(z[:, 1], z[:, 0])
        theta1 = np.arctan2(z1[:, 1], z1[:, 0])
        theta2 = np.arctan2(z2[:, 1], z2[:, 0])

        theta[theta < 0] += 2 * np.pi
        theta1[theta1 < 0] += 2 * np.pi
        theta2[theta2 < 0] += 2 * np.pi

        z_itp = scipy.interpolate.interp1d(theta, z_mag, kind='cubic', fill_value="extrapolate")
        z1_itp = scipy.interpolate.interp1d(theta1, z1_mag, kind='cubic', fill_value="extrapolate")
        z2_itp = scipy.interpolate.interp1d(theta2, z2_mag, kind='cubic', fill_value="extrapolate")

        pz_itp = scipy.interpolate.interp1d(theta, pz_mag, kind='cubic', fill_value="extrapolate")
        pz1_itp = scipy.interpolate.interp1d(theta1, pz1_mag, kind='cubic', fill_value="extrapolate")
        pz2_itp = scipy.interpolate.interp1d(theta2, pz2_mag, kind='cubic', fill_value="extrapolate")

        sampling_angles = np.linspace(0.0, theta_range - 1E-9, 1000)

        z_sampled = z_itp(sampling_angles)
        z1_sampled = z1_itp(sampling_angles)
        z2_sampled = z2_itp(sampling_angles)

        pz_sampled = pz_itp(sampling_angles)
        pz1_sampled = pz1_itp(sampling_angles)
        pz2_sampled = pz2_itp(sampling_angles)

        z1 = z1_sampled - z_sampled
        z2 = z2_sampled - z_sampled
        pz1 = pz1_sampled - pz_sampled
        pz2 = pz2_sampled - pz_sampled

        tmz = np.array([[cx * z1[-1], cx * z2[-1]],
                        [cpx * pz1[-1], cpx * pz2[-1]]])

        # plt.plot(x2)
        # plt.show()
        # Tune calculations
        # cos_sigmap = 0.5 * (tm[0, 0] + tm[1, 1])
        # sigma = np.arccos(cos_sigmap) * (-1)**nr + nr * np.pi
        # nu_r = sigma / (2*np.pi)

        cos_sigma = 0.5 * (tmr[0, 0] + tmr[1, 1])
        sigma = np.arccos(cos_sigma)
        nu_r = sigma / (np.pi / 2.0)

        cos_sigma = 0.5 * (tmz[0, 0] + tmz[1, 1])
        sigma = np.arccos(cos_sigma)
        nu_z = sigma / (2.0 * np.pi)

        # TODO: Check these values with another program -PW

        print("nu_r: {}".format(nu_r))
        print("nu_z: {}".format(nu_z))

        return nu_r, nu_z

    if retrack:
        r, vr = cr_track(ion=ion,
                         r_init=r_init,
                         v_init=vr_init,
                         symmetry="full",
                         end_type="termination",
                         input_bfield=cr._params_analytic["bf_itp"],
                         maxsteps=100000,
                         dt=5e-11,
                         omit_e=True,
                         omit_b=False)

        return r, vr

    return r_init, vr_init


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
    _vx = 0.5 * ion_vel * (np.sin((2*analytic_vars["k"] + 1) * analytic_vars["b"]) +
                           np.sin(-(2*analytic_vars["k"] - 1) * analytic_vars["b"]))

    _y = -0.5 * h * (np.sin(cp * analytic_vars["b"]) / cp +
                     np.sin(-cm * analytic_vars["b"]) / cm)
    _vy = 0.5 * ion_vel * (np.cos((2*analytic_vars["k"] + 1) * analytic_vars["b"]) -
                           np.cos(-(2*analytic_vars["k"] - 1) * analytic_vars["b"]))

    _z = - h * (1.0 - np.sin(analytic_vars["b"]))
    _vz = ion_vel * np.sqrt(1 - np.sin(analytic_vars["b"])**2)

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

    # If there is a known shift, apply it now...
    # TODO: Commented this out due to possible shifting error -PW
    # if si._variables_track["shift"] is not None:
    #     analytic_vars["trj_design"] += si._variables_track["shift"]

    # TODO: This is a work in progress
    # if si._variables_optimization["x_rot"] is not None and si._params_exp["y_opt"]:
    #     xrot = np.array([[1.0, 0.0, 0.0],
    #                     [0.0, np.cos(si._variables_optimization["x_rot"]),
    #                      -np.sin(si._variables_optimization["x_rot"])],
    #                     [0.0, np.sin(si._variables_optimization["x_rot"]),
    #                      np.cos(si._variables_optimization["x_rot"])]])
    #     for i in range(analytic_params["ns"]):
    #         analytic_vars["trj_design"][i, :] = np.matmul(xrot, analytic_vars["trj_design"][i, :])

    if analytic_params["rotation"] != 0.0:
        for i in range(analytic_params["ns"]):
            analytic_vars["trj_design"][i, :] = np.matmul(analytic_vars["rot"],
                                                          analytic_vars["trj_design"][i, :])

    print("Done!")

    if si.debug:
        print("Design Trajectory:")
        print(analytic_vars["trj_design"])
        print("")

    si.analytic_parameters = analytic_params
    si.analytic_variables = analytic_vars

    return analytic_vars["trj_design"]


def generate_numerical_trajectory(si, bf=None, nsteps=100000, dt=1e-12):
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
    step = int(np.floor((ns / 50)))
    interval = [j for j in range(i_init, i_final, step)]
    interval.append(i_final)

    b = _b[interval]
    r = _r[interval, :]
    v = _v[interval, :]

    # Redefine the analytical variables (will change name eventually)
    analytic_vars["trj_design"] = r
    analytic_vars["trj_vel"] = v
    analytic_vars["b"] = b

    # If there is a known shift, apply it now...
    # TODO: Commented this out due to possible shifting error -PW
    # if si._variables_track["shift"] is not None:
    #     analytic_vars["trj_design"] += si._variables_track["shift"]

    analytic_params["ns"] = len(r[:, 0])

    if analytic_params["rotation"] != 0.0:
        for i in range(analytic_params["ns"]):
            analytic_vars["trj_design"][i, :] = np.matmul(analytic_vars["rot"],
                                                          analytic_vars["trj_design"][i, :])
            analytic_vars["trj_vel"][i, :] = np.matmul(analytic_vars["rot"],
                                                       analytic_vars["trj_vel"][i, :])

    si.analytic_parameters = analytic_params
    si.analytic_variables = analytic_vars

    return r, v


def calculate_orbit_center(k, kp, height):
    xc = height * ((1 - 2 * k * np.sin(k * np.pi)) / (1 - 4 * k ** 2) - np.sin(k * np.pi) / (2 * k - kp))
    yc = height * (2 * k / (1 - 4 * k ** 2) + 1 / (2 * k - kp)) * np.cos(k * np.pi)

    return xc, yc


def central_region_simple_track(cr):
    pass
