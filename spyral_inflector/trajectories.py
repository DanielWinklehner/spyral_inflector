from dans_pymodules import *
# import multiprocessing as mp
import time


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
    ef = efield1(r[0])
    bf = bfield1(r[0])
    _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

    # Track for n steps
    for i in range(nsteps):
        ef = efield1(r[i])
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


def cr_track(ion, r_init, v_init, symmetry="full", end_type="termination", input_bfield=None,
             maxsteps=20000, dt=1e-11, omit_e=True, omit_b=False):
    from .central_region import CRSegment

    r_start = r_init
    v_start = v_init

    if symmetry == "quarter":
        segment = CRSegment(ra=np.array([0.0, 0.0]), rb=np.array([0.0, 1.0]), phase=0.0)
    elif symmetry == "half":
        segment = CRSegment(ra=np.array([0.0, 0.0]), rb=np.array([-1.0, 0.0]), phase=0.0)
    else:
        segment = CRSegment(ra=np.array([0.0, 0.0]), rb=np.array([5.0, 0.0]), phase=0.0)

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
        intersected = False
        i = 0
        while not intersected and i < maxsteps:
            ef = efield1(r[i])
            bf = bfield1(r[i])

            r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)
            intersected = segment.check_intersection(r[i, :2], r[i + 1, :2])
            # v[i + 1] = cr.dee_crossing(r[i], r[i + 1], v[i + 1], dt * i)
            i += 1

        # print("Tracked through {} steps.".format(i))
        # track_vars["trj_tracker"] = r[:i, :]
        # cr.track_variables = track_vars

        return r[:i], v[:i]

    return 0


def gordon_algorithm(cr, energy_mev=1.0):

    from dans_pymodules import MyColors
    import scipy.interpolate

    colors = MyColors()

    freq_orbit = cr.rf_frequency / cr.harmonic
    omega_orbit = 2.0 * np.pi * freq_orbit
    ion = IonSpecies("H2_1+", energy_mev)

    a = clight / omega_orbit
    b = omega_orbit / ion.q_over_m()

    total_momentum = ion.gamma() * ion.mass_kg() * ion.v_m_per_s()
    aprime = total_momentum / (ion.beta() * ion.gamma())

    dri, dpri = 0.0, 0.0
    initial_r = ion.beta() * a
    initial_pr = 0.0
    err = 1
    it = 0

    while err > 1e-3:
        print("Iteration {}:".format(it), end='')
        initial_r += dri
        # print("* Initial r = {}".format(initial_r))
        initial_pr += dpri
        # print("* Intial pr = {}".format(initial_pr))

        r_init = np.array([initial_r, 1E-11, 0.0])
        pr_init = np.array([initial_pr, np.sqrt(total_momentum ** 2 - initial_pr ** 2), 0.0])
        # print("pr_init", pr_init)
        vr_init = pr_init / (ion.gamma() * ion.mass_kg())

        # Trial orbit
        r, vr = cr_track(ion=ion,
                         r_init=r_init,
                         v_init=vr_init,
                         symmetry="quarter",
                         end_type="termination",
                         input_bfield=cr._params_analytic["bf_itp"],
                         maxsteps=50000,
                         dt=1e-11,
                         omit_e=True,
                         omit_b=False)

        r_mag = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2)
        pr = ion.gamma() * ion.mass_kg() * vr
        pr_mag = np.sum(r[:, :2] * pr[:, :2] / np.linalg.norm(r[:, :2]), axis=1)

        r1_init = r_init + np.array([0.001, 0.0, 0.0])
        initial_pr1 = initial_pr + 0.0
        pr1_init = np.array([initial_pr1, np.sqrt(total_momentum ** 2 - initial_pr1 ** 2), 0.0])
        vr1_init = pr1_init / (ion.gamma() * ion.mass_kg())

        r1, vr1 = cr_track(ion=ion,
                           r_init=r1_init,
                           v_init=vr1_init,
                           symmetry="quarter",
                           end_type="termination",
                           input_bfield=cr._params_analytic["bf_itp"],
                           maxsteps=50000,
                           dt=1e-11,
                           omit_e=True,
                           omit_b=False)

        r1_mag = np.sqrt(r1[:, 0] ** 2 + r1[:, 1] ** 2)
        pr1 = ion.gamma() * ion.mass_kg() * vr1
        pr1_mag = np.sum(r1[:, :2] * pr1[:, :2] / np.linalg.norm(r1[:, :2]), axis=1)

        r2_init = r_init
        initial_pr2 = initial_pr + total_momentum * 1E-4
        pr2_init = np.array([initial_pr2, np.sqrt(total_momentum ** 2 - initial_pr2 ** 2), 0.0])
        vr2_init = pr2_init / (ion.gamma() * ion.mass_kg())

        r2, vr2 = cr_track(ion=ion,
                           r_init=r2_init,
                           v_init=vr2_init,
                           symmetry="quarter",
                           end_type="termination",
                           input_bfield=cr._params_analytic["bf_itp"],
                           maxsteps=50000,
                           dt=1e-11,
                           omit_e=True,
                           omit_b=False)

        r2_mag = np.sqrt(r2[:, 0] ** 2 + r2[:, 1] ** 2)
        pr2 = ion.gamma() * ion.mass_kg() * vr2
        pr2_mag = np.sum(r2[:, :2] * pr2[:, :2] / np.linalg.norm(r2[:, :2]), axis=1)

        theta = np.arctan2(r[:, 1], r[:, 0])
        theta1 = np.arctan2(r1[:, 1], r1[:, 0])
        theta2 = np.arctan2(r2[:, 1], r2[:, 0])

        # theta[theta < 0] += 2 * np.pi
        # theta1[theta1 < 0] += 2 * np.pi
        # theta2[theta2 < 0] += 2 * np.pi

        # _x_itp = scipy.interpolate.interp1d(theta, r[:, 0], kind='cubic', fill_value="extrapolate")
        # _y_itp = scipy.interpolate.interp1d(theta, r[:, 1], kind='cubic', fill_value="extrapolate")

        r_itp = scipy.interpolate.interp1d(theta, r_mag, kind='cubic', fill_value="extrapolate")
        r1_itp = scipy.interpolate.interp1d(theta1, r1_mag, kind='cubic', fill_value="extrapolate")
        r2_itp = scipy.interpolate.interp1d(theta2, r2_mag, kind='cubic', fill_value="extrapolate")

        pr_itp = scipy.interpolate.interp1d(theta, pr_mag, kind='cubic', fill_value="extrapolate")
        pr1_itp = scipy.interpolate.interp1d(theta1, pr1_mag, kind='cubic', fill_value="extrapolate")
        pr2_itp = scipy.interpolate.interp1d(theta2, pr2_mag, kind='cubic', fill_value="extrapolate")

        # sampling_angles = np.linspace(0.0, 2.0 * np.pi - 1E-6, 1000)
        sampling_angles = np.linspace(0.0, 0.5 * np.pi - 1E-6, 1000)

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

        c_x = 1.0 / x1[0]
        c_px = 1.0 / px2[0]

        tm = np.array([[c_x * x1[-1], c_x * x2[-1]],
                       [c_px * px1[-1], c_px * px2[-1]]])
        tm_det = tm[0, 0]*tm[1, 1] - tm[0, 1]*tm[1, 0]

        # print(tm_det)

        ep1 = (r_sampled[-1] - np.sqrt(r_init[0] ** 2 + r_init[1] ** 2)) * c_x
        ep2 = (pr_sampled[-1] - initial_pr) * c_px
        # Positive ep1: Final radius > Initial radius
        # Negative ep1: Final radius < Initial radius
        # Positive ep2: Final radial momentum > Initial radial momentum
        # Negative ep2: Final radial momentum < Initial radial momentum

        # print("ep1")
        # print(r_sampled[-1])
        # print(np.sqrt(r_init[0] ** 2 + r_init[1] ** 2))
        # print(c_x)
        # print("ep2")
        # print(pr_sampled[-1])
        # print(initial_pr)
        # print(c_px)

        dri = (((tm[1, 1] - 1) * ep1 - tm[0, 1] * ep2) / (tm[0, 0] + tm[1, 1] - 1 - tm_det)) / c_x
        dpri = ((-ep2 - c_x * dri * tm[1, 0]) / (tm[1, 1] - 1)) / c_px

        err = (1.0 / initial_r) * np.sqrt(dri**2 + (dpri * c_px / c_x)**2)

        print(" Total Error: {}".format(err))
        print("*** Epsilon 1: {}".format(ep1))
        print("*** Epsilon 2: {}".format(ep2))
        print("*** dri/r = {}".format(dri/initial_r))
        print("*** dpri/ptot = {}\n".format(dpri/total_momentum))

        # plt.plot(r[:, 0], r[:, 1])
        # plt.plot(r1[:, 0], r1[:, 1])
        # plt.plot(r2[:, 0], r2[:, 1])
        # plt.plot(sampling_angles, r_sampled)
        # plt.show()

        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # ax.plot(r[:, 0], r[:, 1], color=colors[0])
        # ax.set_xlim([-0.3, 0.3])
        # ax.set_ylim([-0.3, 0.3])
        # ax.set_aspect(1)
        # ax.set_xlabel('x (m)')
        # ax.set_ylabel('y (m)')
        # plt.show()

        it += 1

    print("### Final SEO Parameters ###")
    print("Energy (MeV/n): {}".format(energy_mev))
    print("# of iterations: {}".format(it))
    print("Final error: {:.6f}".format(err))
    print("Intial r: {:.6f} m".format(initial_r))
    print("Initial pr: {:.6f} kg*m/s\n".format(initial_pr))  # TODO: Fix this output

    r_init = np.array([initial_r, 1e-6, 0.0])
    pr_init = np.array([initial_pr, np.sqrt(total_momentum ** 2 - initial_pr ** 2), 0.0])
    vr_init = pr_init / (ion.gamma() * ion.mass_kg())

    r, vr = cr_track(ion=ion,
                     r_init=r_init,
                     v_init=vr_init,
                     symmetry="full",
                     end_type="termination",
                     input_bfield=cr._params_analytic["bf_itp"],
                     maxsteps=40000,
                     dt=1e-10,
                     omit_e=True,
                     omit_b=False)

    theta = np.arctan2(r[:, 1], r[:, 0])
    print(r)
    # bfield = cr._params_analytic["bf_itp"]
    #
    # bz = []
    # for i in range(len(r[:, 0])):
    #     bz.append(bfield([r[i, 0], r[i, 1], r[i, 2]])[2])
    #
    # plt.plot(np.array(bz))
    # plt.show()

    vmag = np.sqrt(vr[:, 0]**2 + vr[:, 1]**2 + vr[:, 2]**2)
    plt.plot(vmag)
    plt.show()

    # theta[theta < 0] += 2 * np.pi
    # theta1[theta1 < 0] += 2 * np.pi
    # theta2[theta2 < 0] += 2 * np.pi

    sampling_angles = np.linspace(0.0, 2.0 * np.pi - 1E-6, 1000)

    _x_itp = scipy.interpolate.interp1d(theta, r[:, 0], kind='cubic', fill_value="extrapolate")
    _y_itp = scipy.interpolate.interp1d(theta, r[:, 1], kind='cubic', fill_value="extrapolate")
    x, y = _x_itp(sampling_angles), _y_itp(sampling_angles)
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(x, y, color=colors[0])
    # ax.set_xlim([-0.3, 0.3])
    # ax.set_ylim([-0.3, 0.3])
    # ax.set_aspect(1)
    # ax.set_xlabel('x (m)')
    # ax.set_ylabel('y (m)')
    # plt.show()

    # Returns best r, vr
    return _x_itp(sampling_angles), _y_itp(sampling_angles)


def generate_analytical_trajectory(si):
    if not si.initialized:
        si.initialize()

    analytic_params = si.analytic_parameters
    analytic_vars = si.analytic_variables

    print("Calculating Design Trajectory... ", end="")

    h = analytic_vars["height"]
    tilt = analytic_params["tilt"]  # type: float

    analytic_vars["kp"] = np.tan(np.deg2rad(tilt))  # Tilt parameter
    analytic_vars["k"] = ((h / analytic_vars["r_cyc"]) + analytic_vars["kp"]) / 2.0

    cp = analytic_vars["c+"] = (2.0 * analytic_vars["k"] + 1.0)
    cm = analytic_vars["c-"] = -(2.0 * analytic_vars["k"] - 1.0)

    # --- Trajectory coordinates --- #
    _x = +0.5 * h * ((2.0 / (1.0 - (4.0 * (analytic_vars["k"] ** 2.0)))) -
                     (np.cos(cp * analytic_vars["b"]) / cp) - np.cos(
                -cm * analytic_vars["b"]) / cm)

    _y = -0.5 * h * (np.sin(cp * analytic_vars["b"]) / cp +
                     np.sin(-cm * analytic_vars["b"]) / cm)

    _z = - h * (1.0 - np.sin(analytic_vars["b"]))

    analytic_vars["trj_design"] = np.array([_x, _y, _z]).T

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
    field_val = analytic_vars["ef_design"]
    efield1 = Field(dim=0, field={"x": field_val, "y": 0.0, "z": 0.0})

    if bf is not None:
        bfield1 = bf
    else:
        bfield1 = analytic_params["bf_itp"]

    # initialize the velocity half a step back:
    ef = efield1(_r[0])
    bf = bfield1(_r[0])
    _, _v[0] = pusher.push(_r[0], _v[0], ef, bf, -0.5 * dt)

    # Track for n steps
    # for i in range(nsteps):
    #     ef = efield1(_r[i])
    #     bf = bfield1(_r[i])
    #
    #     _r[i + 1], _v[i + 1] = pusher.push(_r[i], _v[i], ef, bf, dt)
    #
    #     vx, vy, vz = _v[i + 1]
    #     vo = np.sqrt(vx ** 2.0 + vy ** 2.0 + vz ** 2.0)
    #     _b[i + 1] = i * dt * vo / analytic_vars["height"]
    #
    #     # Toprek theory with surgery
    #     Eh = field_val * analytic_vars["kp"] * np.sin(_b[i + 1])
    #     Ehx = -Eh * vy / (np.sqrt(vo ** 2.0 - vz ** 2.0))
    #     Ehy = Eh * vx / (np.sqrt(vo ** 2.0 - vz ** 2.0))
    #
    #     ex = field_val * vx * np.abs(vz) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0)) + Ehx
    #     ey = field_val * vy * np.abs(vz) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0)) + Ehy
    #     ez = -field_val * (vo ** 2.0 - vz ** 2.0) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0))
    #     efield1 = Field(dim=0, field={"x": ex, "y": ey, "z": ez})
    #     if vz < 0:  # Stop when the z-component of the velocity is zero
    #         if si._debug:
    #             print(_r[i + 1, :])  # Print the final position
    #         break

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
