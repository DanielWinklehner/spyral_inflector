from dans_pymodules import *


def track(si, r_start=None, v_start=None, nsteps=10000, dt=1e-12, omit_b=False, omit_e=False):
    # TODO: For now break if r_start or v_start are not given, later get from class properties?
    assert (r_start is not None and v_start is not None), "Have to specify r_start and v_start for now!"

    if si._variables_track["ef_itp"] is None:
        print("No E-Field has been generated. Cannot track!")
        return 1

    pusher = ParticlePusher(si._params_analytic["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    if omit_e:
        efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        efield1 = si._variables_track["ef_itp"]  # type: Field

    if omit_b:
        bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    else:
        bfield1 = si._params_analytic["bf_itp"]  # type: Field

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

    si._variables_track["trj_tracker"] = r

    return r, v


def generate_analytical_trajectory(si):
    if not si._initialized:
        si.initialize()

    print("Calculating Design Trajectory... ", end="")

    h = si._variables_analytic["height"]
    tilt = si._params_analytic["tilt"]  # type: float

    si._variables_analytic["kp"] = np.tan(np.deg2rad(tilt))  # Tilt parameter
    si._variables_analytic["k"] = ((h / si._variables_analytic["r_cyc"]) + si._variables_analytic["kp"]) / 2.0

    cp = si._variables_analytic["c+"] = (2.0 * si._variables_analytic["k"] + 1.0)
    cm = si._variables_analytic["c-"] = -(2.0 * si._variables_analytic["k"] - 1.0)

    # --- Trajectory coordinates --- #
    _x = +0.5 * h * ((2.0 / (1.0 - (4.0 * (si._variables_analytic["k"] ** 2.0)))) -
                     (np.cos(cp * si._variables_analytic["b"]) / cp) - np.cos(
                -cm * si._variables_analytic["b"]) / cm)

    _y = -0.5 * h * (np.sin(cp * si._variables_analytic["b"]) / cp +
                     np.sin(-cm * si._variables_analytic["b"]) / cm)

    _z = - h * (1.0 - np.sin(si._variables_analytic["b"]))

    si._variables_analytic["trj_design"] = np.array([_x, _y, _z]).T

    # Rotation/flip
    if not ((si._variables_analytic["bf_design"] < 0.0) ^ (si._params_analytic["ion"].q() < 0.0)):
        if si._debug:
            print("Flipping direction of cyclotron motion...", end="")
        si._variables_analytic["trj_design"][:, 1] = -si._variables_analytic["trj_design"][:, 1]

    # If there is a known shift, apply it now...
    if si._variables_track["shift"] is not None:
        si._variables_analytic["trj_design"] += si._variables_track["shift"]

    if si._params_analytic["rotation"] != 0.0:
        for i in range(si._params_analytic["ns"]):
            si._variables_analytic["trj_design"][i, :] = np.matmul(si._variables_analytic["rot"],
                                                                   si._variables_analytic["trj_design"][i, :])

    print("Done!")

    if si._debug:
        print("Design Trajectory:")
        print(si._variables_analytic["trj_design"])
        print("")

    return si._variables_analytic["trj_design"]


def generate_numerical_trajectory(si, bf=None, nsteps=100000, dt=1e-12):
    # TODO: Make sure the nsteps and dt are being consistent throughout the code
    if "nsteps" in si._params_track:
        nsteps = si._params_track["nsteps"]
    if "dt" in si._params_track:
        dt = si._params_track["dt"]

    pusher = ParticlePusher(si._params_analytic["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

    tilt = si._params_analytic["tilt"]  # type: float
    si._variables_analytic["kp"] = np.tan(np.deg2rad(tilt))  # Tilt parameter

    r_start = np.array([0.0, 0.0, -si._variables_analytic["height"]])
    v_start = np.array([0.0, 0.0, si._params_analytic["ion"].v_m_per_s()])

    _r = np.zeros([nsteps + 1, 3])
    _v = np.zeros([nsteps + 1, 3])
    _b = np.zeros([nsteps + 1])  # Store the "b" angle for geometry generation
    _r[0, :] = r_start[:]
    _v[0, :] = v_start[:]

    # Create a new electric field, which will be repeatedly re-defined
    field_val = si._variables_analytic["ef_design"]
    efield1 = Field(dim=0, field={"x": field_val, "y": 0.0, "z": 0.0})

    if bf is not None:
        bfield1 = bf
    else:
        bfield1 = si._params_analytic["bf_itp"]

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
    #     _b[i + 1] = i * dt * vo / si._variables_analytic["height"]
    #
    #     # Toprek theory with surgery
    #     Eh = field_val * si._variables_analytic["kp"] * np.sin(_b[i + 1])
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
        _b[i + 1] = i * dt * vo / si._variables_analytic["height"]

        # Toprek theory with surgery
        Eh = field_val * si._variables_analytic["kp"] * np.sin(_b[i + 1])
        Ehx = -Eh * vy / (np.sqrt(vo ** 2.0 - vz ** 2.0))
        Ehy = Eh * vx / (np.sqrt(vo ** 2.0 - vz ** 2.0))

        ex = field_val * vx * np.abs(vz) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0)) + Ehx
        ey = field_val * vy * np.abs(vz) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0)) + Ehy
        ez = -field_val * (vo ** 2.0 - vz ** 2.0) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0))
        efield1 = Field(dim=0, field={"x": ex, "y": ey, "z": ez})
        if vz < 0:  # Stop when the z-component of the velocity is zero
            if si._debug:
                print(_r[i + 1, :])  # Print the final position
            break
        i += 1
    ns = i

    try:
        i_init = np.where(_b >= si._params_analytic["b_lim"][0])[0][0]
        if i_init == 0:
            i_init = 1
    except IndexError:
        i_init = 1

    try:
        i_final = np.where(_b >= si._params_analytic["b_lim"][1])[0][0]
    except IndexError:
        i_final = i

    if si._debug:
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
    si._variables_analytic["trj_design"] = r

    # If there is a known shift, apply it now...
    if si._variables_track["shift"] is not None:
        si._variables_analytic["trj_design"] += si._variables_track["shift"]

    si._params_analytic["ns"] = len(r[:, 0])

    if si._params_analytic["rotation"] != 0.0:
        for i in range(si._params_analytic["ns"]):
            si._variables_analytic["trj_design"][i, :] = np.matmul(si._variables_analytic["rot"],
                                                                   si._variables_analytic["trj_design"][i, :])

    si._variables_analytic["trj_vel"] = v
    si._variables_analytic["b"] = b

    return r, v
