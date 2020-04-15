# from .global_variables import *
from PyPATools.field import Field
from PyPATools.pusher import ParticlePusher
import numpy as np
# import multiprocessing as mp
# import time


def track(si, r_start=None, v_start=None, nsteps=10000, dt=1e-12, omit_b=False, omit_e=False):
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


# def deflection_job(si, particle, j):
#     z_axis = Vector([0.0, 0.0, 1.0])
#     print("Starting new particle process ({})...".format(j))
#     ts = time.time()
#     r, v = si.track(r_start=particle.get_position(),
#                     v_start=particle.get_velocity(),
#                     nsteps=11000,
#                     dt=1e-11)
#     print("Particle {}: Tracking took {:.4f} s.".format(j, time.time() - ts))
#
#     trj_dir = Vector(r[-1] - r[-100])
#     deviation = 90.0 - np.rad2deg(trj_dir.angle_with(z_axis))
#
#     print("Particle {}: Deviation from xy-plane: {:.4f} degrees".format(j, deviation))
#
#
# def deflection_angle_analysis(si, bunch):
#
#     jobs = []
#
#     i = 0
#     for particle in bunch:
#         i += 1
#         p = mp.Process(target=deflection_job, args=(si, particle, i,))
#         jobs.append(p)
#         p.start()


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

    cplus = analytic_vars["c+"] = (2.0 * analytic_vars["k"] + 1.0)
    cminus = analytic_vars["c-"] = -(2.0 * analytic_vars["k"] - 1.0)

    # --- Trajectory coordinates --- #
    _x = +0.5 * h * ((2.0 / (1.0 - (4.0 * (analytic_vars["k"] ** 2.0)))) -
                     (np.cos(cplus * analytic_vars["b"]) / cplus) - np.cos(
                -cminus * analytic_vars["b"]) / cminus)

    _y = -0.5 * h * (np.sin(cplus * analytic_vars["b"]) / cplus +
                     np.sin(-cminus * analytic_vars["b"]) / cminus)

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
        eh = field_val * analytic_vars["kp"] * np.sin(_b[i + 1])
        ehx = -eh * vy / (np.sqrt(vo ** 2.0 - vz ** 2.0))
        ehy = eh * vx / (np.sqrt(vo ** 2.0 - vz ** 2.0))

        ex = field_val * vx * np.abs(vz) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0)) + ehx
        ey = field_val * vy * np.abs(vz) / (vo * np.sqrt(vo ** 2.0 - vz ** 2.0)) + ehy
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

    # If there is a known shift, apply it now...
    # TODO: Commented this out due to possible shifting error -PW
    # if si._variables_track["shift"] is not None:
    #     analytic_vars["trj_design"] += si._variables_track["shift"]

    analytic_params["ns"] = len(r[:, 0])

    if analytic_params["rotation"] != 0.0:
        for i in range(analytic_params["ns"]):
            analytic_vars["trj_design"][i, :] = np.matmul(analytic_vars["rot"],
                                                          analytic_vars["trj_design"][i, :])

    analytic_vars["trj_vel"] = v
    analytic_vars["b"] = b

    si.analytic_parameters = analytic_params
    si.analytic_variables = analytic_vars

    return r, v


def calculate_orbit_center(k, kp, height):

    xc = height * ((1 - 2 * k * np.sin(k * np.pi)) / (1 - 4 * k ** 2) - np.sin(k * np.pi) / (2 * k - kp))
    yc = height * (2 * k / (1 - 4 * k ** 2) + 1 / (2 * k - kp)) * np.cos(k * np.pi)

    return xc, yc
