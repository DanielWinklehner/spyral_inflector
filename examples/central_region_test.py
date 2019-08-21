from spyral_inflector import *
import numpy as np

np.random.seed(137)
#
h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.4 / h2p.a())

# si = SpiralInflector(ion=h2p,
#                      method="analytical",
#                      solver="bempp",
#                      volt=12000,
#                      gap=18e-3,
#                      tilt=20.0,
#                      aspect_ratio=2.5,
#                      dx=10e-3,
#                      sigma=1.5E-3,
#                      ns=60,
#                      debug=False)
#
# si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))
# si.load_bfield(bfield='/home/philip/Downloads/Bx_By_Bz_of_x_y_z.table', bf_scale=1E-4, spatial_unit="cm")
#
# si.initialize()
#
# si.set_parameter(key="h", value=0.01)  # Mesh characteristic length
#
# si.generate_geometry()
# si.generate_meshed_model()
# # si_exit_parameter_space(si)
#
# r_start = si.analytic_variables["trj_design"][-1, :]
# v_start = si.analytic_variables["trj_vel"][-1, :]
# print(r_start / np.linalg.norm(r_start))
# print(v_start / h2p.v_m_per_s())
# cr_start_angle = np.pi / 4.0  # Angular position of the spiral inflector
# si.apply_rotation(angle=-np.arctan2(r_start[1], r_start[0]) + cr_start_angle + 2*np.pi * (np.arctan2(r_start[1], r_start[0]) < 0), angle_unit='rad')

cr = CentralRegion(r_cr=[0.1, 0.3],
                   dee_voltage=70e3,
                   dee_opening_angle=35,
                   rf_phase=0.0,
                   ion=h2p)

# cr.set_inflector(si)
cr.initialize(xi=np.array([0.05, 0.0, 0.0]),
              vi=h2p.v_m_per_s() * np.array([0.0, 1.0, 0.0]))

# cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

# RFQ-DIP_TestCyclotron_MainField.table

# cr.load_bfield(bfield='/home/philip/work/C-44_AIMA_Poles_wVP_large_r_z=pm1cm.comsol',
#                bf_scale=10.36,
#                spatial_unit="m",
#                extents=np.array([[-0.3, 0.3], [-0.3, 0.3], [-0.01, 0.01]]),
#                extents_dims=["X", "Y", "Z"])

cr.load_bfield(bfield='/home/philip/Downloads/RFQ-DIP_TestCyclotron_MainField.table',
               bf_scale=1E-4,
               spatial_unit="cm",
               extents=np.array([[-30.0, 30.0], [-30.0, 30.0], [-2.0, 2.0]]),
               extents_dims=["X", "Y", "Z"])

# bf = cr.analytic_parameters["bf_itp"]
# xr = np.linspace(0.0, 0.3, 100)
# bl = np.zeros_like(xr)
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# for i, x in enumerate(xr):
#     bl[i] = np.abs(bf(np.array([x, 0.0, 0.0]))[2])
#
# ax.plot(xr, bl, color=colors[0])
# cr.load_bfield(bfield='/home/philip/work/C-44_AIMA_Poles_wVP_large_r_z=pm1cm.comsol',
#                bf_scale=10.36,
#                spatial_unit="m",
#                extents=np.array([[-0.3, 0.3], [-0.3, 0.3], [-0.01, 0.01]]),
#                extents_dims=["X", "Y", "Z"])
# bf = cr.analytic_parameters["bf_itp"]
#
# xr = np.linspace(0.0, 0.3, 100)
# bl = np.zeros_like(xr)
# for i, x in enumerate(xr):
#     bl[i] = np.abs(bf(np.array([x, 0.0, 0.0]))[2])
# ax.plot(xr, bl, color=colors[1])
# ax.set_xlabel('r (m)')
# ax.set_ylabel(r'$\vert B \vert$ (T)')
# ax.grid(True)
# ax.set_xlim([0, 0.3])
# ax.set_ylim([0, 1.5])
# plt.show()

gap = 0.06
thickness = 0.025
cl = 0.05

my_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                     r_init=0.05,
                     char_len=0.05,
                     gap=gap,
                     thickness=thickness)
my_dee.initialize()
my_dee.rotate(45, angle_unit="deg")
my_second_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                            r_init=0.05,
                            char_len=0.05,
                            gap=gap,
                            thickness=thickness)
my_second_dee.initialize()
my_second_dee.rotate(90 + 45, angle_unit="deg")

my_third_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                           r_init=0.05,
                           char_len=0.05,
                           gap=gap,
                           thickness=thickness)
my_third_dee.initialize(bottom_angle_offset=0, angle_unit="deg")
my_third_dee.rotate(180 + 45, angle_unit="deg")

my_fourth_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                            r_init=0.05,
                            char_len=0.05,
                            gap=gap,
                            thickness=thickness)
my_fourth_dee.initialize()
my_fourth_dee.rotate(270 + 45, angle_unit="deg")

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlim([-3, 3])
# ax.set_ylim([-3, 3])
# ax.grid(True)
# ax.set_aspect(1)

dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]
# dees = [my_dee]
cr._abstract_dees = dees

for dee in dees:
    dee.make_transforms()
    for k in range(10):
        if k == 0 and dee is my_third_dee:
            dee.next_bottom_segment(angle_offset=0, angle_unit='deg')
        else:
            dee.next_bottom_segment()
        dee.next_top_segment()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([-0.3, 0.3])
ax.set_ylim([-0.3, 0.3])
ax.set_aspect(1)
ax.grid(True)
cr.plot_dees(ax=ax, show=False)
plt.show()

orbit_h2p = IonSpecies("H2_1+", 0.1)

# initial_r, v_angle = orbit_finder(cr, energy_mev=0.1, verbose=False)

# r_init = np.array([1.17023418e-1 - 0.0085, 1e-12, 0.0])
# v_init = np.array([0.0,
#                    orbit_h2p.v_m_per_s(),
#                    0.0])
r_init = np.array([1e-12, 1.17023418e-1 - 0.0085, 0.0])
v_init = np.array([-orbit_h2p.v_m_per_s(),
                   0.0,
                   0.0])

print(r_init)
print(v_init)

# r_init = np.array([0.0, 1.75*0.0883, 0.0])
# v_init = orbit_h2p.v_m_per_s() * np.array([-1.0, 0.0, 0.0])

# theta = np.arctan2(r[:, 1], r[:, 0])
# r_init = r[theta > np.pi / 2.0][0, :]
# v_init = v[theta > np.pi / 2.0][0, :]

# fig = plt.figure()
# ax = fig.add_subplot(111)

# cr.make_dees(dees, n=4, voltage=cr.dee_voltage, gap=gap, thickness=thickness)
# cr.make_dummy_dees(gap=gap, thickness=thickness)
# cr.plot_dees(ax=ax, show=True)
# cr.solve_bempp()
# calculate_potential(cr,
#                     limits=((-0.25, 0.25), (-0.25, 0.25), (-0.01, 0.01)),
#                     res=0.001,
#                     domain_decomp=(4, 4, 4),
#                     overlap=0)
#
# calculate_efield_bempp(cr)
# print("Importing electric field... ", end="")
# with open('dee_35deg_ef_itp.pickle', 'rb') as f:
#     cr.numerical_variables["ef_itp"] = pickle.load(f)
#     # pickle.dump(cr.numerical_variables["ef_itp"], f)
# print("Done!")

omega_rf = 2.0 * np.pi * cr.rf_freq
rf_phase = cr.rf_phase
dt = 1e-11
maxsteps = int(4.0 * 2.0 * np.pi / (dt * omega_rf))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([-0.2, 0.2])
ax.set_ylim([-0.2, 0.2])
ax.set_aspect(1)
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.grid(True)

# res = modulated_track(cr=cr,
#                       r_start=r_init,
#                       v_start=v_init,
#                       nsteps=maxsteps * 1,
#                       dt=1e-11,
#                       freq=cr.rf_freq,
#                       phase=np.deg2rad(180 - 35.0 / 2.0 - 60))
# trj = res[0]
# ax.plot(trj[:, 0], trj[:, 1], color=colors[0])

res2 = simple_tracker(cr,
                      r_start=r_init,
                      v_start=v_init,
                      nturns=1,
                      dt=1e-11,
                      phase=np.deg2rad(180 - 35.0 / 2.0 - 60))

# res_list = [res, res2]

trj = res2[0]
ax.plot(trj[:, 0], trj[:, 1], color=colors[0])

for dee in dees:
    dee.plot_segments(show=False, ax=ax)

# for res in res_list:
#     trj = res[0]
#     ax.plot(trj[:, 0], trj[:, 1])



plt.show()
#
# for res in res_list:
#     r = res[0]
#     ef = cr.numerical_variables["ef_itp"]
#     _ef = np.zeros([maxsteps*4, 3])
#     for i in range(maxsteps*4):
#         pos = r[i, :]
#         _ef[i, :] = ef(pos)
#
#     # plt.plot(_ef[:, 0])
#     # plt.plot(_ef[:, 1], '--')
#     t = np.linspace(0, 4, len(_ef[:, 0]))
#     plt.plot(t, np.sqrt(_ef[:, 0]**2 + _ef[:, 1]**2))
#
# plt.show()
#
# # for res in res_list:
#     # v = res[1]
#     # vmag = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2)
#     # plt.plot(vmag)
#
# # plt.show()
#
# for res in res_list:
#     trj = res[0]
#     plt.plot(trj[:, 2])
# plt.show()