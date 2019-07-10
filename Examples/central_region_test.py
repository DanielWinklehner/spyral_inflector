from spyral_inflector import *
import numpy as np

np.random.seed(137)

h2p = IonSpecies("H2_1+", 0.035)
# h2p.calculate_from_energy_mev(0.3 / h2p.a())
#
si = SpiralInflector(ion=h2p,
                     method="numerical",
                     solver="bempp",
                     volt=12000,
                     gap=18e-3,
                     tilt=27.0,
                     aspect_ratio=2.5,
                     dx=10e-3,
                     sigma=1.5E-3,
                     ns=60,
                     debug=False)

# si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))
si.load_bfield(bfield='/home/philip/Downloads/Bx_By_Bz_of_x_y_z.table', bf_scale=1E-4, spatial_unit="cm")

si.initialize()

si.set_parameter(key="h", value=0.01)  # Mesh characteristic length

si.generate_geometry()
si.generate_meshed_model()

r_start = si.analytic_variables["trj_design"][-1, :]
v_start = si.analytic_variables["trj_vel"][-1, :]

cr = CentralRegion(r_cr=[0.1, 0.3],
                   dee_voltage=70e3,
                   dee_opening_angle=42.5,
                   rf_phase=-42.5 / 2.0,
                   ion=h2p)

cr.set_inflector(si)

cr.initialize(xi=r_start,
              vi=v_start)

cr.track_from_si()

cr.load_bfield(bfield='/home/philip/Downloads/RFQ-DIP_TestCyclotron_MainField.table',
               bf_scale=1E-4,
               spatial_unit="cm",
               extents=np.array([[-30.0, 30.0], [-30.0, 30.0], [-2.0, 2.0]]),
               extents_dims=["X", "Y", "Z"])

trj = si.analytic_variables["trj_design"]
plt.plot(trj[:, 0], trj[:, 1])

r, v = cr.track_from_si(nsteps=10000, dt=1e-11)

plt.plot(r[:, 0], r[:, 1])
plt.show()

# r, v = cr_track(h2p,
#                 r_init=r_start,
#                 v_init=v_start,
#                 end_type="steps",
#                 dt=5e-11,
#                 maxsteps=2000,
#                 input_bfield=cr.analytic_parameters["bf_itp"])
#
# plt.plot(r[:, 0], r[:, 1])
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# energy_list = [0.125, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
#
# for E in energy_list:
#
#     h2p = IonSpecies("H2_1+", E)
#
#     initial_r, v_angle = orbit_finder(cr, E, verbose=False)
#     r_init = np.array([initial_r, 1e-12, 0.0])
#
#     v_init = np.array([h2p.v_m_per_s() * np.sin(v_angle),
#                        h2p.v_m_per_s() * np.cos(v_angle),
#                        0.0])
#
#     r, v = cr_track(h2p,
#                     r_init=r_init,
#                     v_init=v_init,
#                     dt=5e-11,
#                     nterms=1,
#                     input_bfield=cr.analytic_parameters["bf_itp"])
#
#     ax.plot(r[:, 0], r[:, 1], color='k', linewidth=1)
#
# cr.plot_bfield(nx=300, ny=300, fig=fig, ax=ax)
#
# ax.grid(True)
# ax.set_xlim([-0.3, 0.3])
# ax.set_ylim([-0.3, 0.3])
# ax.set_aspect(1)
# ax.set_xlabel('x (m)')
# ax.set_ylabel('y (m)')
#
# plt.savefig("eq_orbits.png")
# # plt.show()

gap = 0.056
thickness = 0.025
cl = 0.05

my_dee = AbstractDee(opening_angle=42.5,
                     r_init=0.05,
                     char_len=0.05,
                     gap=gap,
                     thickness=thickness)
my_dee.initialize()

my_second_dee = AbstractDee(opening_angle=42.5,
                            r_init=0.05,
                            char_len=0.05,
                            gap=gap,
                            thickness=thickness)
my_second_dee.initialize()
my_second_dee.rotate(90, angle_unit="deg")

my_third_dee = AbstractDee(opening_angle=42.5,
                           r_init=0.05,
                           char_len=0.05,
                           gap=gap,
                           thickness=thickness)
my_third_dee.initialize()
my_third_dee.rotate(180, angle_unit="deg")

my_fourth_dee = AbstractDee(opening_angle=42.5,
                            r_init=0.05,
                            char_len=0.05,
                            gap=gap,
                            thickness=thickness)
my_fourth_dee.initialize()
my_fourth_dee.rotate(270, angle_unit="deg")

dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]

cr.make_dees(dees, n=5, voltage=cr.dee_voltage, gap=gap, thickness=thickness)
cr.make_dummy_dees(gap=gap, thickness=thickness)
# cr.plot_dees(ax=ax, show=False)
cr.show(show_screen=True)

# cr.solve_bempp()
# calculate_potential(cr,
#                     limits=((-0.25, 0.25), (-0.25, 0.25), (-0.005, 0.005)),
#                     res=0.0025,
#                     domain_decomp=(4, 4, 4),
#                     overlap=0)
#
# calculate_efield_bempp(cr)
#
# with open('cr_ef_itp.pickle', 'wb') as f:
#     pickle.dump(cr.numerical_variables["ef_itp"], f)