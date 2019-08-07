from spyral_inflector import *
import numpy as np

np.random.seed(137)
#
h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.2 / h2p.a())
#
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
# # si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))
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
#
cr = CentralRegion(r_cr=[0.1, 0.3],
                   dee_voltage=70e3,
                   dee_opening_angle=42.5,
                   rf_phase=-42.5 / 2.0,
                   ion=h2p)

# cr.set_inflector(si)
cr.initialize(xi=np.array([0.05, 0.0, 0.0]),
              vi=h2p.v_m_per_s() * np.array([0.0, 1.0, 0.0]))

cr.load_bfield(bfield='/home/philip/work/C-44_AIMA_Poles_wVP_large_r_z=pm1cm.comsol',
               bf_scale=10.36,
               spatial_unit="m",
               extents=np.array([[-0.3, 0.3], [-0.3, 0.3], [-0.01, 0.01]]),
               extents_dims=["X", "Y", "Z"])

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

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlim([-3, 3])
# ax.set_ylim([-3, 3])
# ax.grid(True)
# ax.set_aspect(1)

t = np.deg2rad(np.linspace(0, 360, 360))
r = np.array([0.07 * np.sin(t), 0.07 * np.cos(t), np.ones(360)])

dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]
# dees = [my_dee]
for dee in dees:
    dee.make_transforms()

# cr.make_dees(dees, n=2, voltage=cr.dee_voltage, gap=gap, thickness=thickness)
# cr.make_dummy_dees(gap=gap, thickness=thickness)
#
# cr.plot_dees(ax=ax, show=False)
cr._abstract_dees = dees
simple_tracker(cr,
               r_start=np.array([0.05*np.cos(-np.pi / 4.0),
                                 0.05*np.sin(-np.pi / 4.0),
                                 0.0]),
               v_start=h2p.v_m_per_s() * np.array([np.cos(-np.pi / 4.0 + np.pi / 2.0),
                                                   np.sin(-np.pi / 4.0 + np.pi / 2.0),
                                                   0.0]))

# cr.show(show_screen=True)

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