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

# cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))
cr.load_bfield(bfield='/home/philip/Downloads/RFQ-DIP_TestCyclotron_MainField.table',
               bf_scale=1E-4,
               spatial_unit="cm",
               extents=np.array([[-30.0, 30.0], [-30.0, 30.0], [-2.0, 2.0]]),
               extents_dims=["X", "Y", "Z"])

gap = 0.06
thickness = 0.025
cl = 0.05

my_dee = AbstractDee(opening_angle=cr.dee_opening_angle,  # Opening angle of the dee
                     r_init=0.05,  # Initial r, where the dee electrode starts
                     char_len=cl,  # Characteristic length of dee segments
                     gap=gap,  # Gap between dee and dummy dee
                     thickness=thickness)  # Thickness of the dee electrodes
my_dee.initialize()
my_dee.rotate(45, angle_unit="deg")
my_second_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                            r_init=0.05,
                            char_len=cl,
                            gap=gap,
                            thickness=thickness)
my_second_dee.initialize()
my_second_dee.rotate(90 + 45, angle_unit="deg")

my_third_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                           r_init=0.05,
                           char_len=cl,
                           gap=gap,
                           thickness=thickness)
my_third_dee.initialize(bottom_angle_offset=0, angle_unit="deg")
my_third_dee.rotate(180 + 45, angle_unit="deg")

my_fourth_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                            r_init=0.05,
                            char_len=cl,
                            gap=gap,
                            thickness=thickness)
my_fourth_dee.initialize()
my_fourth_dee.rotate(270 + 45, angle_unit="deg")

dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]

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

cr.make_dees(dees, n=4, voltage=cr.dee_voltage, gap=gap, thickness=thickness)
cr.make_dummy_dees(gap=gap, thickness=thickness)

cr.solve_bempp()
cr.calculate_potential(limits=((-0.25, 0.25), (-0.25, 0.25), (-0.01, 0.01)),
                       res=0.001,
                       domain_decomp=(4, 4, 4),
                       overlap=0)

cr.calculate_efield()

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

res = modulated_track(cr=cr,
                      r_start=r_init,
                      v_start=v_init,
                      nsteps=maxsteps * 1,
                      dt=1e-11,
                      freq=cr.rf_freq,
                      phase=np.deg2rad(180 - 35.0 / 2.0 - 60))

trj = res[0]
ax.plot(trj[:, 0], trj[:, 1], color=colors[0])

res2 = simple_tracker(cr,
                      r_start=r_init,
                      v_start=v_init,
                      nturns=1,
                      dt=1e-11,
                      phase=np.deg2rad(180 - 35.0 / 2.0 - 60))

trj = res2[0]
ax.plot(trj[:, 0], trj[:, 1], color=colors[0])

for dee in dees:
    dee.plot_segments(show=False, ax=ax)

plt.plot()