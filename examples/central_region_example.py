from spyral_inflector import *
import numpy as np

h2p = IonSpecies("H2_1+", 0.1)

cr = CentralRegion(r_cr=[0.1, 0.3],  # TODO: Deprecated
                   dee_voltage=70e3,  # Voltage applied to the dees
                   dee_opening_angle=35,
                   rf_phase=0.0,  # Initial rf phase
                   ion=h2p)

cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.0}))

gap = 0.06
thickness = 0.025
cl = 0.05

my_dee = AbstractDee(opening_angle=cr.dee_opening_angle,  # Opening angle of the dee
                     r_init=0.05,  # Initial r, where the dee electrode starts
                     char_len=cl,  # Characteristic length of dee segments
                     gap=gap,  # Gap between dee and dummy dee
                     thickness=thickness)  # Thickness of the dee electrodes
my_dee.initialize()  # Must be initialized before rotating
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

# Pack the dees in to a list
dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]

# Find initial r and v
initial_r, v_angle = orbit_finder(cr, energy_mev=0.1, verbose=True)

r_init = np.array([initial_r, 0.0, 0.0])
v_init = h2p.v_m_per_s() * np.array([np.sin(v_angle), np.cos(v_angle), 0.0])

# Make the dees and dummy dees with the defined voltages, gaps, and thicknesses
cr.make_dees(dees, n=4, voltage=cr.dee_voltage, gap=gap, thickness=thickness)
cr.make_dummy_dees(gap=gap, thickness=thickness)

# Calculate the electrostatic potential and fields
cr.solve_bempp()
cr.calculate_potential(limits=((-0.25, 0.25), (-0.25, 0.25), (-0.01, 0.01)),
                       res=0.005,
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

# Track using the full 3D field map calculated with bempp
res = modulated_track(cr=cr,
                      r_start=r_init,
                      v_start=v_init,
                      nsteps=maxsteps * 1,
                      dt=1e-11,
                      freq=cr.rf_freq,
                      phase=np.deg2rad(180 - 35.0 / 2.0 - 60))

trj = res[0]
ax.plot(trj[:, 0], trj[:, 1], color=colors[0])

# Track with the 2D approximations
res2 = simple_tracker(cr,
                      r_start=r_init,
                      v_start=v_init,
                      nturns=1,
                      dt=1e-11,
                      phase=np.deg2rad(180 - 35.0 / 2.0 - 60))

trj = res2[0]
ax.plot(trj[:, 0], trj[:, 1], '--', color=colors[1])

for dee in dees:
    dee.plot_segments(show=False, ax=ax)

plt.show()