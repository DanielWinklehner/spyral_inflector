from spyral_inflector import *
import numpy as np

# h2p = IonSpecies("H2_1+", 0.035)
# h2p.calculate_from_energy_mev(0.07 / h2p.a())
#
# si = SpiralInflector(ion=h2p,
#                      method="numerical",
#                      solver="bempp",
#                      volt=12000,
#                      gap=18e-3,
#                      tilt=27.0,
#                      aspect_ratio=2.5,
#                      dx=10e-3,
#                      sigma=1.5E-3,
#                      ns=60,
#                      rotation=-12.5,
#                      debug=False)
#
# # si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))
# si.load_bfield(bfield='/home/philip/Downloads/Bx_By_Bz_of_x_y_z.table', bf_scale=1E-4, spatial_unit="cm")
#
# si.initialize()
#
# si.set_parameter(key="h", value=0.005)  # Mesh characteristic length
#
# si.generate_geometry()

cr = CentralRegion()
# cr.initialize()
cr.load_bfield(bfield='/home/philip/Downloads/RFQ-DIP_TestCyclotron_MainField.table',
               bf_scale=1E-4,
               spatial_unit="cm",
               extents=np.array([[-30.0, 30.0], [-30.0, 30.0], [-2.0, 2.0]]),
               extents_dims=["X", "Y", "Z"])

bfield = cr._params_analytic["bf_itp"]

energy_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

tunes = []
xl, yl = [], []
for E in energy_list:
    x, y, t = gordon_algorithm(cr, energy_mev=E, symmetry_mode="quarter")
    xl.append(x)
    yl.append(y)
    xavg = np.mean(x)
    yavg = np.mean(y)
    tunes.append(t)

fig = plt.figure()
ax = fig.add_subplot(111)

for x, y in zip(xl, yl):
    ax.plot(x, y, color='k', linewidth=1)

ax.set_xlim([-0.3, 0.3])
ax.set_ylim([-0.3, 0.3])
ax.set_aspect(1)
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.grid(True)

plt.show()

"""

#
# h2p = IonSpecies("H2_1+", 0.035)
# h2p.calculate_from_energy_mev(0.07 / h2p.a())
#
# si = SpiralInflector(ion=h2p,
#                      method="numerical",
#                      solver="bempp",
#                      volt=12000,
#                      gap=18e-3,
#                      tilt=27.0,
#                      aspect_ratio=2.5,
#                      dx=10e-3,
#                      sigma=1.5E-3,
#                      ns=60,
#                      rotation=-12.5,
#                      debug=False)
#
# # si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))
# si.load_bfield(bfield='/home/philip/Downloads/Bx_By_Bz_of_x_y_z.table', bf_scale=1E-4, spatial_unit="cm")
#
# si.initialize()
#
# si.set_parameter(key="h", value=0.005)  # Mesh characteristic length
#
# si.generate_geometry()

cr = CentralRegion()
# cr.initialize()
cr.load_bfield(bfield='/home/philip/Downloads/RFQ-DIP_TestCyclotron_MainField.table',
               bf_scale=1E-4,
               spatial_unit="cm",
               extents=np.array([[-30.0, 30.0], [-30.0, 30.0], [-2.0, 2.0]]),
               extents_dims=["X", "Y", "Z"])

bfield = cr._params_analytic["bf_itp"]

import scipy.interpolate

sampling_angles = np.linspace(0.0, 2.0 * np.pi, 1000)
sampling_radii = np.linspace(0.0, 0.3, 1000)

new_bfield = np.zeros([1000, 1000])

for i, theta in enumerate(sampling_angles):
    for j, r in enumerate(sampling_radii):
        new_bfield[i, j] = bfield(np.array([r * np.cos(theta), r * np.sin(theta), 0.0]))[2]

# new_bfield_itp = scipy.interpolate.interp2d(sampling_angles, sampling_radii, new_bfield)
new_bfield_itp = scipy.interpolate.RegularGridInterpolator(points=[sampling_angles,
                                                                   sampling_radii],
                                                           values=new_bfield)
with open("polar_bfield.pickle", "wb") as f:
    pickle.dump(new_bfield_itp, f)

exit()

h2p = IonSpecies("H2_1+", 0.5)
h2p.calculate_from_energy_mev(1.0 / h2p.a())

gordon_algorithm(cr, energy_mev=0.5)

exit()

#
# fig = plt.figure()
# ax1 = fig.add_subplot(221)
# ax2 = fig.add_subplot(222)
# ax3 = fig.add_subplot(212)
#
# m, n = 1, 7
#
# v_angle_lims = [-1, 1]
# v_angle_lims = np.deg2rad(v_angle_lims)
#
# # v_angle = np.linspace(v_angle_lims[0], v_angle_lims[1], m)
#
# x_lims = [0.19462, 0.19463]
# x_init = np.linspace(x_lims[0], x_lims[1], n)
#
# # v_angle = np.linspace(v_angle_lims[0], v_angle_lims[1], m)
# v_angle = [0.0]
#
# rdiff = 1E20
# best_x_init = 0.0
# best_v_angle = 0.0
#
# _ts = time.time()

# iteration = 0
# while rdiff > 1E-8:
#     res = np.zeros([m, n])
#     for j in range(m):
#         for i in range(n):
#             vx_init = h2p.v_m_per_s() * np.sin(v_angle[j])
#             vy_init = h2p.v_m_per_s() * np.cos(v_angle[j])
#
#             r_init = np.array([x_init[i], 1E-9, 0.0])
#             v_init = np.array([vx_init, vy_init, 0.0])
#
#             vr_init = np.dot(r_init / np.linalg.norm(r_init), v_init)
#
#             r, v = cr_track(ion=h2p,
#                             r_init=r_init,
#                             v_init=v_init,
#                             symmetry="full",
#                             input_bfield=cr._params_analytic["bf_itp"],
#                             maxsteps=3000, dt=5e-11, omit_e=True,
#                             omit_b=False)
#
#             _r = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2)
#             # vr = np.dot(r, v)
#             vr = np.sum(r * v / np.linalg.norm(r), axis=1)
#
#             # plt.plot(_r - np.linalg.norm(r_init), vr - vr_init)
#             t = np.linspace(0, len(_r) - 1, len(_r)) * 5e-11
#             T = 1.0 / (32.8E6 / 4)
#
#             dr = np.abs(_r[-1] - np.linalg.norm(r_init))
#
#             res[j, i] = dr
#
#             # ax1.plot(r[:, 0], r[:, 1])
#             # ax2.plot(_r - np.linalg.norm(r_init), (vr - vr_init) / 1E5)
#             # ax3.plot(t / T, _r - np.linalg.norm(r_init))
#
#     jmin, imin = np.unravel_index(res.argmin(), res.shape)
#
#     best_x_init = x_init[imin]
#     best_v_angle = v_angle[jmin]
#
#     rdiff = np.min(res)
#
#     print("Iteration: {}, best x_init: {}, best rdiff: {}".format(iteration, best_x_init, rdiff))
#     iteration += 1
#
#     xr = x_lims[1] - x_lims[0]
#     x_lims = [x_init[imin]*(1 - 0.25/iteration) + 1e-8, x_init[imin]*(1 + 0.25/iteration)]
#     x_init = np.linspace(x_lims[0], x_lims[1], n)
#     print("New limits: ", x_lims)
#
# print("Orbit finder took {:.3f} s.".format(time.time() - _ts))
# print("Best values: {:.6f}, {:.6f}.".format(best_x_init, best_v_angle))
# print("Results: {:.9f} m.".format(rdiff))
#
# exit()
#
# for j in range(m):
#     for i in range(n):
#         vx_init = h2p.v_m_per_s() * np.sin(v_angle[j])
#         vy_init = h2p.v_m_per_s() * np.cos(v_angle[j])
#
#         r_init = np.array([x_init[i], 1E-9, 0.0])
#         v_init = np.array([vx_init, vy_init, 0.0])
#
#         vr_init = np.dot(r_init / np.linalg.norm(r_init), v_init)
#
#         r, v = cr_track(ion=h2p,
#                         r_init=r_init,
#                         v_init=v_init,
#                         symmetry="full",
#                         input_bfield=cr._params_analytic["bf_itp"],
#                         maxsteps=3000, dt=5e-11, omit_e=True,
#                         omit_b=False)
#
#         _r = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2)
#         # vr = np.dot(r, v)
#         vr = np.sum(r * v / np.linalg.norm(r), axis=1)
#
#         # plt.plot(_r - np.linalg.norm(r_init), vr - vr_init)
#         t = np.linspace(0, len(_r) - 1, len(_r)) * 5e-11
#         T = 1.0 / (32.8E6 / 4)
#
#         dr = np.abs(_r[-1] - np.linalg.norm(r_init))
#
#         ax1.plot(r[:, 0], r[:, 1])
#         ax2.plot(_r - np.linalg.norm(r_init), (vr - vr_init) / 1E5)
#         ax3.plot(t / T, _r - np.linalg.norm(r_init))
#
# ax1.set_xlabel('x (m)')
# ax1.set_ylabel('y (m)')
# ax1.set_xlim([-0.3, 0.3])
# ax1.set_ylim([-0.3, 0.3])
# ax1.set_aspect(1)
# ax1.grid(True)
#
# ax2.set_xlabel(r'$r - r_{init}$ (cm)')
# ax2.set_ylabel(r'$!v_{r} - v_{r, init})/10^{5}$ (m/s)')
# ax2.grid(True)
#
# ax3.set_xlabel('t / T')
# ax3.set_ylabel(r'$r - r_{init}$')
# ax3.grid(True)
#
# fig.tight_layout()
#
# plt.show()
#
# exit()

### Track equilibrium orbit

# best_v_angle = 8.31216e-09
# best_x_init = 0.194625
#
# vx_init = h2p.v_m_per_s() * np.sin(best_v_angle)
# vy_init = h2p.v_m_per_s() * np.cos(best_v_angle)
#
# r_init_eq = np.array([best_x_init, 1E-9, 0.0])
# v_init = np.array([vx_init, vy_init, 0.0])
#
# vr_init = np.dot(r_init_eq / np.linalg.norm(r_init_eq), v_init)
#
# r_eq, v_eq = cr_track(ion=h2p,
#                       r_init=r_init_eq,
#                       v_init=v_init,
#                       symmetry="full",
#                       end_type="termination",
#                       input_bfield=cr._params_analytic["bf_itp"],
#                       maxsteps=50000,
#                       dt=1e-11,
#                       omit_e=True,
#                       omit_b=False)
#
# _r_eq = np.sqrt(r_eq[:, 0] ** 2 + r_eq[:, 1] ** 2)
#
# import scipy.interpolate
# theta_eq = np.arctan2(r_eq[:, 1], r_eq[:, 0])
#
# theta_eq[theta_eq < 0] += 2*np.pi
#
# theta_eq_itp = scipy.interpolate.interp1d(theta_eq, _r_eq, kind='cubic', fill_value="extrapolate")
# sampling_angles = np.linspace(0.0, 2.0*np.pi-1E-9, 1000)
# r_eq_sampled = theta_eq_itp(sampling_angles)
#
# # vr = np.dot(r, v)
# vr_eq = np.sum(r_eq * v_eq / np.linalg.norm(r_eq), axis=1)
#
# # plt.plot(_r - np.linalg.norm(r_init), vr - vr_init)
#
# t = np.linspace(0, len(_r_eq) - 1, len(_r_eq)) * 2.5e-11
# T = 1.0 / (32.8E6 / 4)
#
# ax1.plot(r_eq[:, 0], r_eq[:, 1])
# ax2.plot(1E2 * (_r_eq - np.linalg.norm(r_init_eq)), (vr_eq - vr_init) / 1E5)
#
# ax3.plot(sampling_angles, r_eq_sampled, color='r')
#
# ####
#
# vx_init = h2p.v_m_per_s() * np.sin(best_v_angle + 0.001)
# vy_init = h2p.v_m_per_s() * np.cos(best_v_angle + 0.001)
#
# r_init = np.array([best_x_init, 1E-9, 0.0])
# v_init = np.array([vx_init, vy_init, 0.0])
#
# vr_init = np.dot(r_init / np.linalg.norm(r_init), v_init)
#
# r, v = cr_track(ion=h2p,
#                 r_init=r_init,
#                 v_init=v_init,
#                 symmetry="full",
#                 end_type="termination",
#                 input_bfield=cr._params_analytic["bf_itp"],
#                 maxsteps=50000,
#                 dt=1e-11,
#                 omit_e=True,
#                 omit_b=False)
#
# _r = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2)
#
# theta = np.arctan2(r[:, 1], r[:, 0])
#
# theta[theta < 0] += 2*np.pi
#
# theta_itp = scipy.interpolate.interp1d(theta, _r, kind='cubic', fill_value="extrapolate")
# r_sampled = theta_itp(sampling_angles)
#
# vr = np.sum(r * v / np.linalg.norm(r), axis=1)
#
# ax1.plot(r[:, 0], r[:, 1])
# ax2.plot(1E2 * (_r - np.linalg.norm(r_init)), (vr - vr_init) / 1E5)
#
# ax3.plot(sampling_angles, r_sampled, color='g')
#
# ####
#
# ax3.plot(sampling_angles, r_sampled - r_eq_sampled, 'b')
#
# ax1.set_xlabel('x (m)')
# ax1.set_ylabel('y (m)')
# ax1.set_xlim([-0.3, 0.3])
# ax1.set_ylim([-0.3, 0.3])
# ax1.set_aspect(1)
# ax1.grid(True)
#
# ax2.set_xlabel(r'$r - r_{init}$ (cm)')
# ax2.set_ylabel(r'$!v_{r} - v_{r, init})/10^{5}$ (m/s)')
# ax2.grid(True)
#
# ax3.set_xlabel('Theta (rad)')
# ax3.set_ylabel(r'$r - r_{init}$')
# ax3.grid(True)
#
# fig.tight_layout()
#
# plt.show()
#
# import scipy.signal as signal
#
# sig = (r_sampled - np.linalg.norm(r_init)) - (r_eq_sampled - np.linalg.norm(r_init_eq))
# f = np.linspace(1E-9, 15.0, 1000)
#
# pgram = signal.lombscargle(sampling_angles, sig, f, normalize=True)
#
# fig = plt.figure()
# plt.subplot(2, 1, 1)
# plt.plot(sampling_angles, sig, 'b')
# # plt.plot(sampling_angles, r_eq_sampled, 'g')
# # plt.plot(sampling_angles, r_sampled, 'k')
#
#
# plt.subplot(2, 1, 2)
# plt.plot(f, pgram, 'r')
#
# plt.show()
#
# print("Frequency at maximum: {:.4f}".format(f[np.unravel_index(pgram.argmax(), pgram.shape)]))

# dr = np.linalg.norm(r[-1, :] - r_init)
# dvy = np.linalg.norm(v[-1, 0] - v_init[0])
#
# print(r[-1, :])
# print(v[-1, :])
#
# print(r_init)
# print(v_init)
#
# print(dr)
# print(dvy)
# print(np.sqrt(dr**2 + (dvy*1E-5)**2))

# thetas = np.linspace(0, 7 * np.pi / 4, 8) + np.pi / 8.0
# r_init, r_final = 0.05, 1
#
# ra = np.array([r_init * np.cos(thetas), r_init * np.sin(thetas)])
# rb = np.array([r_final * np.cos(thetas), r_final * np.sin(thetas)])

# cr._segments.append(CRSegment(ra=ra[:, 0], rb=rb[:, 0], phase=np.pi))
# cr._segments.append(CRSegment(ra=ra[:, 1], rb=rb[:, 1], phase=0))
# cr._segments.append(CRSegment(ra=ra[:, 2], rb=rb[:, 2], phase=np.pi))
# test_seg = cr._segments[0]
# cr._segments.append(CRSegment(ra=ra[:, 3], rb=rb[:, 3], phase=0))
# cr._segments.append(CRSegment(ra=ra[:, 4], rb=rb[:, 4], phase=np.pi))
# cr._segments.append(CRSegment(ra=ra[:, 5], rb=rb[:, 5], phase=0))
# cr._segments.append(CRSegment(ra=ra[:, 6], rb=rb[:, 6], phase=np.pi))
# cr._segments.append(CRSegment(ra=ra[:, 7], rb=rb[:, 7], phase=0))


# r, v = cr.track(nsteps=25000, omit_e=True)
# r, v = cr.track_segment(end_segment=test_seg)
# # print(r)
# cr.plot_segments()
# cr.plot_bfield(100, 100)

plt.clf()

bfield = cr._params_analytic["bf_itp"]

radial_positions = np.linspace(0.0, 0.25, 126)  # 0 cm to 25 cm in 2 mm increments
# pos = np.zeros([126, 3])
# pos[:, 0] = radial_positions
B_avg = []
for k, p in enumerate(radial_positions):
    th = np.deg2rad(np.linspace(0, 2 * np.pi, 360, endpoint=False))
    _r = np.zeros([360, 3])
    _r[:, 0] = p * np.cos(th)
    _r[:, 1] = p * np.sin(th)
    B = bfield(_r)
    bavg = np.mean(B[2])
    B_avg.append(bavg)

B_avg = np.array(B_avg)
dBz = np.gradient(B_avg) * 2E3

import scipy.interpolate

B_avg_itp = scipy.interpolate.interp1d(B_avg, radial_positions, kind='cubic', fill_value="extrapolate")
dBz_itp = scipy.interpolate.interp1d(dBz, radial_positions, kind='cubic', fill_value="extrapolate")
sampling_positions = np.linspace(0.0, 0.25, 1001)

# sampling_angles = np.linspace(0.0, 2.0*np.pi-1E-9, 1000)
# r_eq_sampled = theta_eq_itp(sampling_angles)
k = (sampling_positions / B_avg_itp(sampling_positions)) * dBz_itp(sampling_positions)
# k = (radial_positions / B_avg) * dBz * 2E3

plt.plot(sampling_positions, np.sqrt(1 + k))
plt.show()

"""
