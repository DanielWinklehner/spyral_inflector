from spyral_inflector import *

h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.07 / h2p.a())

si = SpiralInflector(ion=h2p,
                     method="numerical",
                     volt=12000,
                     gap=20e-3,
                     tilt=27.0,
                     dx=10e-3,
                     sigma=1.5E-3,
                     ns=60,
                     debug=False)

si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

si.initialize()

# si.generate_geometry()
# draw_geometry(si, freq=10, show=True)

si.set_parameter(key="h", value=0.005)  # Mesh characteristic length
si.set_parameter(key="make_aperture", value=False)
si.set_parameter(key="aperture_params", value={"thickness": 4e-3,
                                               "radius": 50e-3,
                                               "length": 45e-3,
                                               "width": 18e-3,
                                               "top_distance": 5e-3,
                                               "bottom_distance": 10e-3,
                                               "voltage": 0.0})
si.set_parameter(key="make_cylinder", value=False)
si.set_parameter(key="cylinder_params", value={"radius": 120e-3,
                                               "zmin": -150e-3,
                                               "zmax": 80e-3,
                                               "voltage": 0.0})

si.generate_meshed_model()
si.solve_bempp()

ts = time.time()
si.optimize_fringe(maxiter=2, tol=0.02, res=0.005)
print("Optimizing took {:.4f} s".format(time.time() - ts))

ts = time.time()
print("Calculating electric field...")

si.calculate_potential(res=0.002,
                       limits=((-0.08, 0.08), (-0.08, 0.08), (-0.12, 0.05)),
                       domain_decomp=(3, 3, 3))

si.calculate_efield()
print("Calculating field took {:.4f} s".format(time.time() - ts))

with open('timing.txt', 'a') as outfile:
    outfile.write("Generating electric field took {:.4f} s\n".format(time.time() - ts))

# Demonstrating saving/loading
si.save('my_inflector.pickle')

new_inflector = SpiralInflector()
new_inflector.load('my_inflector.pickle')

new_inflector.calculate_potential()
new_inflector.calculate_efield()

ts = time.time()

new_inflector.track(r_start=np.array([0.0, 0.0, -0.15]),
                    v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
                    nsteps=15000,
                    dt=1e-11)

print("Tracking took {:.4f} s".format(time.time() - ts))

with open('timing.txt', 'a') as outfile:
    outfile.write("Tracking took {:.4f} s\n".format(time.time() - ts))

# ts = time.time()
#
# fast_track(si,
#            r_start=np.array([0.0, 0.0, -0.15]),
#            v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
#            nsteps=15000,
#            dt=1e-11)
#
# print("Fast tracking took {:.4f} s".format(time.time() - ts))
#
# with open('timing.txt', 'a') as outfile:
#     outfile.write("Fast tracking took {:.4f} s\n".format(time.time() - ts))

# ts = time.time()
#
# fast_track_with_termination(si,
#                             r_start=np.array([0.0, 0.0, -0.15]),
#                             v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
#                             nsteps=15000,
#                             dt=1e-11)
#
# print("Fast tracking with termination took {:.4f} s".format(time.time() - ts))
#
# with open('timing.txt', 'a') as outfile:
#     outfile.write("Fast tracking with termination took {:.4f} s\n".format(time.time() - ts))

new_inflector.draw_geometry(show=True, filename='auto', aux_trajectories=None)
# export_electrode_geometry(si, fname='electrode_macro.ivb')
# export_aperture_geometry(si, fname='aperture_macro.ivb')
# save_geo_files(si)
