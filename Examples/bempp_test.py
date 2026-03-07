from spyral_inflector import *

h2p = ParticleDistribution(species=IonSpecies("H2_1+"))
h2p.set_mean_energy_z_mev(0.07)

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

# si.load_bfield()
si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

si.initialize()

# si.generate_geometry()
# draw_geometry(si, freq=10, show=True)

si.set_parameter(key="h", value=0.005)  # Mesh characteristic length
si.set_parameter(key="make_aperture", value=True)
si.set_parameter(key="aperture_params", value={"thickness": 4e-3,
                                               "radius": 50e-3,
                                               "length": 45e-3,
                                               "width": 18e-3,
                                               "top_distance": 5e-3,
                                               "bottom_distance": 10e-3,
                                               "hole_type": "rectangle",
                                               "voltage": 0.0})

si.set_parameter(key="make_cylinder", value=False)
si.set_parameter(key="cylinder_params", value={"radius": 120e-3,
                                               "zmin": -250e-3,
                                               "zmax": 80e-3,
                                               "voltage": 0.0})

si.set_parameter(key="make_housing",
                 value=True)
si.set_parameter(key="housing_params",
                 value={"zmin": -0.12,
                        "zmax": 0.03,
                        "span": True,
                        "gap": 10E-3,
                        "thickness": 5E-3,
                        "voltage": 0.0,
                        "experimental": True})

si.generate_geometry()
si.generate_meshed_model()

ts = time.time()
si.solve()
print("Solving took {:.4f} s".format(time.time() - ts))

# ts = time.time()
# si.optimize_fringe(maxiter=2, tol=0.02, res=0.005)
# print("Optimizing took {:.4f} s".format(time.time() - ts))

ts = time.time()
print("Calculating electric potential/field...")
si.calculate_potential(res=0.005,
                       limits=((-0.08, 0.08), (-0.08, 0.08), (-0.12, 0.05)),
                       domain_decomp=(3, 3, 3))
si.calculate_efield()
print("Calculating field took {:.4f} s".format(time.time() - ts))

ts = time.time()
bunch = ParticleDistribution.generate_distribution(IonSpecies("H2_1+"),
                                                   type=['gaussian', 'gaussian', 'gaussian'],
                                                   s_direction='z',
                                                   n_particles=200,
                                                   correlation_matrix=np.eye(6),
                                                   sigma_x=2e-3,
                                                   sigma_px=1e-20,
                                                   sigma_y=2e-3,
                                                   sigma_py=1e-20,
                                                   sigma_z=4e-3,
                                                   sigma_pz=1e-20,
                                                   cutoff_x=3,
                                                   cutoff_px=3,
                                                   cutoff_y=3,
                                                   cutoff_py=3,
                                                   cutoff_z=3,
                                                   cutoff_pz=3
                                                   )

bunch.set_centroid(0.0, 0.0, -0.13)
bunch.set_mean_energy_z_mev(0.07)

r, v, active = si.fast_track_batch_with_termination(r_start=bunch.x_vec,
                                                    v_start=bunch.v_vec,
                                                    nsteps=1500,
                                                    dt=1e-10)
print("Tracking took {:.4f} s".format(time.time() - ts))

r = np.swapaxes(r, 0, 1)
si.draw_geometry(freq=50, show=True, aux_trajectories=r)

# print(si.electrode_geometry_macro())
# print(si.aperture_geometry_macro())
