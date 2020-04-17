from spyral_inflector import *

h2p = ParticleDistribution(species=IonSpecies("H2_1+"))
h2p.set_mean_energy_z_mev(0.07)

si = SpiralInflector(ion=h2p,
                     method="numerical",
                     solver="fenics",
                     volt=12000,
                     gap=20e-3,
                     tilt=27.0,
                     aspect_ratio=2.5,
                     dx=10e-3,
                     sigma=1.5E-3,
                     ns=40,
                     debug=False)

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

si.set_parameter(key="make_cylinder", value=True)
si.set_parameter(key="cylinder_params", value={"radius": 150e-3,
                                               "zmin": -250e-3,
                                               "zmax": 80e-3,
                                               "voltage": 0.0})

si.set_parameter(key="make_housing",
                 value=True)
si.set_parameter(key="housing_params",
                 value={"zmin": -0.12,
                        "zmax": 0.03,
                        "span": True,  # zmin, zmax are ignored if span is True
                        "gap": 12.5E-3,
                        "thickness": 7.5E-3,
                        "voltage": 0.0,
                        "experimental": True})

si.generate_geometry()
si.generate_meshed_model()

si.solve()
# si.calculate_efield()

ts = time.time()

si.optimize_fringe(maxiter=3, tol=0.02, res=0.005)
print("Optimizing took {:.4f} s".format(time.time() - ts))

ts = time.time()

si.track(r_start=np.array([0.0, 0.0, -0.13]),
         v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
         nsteps=15000,
         dt=1e-11)

print("Tracking took {:.4f} s".format(time.time() - ts))

si.draw_geometry(show=True)
