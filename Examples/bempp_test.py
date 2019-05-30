from spyral_inflector import *

h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.07 / h2p.a())

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

si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

si.initialize()

# si.generate_geometry()
# draw_geometry(si, freq=10, show=True)

si.set_parameter(key="h", value=0.01)  # Mesh characteristic length
si.set_parameter(key="make_aperture", value=False)
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
                 value=False)
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

si.solve()

# ts = time.time()
#
# si.optimize_fringe(maxiter=3, tol=0.02, res=0.005)
# print("Optimizing took {:.4f} s".format(time.time() - ts))
#
ts = time.time()

print("Calculating electric field...")
si.calculate_potential(res=0.005,
                       limits=((-0.08, 0.08), (-0.08, 0.08), (-0.12, 0.05)),
                       domain_decomp=(3, 3, 3))

si.calculate_efield()
print("Calculating field took {:.4f} s".format(time.time() - ts))

ts = time.time()

si.track(r_start=np.array([0.0, 0.0, -0.13]),
         v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
         nsteps=15000,
         dt=1e-11)

print("Tracking took {:.4f} s".format(time.time() - ts))

si.draw_geometry(show=True)
