import pytest
from spyral_inflector import *


def test_analytical_spyral_inflector():
    h2p = IonSpecies("H2_1+", 0.035)
    h2p.calculate_from_energy_mev(0.07 / h2p.a())

    # Create spiral inflector object with the follow parameters
    si = SpiralInflector(ion=h2p,
                         method="analytical",
                         volt=12000,  # Electrode voltages
                         gap=18e-3,  # Electrode gap distance
                         tilt=0,  # Tilt angle
                         aspect_ratio=2.5,  # Aspect ratio for the electrodes
                         dx=10e-3,  # Electrode thickness
                         sigma=0,  # v-shape parameter
                         ns=60,  # Number of trajectory points
                         debug=False)

    # Load a constant magnetic field pointing in the negative z direction
    si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

    # Initialize the spiral inflector and generate the geometry
    si.initialize()

    si.generate_geometry()

    # Save the geometry as an AutoDesk Inventor macro
    si.electrode_geometry_macro(fname='electrode_macro.ivb')


def test_numerical_spyral_inflector():

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

    si.set_parameter(key="h", value=0.008)  # Mesh characteristic length
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

    ts = time.time()

    si.optimize_fringe(maxiter=1, tol=1.0, res=0.008)
    print("Optimizing took {:.4f} s".format(time.time() - ts))

    ts = time.time()

    print("Calculating electric field...")
    si.calculate_potential(res=0.008,
                           limits=((-0.08, 0.08), (-0.08, 0.08), (-0.12, 0.05)),
                           domain_decomp=(3, 3, 3))

    si.calculate_efield()
    print("Calculating field took {:.4f} s".format(time.time() - ts))

    ts = time.time()

    si.track(r_start=np.array([0.0, 0.0, -0.13]),
             v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
             nsteps=1500,
             dt=1e-10)

    print("Tracking took {:.4f} s".format(time.time() - ts))

    si.draw_geometry(show=True)
