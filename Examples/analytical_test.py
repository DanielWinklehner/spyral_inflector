from spyral_inflector import *

# Define the ion species with energy 35 keV/n
h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.07 / h2p.a())

# Create spiral inflector object with the follow parameters
si = SpiralInflector(ion=h2p,
                     method="analytical",
                     volt=12000,        # Electrode voltages
                     gap=18e-3,         # Electrode gap distance
                     tilt=0,            # Tilt angle
                     aspect_ratio=2.5,  # Aspect ratio for the electrodes
                     dx=10e-3,          # Electrode thickness
                     sigma=0,           # v-shape parameter
                     ns=60,             # Number of trajectory points
                     debug=False)

# Load a constant magnetic field pointing in the negative z direction
si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

# Initialize the spiral inflector and generate the geometry
si.initialize()

si.generate_geometry()

# Save the geometry as an AutoDesk Inventor macro
si.electrode_geometry_macro(fname='electrode_macro.ivb')