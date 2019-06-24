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

# cr.load_bfield(bfield='/home/philip/Downloads/RFQ-DIP_TestCyclotron_MainField.table',
#                bf_scale=1E-4,
#                spatial_unit="cm",
#                extents=np.array([[-30.0, 30.0], [-30.0, 30.0], [-2.0, 2.0]]),
#                extents_dims=["X", "Y", "Z"])

cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

