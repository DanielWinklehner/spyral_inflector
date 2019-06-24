from spyral_inflector import *
import numpy as np

h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.07 / h2p.a())
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

cr = CentralRegion(r_cr=0.5)
cr.initialize(xi=np.array([0.05, 0.0, 0.0]),
              vi=np.array([0.0, h2p.v_m_per_s(), 0.0]))

my_dee = Dee()
my_dee.initialize()

my_second_dee = Dee()
my_second_dee.initialize()
my_second_dee.rotate(90, angle_unit="deg")

my_third_dee = Dee()
my_third_dee.initialize()
my_third_dee.rotate(180, angle_unit="deg")

my_fourth_dee = Dee()
my_fourth_dee.initialize()
my_fourth_dee.rotate(270, angle_unit="deg")

cr.add_dee(my_dee)
cr.add_dee(my_second_dee)
cr.add_dee(my_third_dee)
cr.add_dee(my_fourth_dee)
cr.plot_dees()

thetas = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]

for theta in thetas:
    for dee in cr._dees:
        dee.next_top_segment(angle_offset=theta)
        dee.next_bottom_segment(angle_offset=theta)
cr.plot_dees()

# cr.load_bfield(bfield='/home/philip/Downloads/RFQ-DIP_TestCyclotron_MainField.table',
#                bf_scale=1E-4,
#                spatial_unit="cm",
#                extents=np.array([[-30.0, 30.0], [-30.0, 30.0], [-2.0, 2.0]]),
#                extents_dims=["X", "Y", "Z"])

# cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

