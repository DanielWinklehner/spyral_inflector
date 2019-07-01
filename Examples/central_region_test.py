from spyral_inflector import *
import numpy as np

np.random.seed(137)

h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.07 / h2p.a())
#
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
                     rotation=-12.5,
                     debug=False)

si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 1.04}))
# si.load_bfield(bfield='/home/philip/Downloads/Bx_By_Bz_of_x_y_z.table', bf_scale=1E-4, spatial_unit="cm")
#
si.initialize()
#
si.set_parameter(key="h", value=0.01)  # Mesh characteristic length
#
si.generate_geometry()
si.generate_meshed_model()

fig = plt.figure()
ax = fig.add_subplot(111)

cr = CentralRegion(r_cr=0.3,
                   dee_voltage=70e3,
                   dee_opening_angle=42.5,
                   rf_phase=-42.5/2.0,
                   ion=h2p)

cr.set_inflector(si)

cr.initialize(xi=np.array([0.15, 0.0, 0.0]),
              vi=np.array([0.0, h2p.v_m_per_s(), 0.0]))

# cr.load_bfield(bfield='/home/philip/Downloads/RFQ-DIP_TestCyclotron_MainField.table',
#                bf_scale=1E-4,
#                spatial_unit="cm",
#                extents=np.array([[-30.0, 30.0], [-30.0, 30.0], [-2.0, 2.0]]),
#                extents_dims=["X", "Y", "Z"])

# cr.load_bfield(bfield='/home/philip/work/C-44_AIMA_Poles_wVP_large_r.comsol',
#                bf_scale=10.36,
#                spatial_unit="m",
#                extents=np.array([[-0.3, 0.3], [-0.3, 0.3], [-0.005, 0.005]]),
#                extents_dims=["X", "Y", "Z"])

cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))


# _r, _v = gordon_algorithm(cr,
#                           energy_mev=0.5,
#                           symmetry_mode="full",
#                           starting_angle=0.0,
#                           retrack=True,
#                           verbose=True)
# xstart = _r[0, :]
# vstart = _v[0, :]

# ax.plot(_r[:, 0], _r[:, 1], color='k', linewidth=1)

cr.make_dees(n=5, voltage=cr.dee_voltage, gap=0.056, thickness=0.025, cl=0.05)
cr.make_dummy_dees(gap=0.056, thickness=0.025)
# cr.plot_dees(ax=ax, show=False)
cr.show(show_screen=True)

# cr.solve_bempp()
# calculate_potential(cr,
#                     limits=((-0.25, 0.25), (-0.25, 0.25), (-0.005, 0.005)),
#                     res=0.0025,
#                     domain_decomp=(4, 4, 4),
#                     overlap=0)
#
# calculate_efield_bempp(cr)
#
# with open('cr_ef_itp.pickle', 'wb') as f:
#     pickle.dump(cr.numerical_variables["ef_itp"], f)

exit()


