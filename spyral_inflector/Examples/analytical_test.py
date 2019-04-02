from spyral_inflector import *

h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.07 / h2p.a())

si = SpiralInflector(ion=h2p,
                     method="analytical",
                     volt=12000,
                     gap=18e-3,
                     tilt=27.0,
                     dx=10e-3,
                     sigma=1.5E-3,
                     ns=60,
                     debug=False,
                     rotation=45.0)

si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

si.initialize()

si.generate_geometry()

si.set_parameter(key="h", value=0.005)

si.generate_meshed_model()

si.draw_geometry(show=True, filename='auto')
si.export_electrode_geometry(fname='electrode_macro.ivb')
si.save_geo_files()
