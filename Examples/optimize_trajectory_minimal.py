"""Minimal use of optimize_trajectory(): inflector with entrance aperture and housing in
a uniform field, the exit truncation held fixed, the other knobs solved.

Takes a few minutes (each evaluation meshes, solves and tracks).
"""
import numpy as np
from PyPATools.beam import ParticleDistribution
from PyPATools.field import Field
from PyPATools.species import IonSpecies

from spyral_inflector import *  # noqa: F401,F403

ion = ParticleDistribution(species=IonSpecies("H2_1+"))
ion.set_mean_energy_z_mev(0.070)

si = SpiralInflector(ion=ion, method="numerical", solver="bempp",
                     volt=12000.0, gap=0.019, tilt=31.0, dx=0.01, sigma=0.0022,
                     vee_shape="parabolic", ns=100, aspect_ratio=2.4,
                     gammaAng=5.0, anglingAng=11.0)
si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -0.901}))
si.initialize()
si.set_parameter(key="h", value=0.005)
si.set_parameter(key="make_aperture", value=True)
si.set_parameter(key="aperture_params", value={"thickness": 4e-3, "radius": 50e-3, "length": 40e-3,
                                               "width": 15e-3, "top_distance": 5e-3,
                                               "bottom_distance": 10e-3, "hole_type": "rectangle",
                                               "voltage": 0.0})
si.set_parameter(key="make_housing", value=True)
si.set_parameter(key="housing_params", value={"zmin": -0.12, "zmax": 0.03, "span": True, "gap": 6e-3,
                                              "thickness": 4e-3, "voltage": 0.0, "experimental": True})
si.generate_geometry()

# Knobs: 0 entrance truncation [deg], 1 exit truncation [deg], 2 axial shift dz [m],
# 3 voltage scale (solved internally when every other conductor is grounded).
# fixed={1: 1.0} holds the exit truncation at 1 deg; the rest is free.
result = si.optimize_trajectory(maxiter=12, res=0.005, solver="auto", fixed={1: 1.0},
                                bounds=((0.0, 15.0), (0.0, 15.0), (-15e-3, 15e-3), (0.9, 1.1)))

print("status:", result["status"])
print("entrance truncation {:.3f} deg, exit truncation {:.3f} deg (fixed), dz {:+.3f} mm, "
      "voltage {:.1f} V (scale {:.4f})".format(result["db_entrance"], result["db_exit"],
                                                1e3 * result["dz"], result["voltage"], result["volt_scale"]))
angle, centre, z_off, width = result["residual_final"]
print("test particle: exit angle {:+.3f} deg, centering {:+.3f} mm, z offset {:+.3f} mm, "
      "width offset {:+.3f} mm".format(angle, 1e3 * centre, 1e3 * z_off, 1e3 * width))
print("evaluations:", result["n_evaluations"])
for h in result["history"]:
    print("  {:>6s} {:2d}: knobs {} -> residual norm {:.2f}".format(
        h["tag"], h["index"], np.round(h["x"], 4), h["norm"]))
