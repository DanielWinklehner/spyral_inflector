"""IsoDAR HCHC-60 spiral inflector with entrance aperture, housing and a quadrupole
doublet, optimized for the fringe fields, then exported as STEP files.

Follows Scripts/GenerateInflectorWithQuadrupoles.py of the deck. Comments marked
CHANGED say where and why this differs.
"""
import csv
import os
import re

import numpy as np
from PyPATools.beam import ParticleDistribution
from PyPATools.field import Field
from PyPATools.species import IonSpecies

from spyral_inflector import *  # noqa: F401,F403

# CHANGED: the two recursive electrode finders and the OCC-based export are gone.
# The electrodes are simply si.numerical_variables["objects"].electrodes, and
# PyElectrode.export() now goes through gmsh with the electrode's rotation and
# translation applied, so the STEP file matches the BEM mesh. The old export wrote
# the untransformed OCC solid, which dropped the axial shift the optimizer applies.

# CHANGED: BFIELD_FILE = None uses a uniform field. The deck's map is exported from
# the magnet model with the bore at negative z and Bz negative at the median plane
# (Scripts/convert_bfield.py); a map with the bore at +z has to be rotated 180 deg
# about x first (Scripts/flip_bfield_z.py), not mirrored.
BFIELD_FILE = None   # e.g. r"...\Fields\HCHC-60_CentralBField_z-40to5cm_1mm.pickle"
OUT_DIR = "hchc60_steps"

##- Injection parameters -##
MEAN_INJECTION_ENERGY = 0.069337  # MeV, mean of the RFQ particle file (was 0.07 nominal)
ION_SPECIES = "H2_1+"

##- Spiral Inflector Parameters -##
METHOD = "numerical"
SOLVER = "bempp"
VOLTAGE = 12000.0      # +V on the anode, -V on the cathode [V]; the optimizer rescales it
PLATE_GAP = 0.019      # Gap between electrodes [m]
TILT = 31.0            # k' tilt of the exit [deg]
THICKNESS = 0.01       # Electrode thickness [m]
V_SHAPE = 0.0022       # V-shape depth [m]
V_TYPE = "parabolic"
SIM_POINTS = 100       # Points along the design trajectory
ASPECT_RATIO = 2.4     # Electrode width / gap
ROTATION = 0.0         # Rotation [deg] of the inflector about z
GAMMA = 5.0            # Exit wedge cut angle [deg]
ANGLING = 11.0         # Plate face angling [deg]
DEBUG = False

H2P = ParticleDistribution(species=IonSpecies(ION_SPECIES))
H2P.set_mean_energy_z_mev(MEAN_INJECTION_ENERGY)

si = SpiralInflector(ion=H2P, method=METHOD, solver=SOLVER, volt=VOLTAGE, gap=PLATE_GAP,
                     tilt=TILT, dx=THICKNESS, sigma=V_SHAPE, vee_shape=V_TYPE, ns=SIM_POINTS,
                     aspect_ratio=ASPECT_RATIO, rotation=ROTATION, debug=DEBUG,
                     gammaAng=GAMMA, anglingAng=ANGLING)

if BFIELD_FILE is None:
    si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -0.901}))
else:
    si.load_bfield(bfield=BFIELD_FILE)

si.initialize()

si.set_parameter(key="h", value=0.005)  # Mesh characteristic length [m]
si.set_parameter(key="make_aperture", value=True)
si.set_parameter(key="aperture_params", value={"thickness": 4e-3,
                                               "radius": 50e-3,
                                               "length": 40e-3,
                                               "width": 15e-3,
                                               "top_distance": 5e-3,
                                               "bottom_distance": 10e-3,
                                               "hole_type": "rectangle",
                                               "voltage": 0.0})
si.set_parameter(key="make_housing", value=True)
si.set_parameter(key="housing_params", value={"zmin": -0.12,
                                              "zmax": 0.03,
                                              "span": True,
                                              "gap": 6e-3,
                                              "thickness": 4e-3,
                                              "voltage": 0.0,
                                              "experimental": True})

# Quadrupoles: a, b hyperbola vertex radii, radius the outer radius, aper_rad the hole
# DIAMETER of the grounded apertures (0.025 -> 12.5 mm holes, just inside the 13 mm
# electrode tips). D0,D1 (x axis) are at +V and D2,D3 (y axis) at -V for quad 1,
# D4..D7 likewise for quad 2.
# CHANGED: voltages [3500, -3500] instead of [3500, 3500]. Equal signs make both quads
# focus the same plane; the doublet (opposite signs) transmitted much more of the RFQ
# bunch. The quads move with the inflector's axial shift (move_quadrupoles=True).
si.set_parameter(key="make_quadrupoles", value=True)
si.set_parameter(key="quadrupole_params", value={"a": 0.013,
                                                 "b": 0.013,
                                                 "radius": 0.04,
                                                 "z_starts": [-0.27, -0.19],
                                                 "lengths": [0.045, 0.045],
                                                 "voltages": [3500, -3500],
                                                 "aper_rad": 0.025})

si.generate_geometry()

# CHANGED: no generate_meshed_model()/solve() here, and optimize_fringe() is replaced.
# optimize_trajectory() meshes and solves itself, for each candidate geometry. It
# sets entrance and exit truncation, the axial shift dz and the electrode voltage so
# that a test particle tracked from the bore through the whole assembly leaves flat,
# centred between the electrodes and on the median plane. The quadrupoles are left
# out while it iterates (they do not affect an on-axis particle) and put back for
# the final solve. The voltage typically comes out ~3 % below VOLTAGE: the V-shaped,
# angled faces enhance the field between them.
# (The erratic 8 to 240 s solves of the old script came from the entrance aperture
# being placed at the median plane by an absolute shift; fixed in geometry.py.)
result = si.optimize_trajectory(maxiter=15, res=0.005, exclude_quadrupoles=True)
print("entrance/exit truncation {:.3f} / {:.3f} deg, dz {:+.2f} mm, voltage {:.1f} V".format(
    result["db_entrance"], result["db_exit"], 1e3 * result["dz"], result["voltage"]))
print("final residuals: exit angle {:+.3f} deg, centering {:+.3f} mm, z {:+.3f} mm, width {:+.3f} mm".format(
    result["residual_final"][0], *(1e3 * np.asarray(result["residual_final"][1:]))))

# Potential and field on a grid for tracking. CHANGED: unequal domain decompositions
# are allowed now; the box must cover the quadrupoles.
si.calculate_potential(res=0.005,
                       limits=((-0.10, 0.10), (-0.10, 0.10), (-0.29, 0.06)),
                       domain_decomp=(3, 3, 3))
si.calculate_efield()

# One particle from the bunch start plane through quads and inflector.
r, v = si.fast_track(r_start=np.array([0.0, 0.0, -0.275]),
                     v_start=np.array([0.0, 0.0, H2P.v_mean_m_per_s]), nsteps=2400, dt=1e-10)
i_exit = int(np.argmax(np.hypot(r[:, 0], r[:, 1]) > 0.06))
print("test particle reaches r = 60 mm at z = {:+.2f} mm, vertical angle {:+.2f} deg".format(
    1e3 * r[i_exit, 2], np.degrees(np.arcsin(v[i_exit, 2] / np.linalg.norm(v[i_exit])))))

# STEP export, one file per electrode (transformations included) plus the voltages,
# which the STEP files do not carry. CHANGED: the merged assembly file goes next to
# the folder, not into it -- scripts that load every .step in a folder as one
# electrode would otherwise see every surface twice.
os.makedirs(OUT_DIR, exist_ok=True)
assembly = si.numerical_variables["objects"]
with open(os.path.join(OUT_DIR, "voltages.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Electrode", "Voltage (V)"])
    for i, e in enumerate(assembly.electrodes.values()):
        name = re.sub(r"[^A-Za-z0-9_-]+", "_", e.name).strip("_") or "electrode_{:03d}".format(i)
        fn = os.path.join(OUT_DIR, "{:03d}_{}.step".format(i, name))
        if e.export(fn) != 0:
            raise RuntimeError("STEP export failed for {}".format(e.name))
        w.writerow([name, "{:.3f}".format(e.voltage)])
        print("wrote {} ({:+.1f} V)".format(fn, e.voltage))
if assembly.export("hchc60_assembly.step") != 0:
    raise RuntimeError("assembly STEP export failed")
print("wrote hchc60_assembly.step")

si.draw_geometry(freq=50, show=False, filename="hchc60_geometry.png", aux_trajectories=[r])
