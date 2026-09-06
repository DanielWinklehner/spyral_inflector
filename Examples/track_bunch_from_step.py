"""Reload an exported inflector (one STEP file per electrode plus voltages.csv, as
written by hchc60_deck.py), solve the field, and track a Gaussian bunch through it
with collision detection. Reports transmission and losses per electrode.

    python track_bunch_from_step.py hchc60_steps
"""
import csv
import os
import sys

import numpy as np
from py_electrodes.py_electrodes import PyElectrode, PyElectrodeAssembly
from PyPATools.beam import ParticleDistribution
from PyPATools.field import Field
from PyPATools.pusher import Pusher
from PyPATools.species import IonSpecies
from PyPATools.trackers import Tracker, Terminator

from spyral_inflector import *  # noqa: F401,F403

STEP_DIR = sys.argv[1] if len(sys.argv) > 1 else "hchc60_steps"
BFIELD = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -0.901})   # or Field.from_file(...)
N_PARTICLES = 2000

# ---------------------------------------------------------------- geometry + voltages
volts = {}
with open(os.path.join(STEP_DIR, "voltages.csv")) as fh:
    for row in csv.DictReader(fh):
        volts[row["Electrode"]] = float(row["Voltage (V)"])

assembly = PyElectrodeAssembly("reloaded inflector")
for fn in sorted(os.listdir(STEP_DIR)):
    if not fn.lower().endswith(".step"):
        continue
    name = os.path.splitext(fn)[0].split("_", 1)[1]
    e = PyElectrode(name=name, voltage=volts.get(name, 0.0))
    if e.generate_from_file(os.path.join(STEP_DIR, fn)) != 0:
        raise RuntimeError("could not load " + fn)
    assembly.add_electrode(e)
for e in assembly.electrodes.values():
    if e.generate_mesh(brep_h=0.005) != 0:
        raise RuntimeError("failed to mesh " + e.name)
mesh = assembly.get_bempp_mesh(brep_h=0.005)
print("{} electrodes, {} triangles".format(len(assembly.electrodes), mesh["elems"].shape[1]))

# ---------------------------------------------------------------- field solve
# The SpiralInflector object only hosts the solver here; its own parameters are not used.
ion = ParticleDistribution(species=IonSpecies("H2_1+"))
ion.set_mean_energy_z_mev(0.070)
si = SpiralInflector(ion=ion, method="numerical", solver="bempp", volt=12000.0, gap=0.019,
                     tilt=31.0, dx=0.01, sigma=0.0022, vee_shape="parabolic", ns=100,
                     aspect_ratio=2.4, gammaAng=5.0, anglingAng=11.0)
si.load_bfield(bfield=BFIELD)
si.initialize()
si.set_parameter(key="h", value=0.005)
si.numerical_variables["objects"] = assembly
si.numerical_variables["full mesh"] = {"verts": mesh["verts"], "elems": mesh["elems"], "domns": mesh["domns"]}
si.solve()
si.calculate_potential(res=0.005, limits=((-0.10, 0.10), (-0.10, 0.10), (-0.29, 0.06)),
                       domain_decomp=(3, 3, 3))
si.calculate_efield()
efield = si.numerical_variables["ef_itp"]

# ---------------------------------------------------------------- bunch
# Momenta are in units of beta*gamma (8.6e-3 for 70 keV H2+): sigma_px = 1.7e-4 is a
# 20 mrad rms divergence, sigma_pz = 1e-5 about 0.2 % rms energy spread.
bunch = ParticleDistribution.generate_distribution(
    IonSpecies("H2_1+"), type=["gaussian", "gaussian", "gaussian"], s_direction="z",
    n_particles=N_PARTICLES, correlation_matrix=np.eye(6),
    sigma_x=2e-3, sigma_px=1.7e-4, sigma_y=2e-3, sigma_py=1.7e-4, sigma_z=5e-3, sigma_pz=1e-5,
    cutoff_x=3, cutoff_px=3, cutoff_y=3, cutoff_py=3, cutoff_z=3, cutoff_pz=3)
bunch.set_centroid(0.0, 0.0, -0.13)
bunch.set_mean_energy_z_mev(0.070)


class ElectrodeCollision(Terminator):
    """Kill particles whose step crosses an electrode surface, and remember where."""

    def __init__(self, assembly):
        self.assembly = assembly
        self.names = [e.name for e in assembly.electrodes.values()]
        self.hit = None

    def update(self, step, r_prev, v_prev, r, v, active, t):
        if self.hit is None:
            self.hit = np.full(len(r), -1, dtype=int)
        idx = np.where(active)[0]
        if idx.size:
            data = self.assembly.segment_intersects_surface(r_prev[idx], r[idx])
            gone = idx[data["hit_mask"]]
            if gone.size:
                self.hit[gone] = data["electrode_ids"][data["hit_mask"]]
                active = active.copy()
                active[gone] = False
        return active


class ExitRadius(Terminator):
    """Retire particles once they pass r_exit, moving outward: they have left the inflector."""

    def __init__(self, r_exit):
        self.r_exit = r_exit
        self.crossed = None

    def update(self, step, r_prev, v_prev, r, v, active, t):
        if self.crossed is None:
            self.crossed = np.zeros(len(r), dtype=bool)
        rad = np.hypot(r[:, 0], r[:, 1])
        outward = (r[:, 0] * v[:, 0] + r[:, 1] * v[:, 1]) > 0
        new = active & ~self.crossed & (rad >= self.r_exit) & outward
        if new.any():
            self.crossed[new] = True
            active = active.copy()
            active[new] = False
        return active


collision = ElectrodeCollision(assembly)
exit_plane = ExitRadius(0.070)
Tracker(Pusher(IonSpecies("H2_1+"), algorithm="rk4_rel"), efield, BFIELD,
        terminators=[exit_plane, collision]).run(bunch, 1e-10, 1600, show_progress=False)

n = N_PARTICLES
print("transmitted: {} of {} ({:.1f} %)".format(exit_plane.crossed.sum(), n, 100 * exit_plane.crossed.mean()))
for i, name in enumerate(collision.names):
    k = int((collision.hit == i).sum())
    if k:
        print("  lost on {:<20s} {:5d} ({:.1f} %)".format(name, k, 100 * k / n))
