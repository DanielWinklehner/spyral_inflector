"""Deck shim: the tracking machinery now lives in spyral_inflector.tracking; this module
keeps the deck scripts' `import track_inflector as ti` interface (module-level STEP_DIR /
PARTICLES / BFIELD defaults and the argument-less load_assembly() / load_particles(n)).
The original script is in legacy/track_inflector.py."""
import os

from PyPATools.species import IonSpecies  # noqa: F401  (ti.IonSpecies)
from spyral_inflector.tracking.deck import *  # noqa: F401,F403
from spyral_inflector.tracking.deck import load_particles as _load_particles, load_step_assembly, load_bfield  # noqa: F401
from spyral_inflector.tracking.hooks import *  # noqa: F401,F403
from spyral_inflector.tracking.hooks import _PD  # noqa: F401
from spyral_inflector.tracking.handoff import (save_handoff_openpmd, save_snapshot_openpmd,  # noqa: F401
                                               PMD_SPECIES_NAMES, _pmd_species)

PROJECT = os.environ.get("SI_DECK") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
STEP_DIR = os.path.join(PROJECT, "Geometry", "final_steps")
PARTICLES = os.path.join(PROJECT, "Particles", "TenThousandRFQParticles.txt")
BFIELD = os.path.join(PROJECT, "Fields", "HCHC-60_CentralBField_zflipped.pickle")

# voltages of the original deck (only used by load_assembly(); collision checks ignore them)
VOLTAGES = {
    "SI_Anode": 12000.0, "SI_Cathode": -12000.0, "Housing": 0.0,
    "Entrance_Aperture": 0.0, "ent0": 0.0, "ext0": 0.0, "ent4": 0.0, "ext4": 0.0,
    "D0": 3500.0, "D1": 3500.0, "D2": -3500.0, "D3": -3500.0,
    "D4": 3500.0, "D5": 3500.0, "D6": -3500.0, "D7": -3500.0,
}


def load_particles(n_target, seed=20260905):
    """The RFQ file named by the module-level PARTICLES, resampled to n_target particles."""
    return _load_particles(PARTICLES, n_target, seed=seed)


def load_assembly():
    """PyElectrodeAssembly from the STEP files in the module-level STEP_DIR."""
    return load_step_assembly(STEP_DIR, VOLTAGES, name="HCHC-60 spiral inflector")
