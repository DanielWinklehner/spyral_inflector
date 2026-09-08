"""Bunch tracking through an exported spiral-inflector geometry.

Pipeline: build_geometry (knobs -> optimized design particle -> STEP files + state)
-> solve_step_assembly (bempp field of the STEP assembly at the final voltages)
-> track_bunch (RFQ bunch, optional PyAMG space charge, openPMD hand-off files)
-> exit_metrics / plots; run_final chains the last three for one geometry.
"""
from .deck import (Z_START, SPECIES, CLIGHT, RF_FREQ_HZ, QUAD1, QUAD2, QUAD_SIGN, SPIRAL_NAMES, QUAD_NAMES,
                   read_voltages, write_voltages, set_quad_voltages, set_spiral_voltage, load_state,
                   load_particles, read_dst, dst_rows, particle_rows, orient_beam, mean_energy_mev,
                   load_step_assembly, rotate_quads, mesh_assembly,
                   drop_electrodes, with_fast_interpolator, load_bfield, superpose_basis)
from .hooks import (ElectrodeCollision, ExitPlane, PlaneCrossing, TrajectoryRecorder, SnapshotRecorder, Envelope,
                    SpaceCharge, plane_axes, continue_design)
from .handoff import save_handoff_openpmd, save_snapshot_openpmd
from .metrics import exit_metrics, fmt as fmt_metrics
from .bem_reload import solve_step_assembly
from .bunch import track_bunch
from .geometry_point import build_geometry, DEFAULT_KNOBS
from .run_final import run_final
from .plots import electrode_meshes, plot_geometry_trajectories, plot_side_view

__all__ = [
    "Z_START", "SPECIES", "CLIGHT", "RF_FREQ_HZ", "QUAD1", "QUAD2", "QUAD_SIGN", "SPIRAL_NAMES", "QUAD_NAMES",
    "read_voltages", "write_voltages", "set_quad_voltages", "set_spiral_voltage", "load_state",
    "load_particles", "read_dst", "dst_rows", "particle_rows", "orient_beam", "mean_energy_mev", "load_step_assembly", "rotate_quads", "mesh_assembly",
    "drop_electrodes", "with_fast_interpolator", "load_bfield", "superpose_basis",
    "ElectrodeCollision", "ExitPlane", "PlaneCrossing", "TrajectoryRecorder", "SnapshotRecorder", "Envelope",
    "SpaceCharge", "plane_axes", "continue_design", "save_handoff_openpmd", "save_snapshot_openpmd",
    "exit_metrics", "fmt_metrics", "solve_step_assembly", "track_bunch", "build_geometry", "DEFAULT_KNOBS", "run_final",
    "electrode_meshes", "plot_geometry_trajectories", "plot_side_view",
]
