"""The workflow plumbing end to end on generated inputs, in a few seconds and without a BEM solve.

    python run_minimal.py --config design.toml --paths paths.toml

1. reads the design config (which includes machine.toml) with the per-machine paths file: every disk location in
   the config is written as ${work}/..., the paths file says where ${work} is on this computer;
2. puts the [machine] settings into the environment, as the chain does for every step;
3. generates the inputs the config names: a B map on a grid (axisymmetric, Maxwell-consistent to second order in r)
   and a synthetic injected beam in the 10-column text format;
4. reads the beam with spyral_inflector.tracking.deck (species, RF and start plane come from the machine settings);
5. designs the inflector analytically for the beam's mean energy in the map's central field and writes the design
   orbit, the electrode surfaces and a summary into the results folder.

The levelling chain of a real design runs the same way from its config: ``spyral-inflector chain --config ...``
(see the README); its steps solve the electrodes with BEM and track with space charge, which takes hours, not seconds.
"""
import argparse
import json
import os

import numpy as np

from spyral_inflector import machine
from spyral_inflector.workflows.config import load_config


def make_bfield(path, b0=-1.0, z0=0.06, half=0.06, step=0.006):
    """Bz = b0 / (1 + (z/z0)^2) on the axis, Br = -(r/2) dBz/dz, Bz corrected by the r^2 term: a centre-bore
    field falling off upstream (deck frame: the beam comes from negative z to the median plane at z = 0)."""
    from PyPATools.field import Field
    x = np.arange(-half, half + step / 2, step)
    y = x.copy()
    z = np.arange(-0.15, 0.03 + step / 2, step)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    g = 1.0 / (1.0 + (Z / z0) ** 2)
    dg = -2.0 * Z / z0 ** 2 * g ** 2
    d2g = (-2.0 / z0 ** 2) * g ** 2 + (-2.0 * Z / z0 ** 2) * 2.0 * g * dg
    r2 = X ** 2 + Y ** 2
    bz = b0 * (g - 0.25 * r2 * d2g)
    br_over_r = -0.5 * b0 * dg
    field = Field.from_arrays(grid={"x": x, "y": y, "z": z}, values={"x": br_over_r * X, "y": br_over_r * Y, "z": bz},
                              label="generated centre-bore field", dim=3, units="m", interpolator_backend="scipy")
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    field.save(path)
    return field


def make_beam(path, n=400, e_mev=0.070, seed=1):
    """A Gaussian bunch in the 10-column text format: x(mm) x'(mrad) y(mm) y'(mrad) z(mm) z'(mrad) phase(deg)
    time(s) energy(MeV) loss; one RF period wide in phase, 1 % rms energy spread."""
    rng = np.random.default_rng(seed)
    rows = np.zeros((n, 10))
    rows[:, 0], rows[:, 2] = rng.normal(0.0, 1.0, n), rng.normal(0.0, 1.0, n)        # mm
    rows[:, 1], rows[:, 3] = rng.normal(0.0, 10.0, n), rng.normal(0.0, 10.0, n)      # mrad
    rows[:, 6] = rng.uniform(-60.0, 60.0, n)                                         # deg
    rows[:, 8] = e_mev * (1.0 + 0.01 * rng.normal(size=n))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    np.savetxt(path, rows, header="x(mm) x'(mrad) y(mm) y'(mrad) z(mm) z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss",
               comments="")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--config", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "design.toml"))
    p.add_argument("--paths", default=None, help="per-machine paths file (default: the search of workflows.paths)")
    p.add_argument("--volt", type=float, default=10000.0, help="electrode voltage [V]")
    p.add_argument("--gap", type=float, default=0.018, help="electrode gap [m]")
    p.add_argument("--tilt", type=float, default=0.0, help="tilt k' [deg]")
    a = p.parse_args(argv)

    cfg = load_config(a.config, a.paths, required={"machine": ("species", "rf_mhz", "z_start_m"),
                                                    "inputs": ("bfield", "particles"), "run": ("results",)})
    os.environ.update(machine.to_env(cfg["machine"]))
    out = cfg["run"]["results"]
    os.makedirs(out, exist_ok=True)
    print("config {} with paths file {}; results in {}".format(cfg.source, cfg.paths.source, out))
    print("machine: {} {} at {} MHz, bunch starting at z = {} m".format(machine.name() or "(unnamed)", machine.species(),
                                                                       machine.rf_mhz(), machine.z_start()))

    from spyral_inflector.tracking.deck import load_bfield, load_particles
    bfield = make_bfield(cfg["inputs"]["bfield"])
    make_beam(cfg["inputs"]["particles"])
    r, v, ion, rows = load_particles(cfg["inputs"]["particles"], 10 ** 6)
    e_mev = float(rows[:, 8].mean())
    b_centre = float(load_bfield(cfg["inputs"]["bfield"])(np.array([[0.0, 0.0, 0.0]]))[0, 2])
    print("beam: {} particles, mean {:.2f} keV, z from {:+.4f} to {:+.4f} m; Bz at the centre {:+.3f} T".format(
        len(r), 1e3 * e_mev, r[:, 2].min(), r[:, 2].max(), b_centre))

    from PyPATools.field import Field
    from PyPATools.particles import ParticleDistribution
    from PyPATools.species import IonSpecies
    from spyral_inflector import SpiralInflector
    ion_d = ParticleDistribution(species=IonSpecies(machine.species()))
    ion_d.set_mean_energy_z_mev(e_mev)
    si = SpiralInflector(ion=ion_d, method="analytical", volt=a.volt, gap=a.gap, tilt=a.tilt, aspect_ratio=2.5,
                         dx=0.01, sigma=0.0, ns=60, debug=False)
    si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": b_centre}))
    si.initialize()
    si.generate_geometry()
    trj = np.asarray(si.analytic_variables["trj_design"])
    np.save(os.path.join(out, "design_orbit.npy"), trj)
    summary = {"machine": cfg["machine"], "n_particles": int(len(r)), "mean_energy_keV": 1e3 * e_mev, "bz_centre_T": b_centre,
               "height_m": float(si.analytic_variables["height"]), "exit_point_m": trj[-1].tolist(),
               "bunch_z_range_m": [float(r[:, 2].min()), float(r[:, 2].max())]}
    with open(os.path.join(out, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    print("design orbit: height {:.1f} mm, exit at ({:+.1f}, {:+.1f}, {:+.1f}) mm -> {}".format(
        1e3 * summary["height_m"], *(1e3 * trj[-1]), os.path.join(out, "summary.json")))
    return summary


if __name__ == "__main__":
    main()
