"""Hand-off position and direction of the mean momentum, every transmitted particle, x-y plane (Baseline frame = machine frame in x, y).
alpha = angle of the mean momentum from the local tangent (azimuthal direction of circulation), positive = outward (v . r_hat > 0).
Also the circulation sense (viewed from +z of the frame the npz is in) and the per-particle alpha spread.

    python -m spyral_inflector.workflows.handoff_direction <bunch npz> [...] [--mass-kg <kg>]
"""
import argparse

import numpy as np

from spyral_inflector import machine


def main(argv=None):
    """Hand-off position and mean-momentum direction (command line: see the module docstring)."""
    p = argparse.ArgumentParser()
    p.add_argument("npz", nargs="+")
    p.add_argument("--mass-kg", type=float, default=None,
                   help="ion mass [kg] for the kinetic energy (default: the machine species' mass from PyPATools)")
    a = p.parse_args(argv)
    if a.mass_kg is None:
        from PyPATools.species import IonSpecies
        a.mass_kg = IonSpecies(machine.species()).mass_kg
    M, QE = a.mass_kg, 1.602176634e-19
    for path in a.npz:
        d = np.load(path, allow_pickle=True)
        out = d["crossed"] & ~d["post_exit_hit"]
        s = d["handoff_state"][out]
        s = s[np.all(np.isfinite(s), axis=1)]
        p = s[:, 0:2]
        v = s[:, 3:5]
        pm, vm = p.mean(axis=0), v.mean(axis=0)
        r = np.linalg.norm(pm)
        rh = pm / r
        vh = vm / np.linalg.norm(vm)
        cross = rh[0] * vh[1] - rh[1] * vh[0]          # + : counter-clockwise circulation seen from +z
        sense = np.sign(cross)
        th = sense * np.array([-rh[1], rh[0]])          # tangent along the circulation
        alpha = np.degrees(np.arctan2(vh @ rh, vh @ th))
        # per particle, each against its own radial direction
        rp = np.linalg.norm(p, axis=1)
        rhp = p / rp[:, None]
        thp = sense * np.column_stack([-rhp[:, 1], rhp[:, 0]])
        vhp = v / np.linalg.norm(v, axis=1)[:, None]
        ap = np.degrees(np.arctan2(np.sum(vhp * rhp, axis=1), np.sum(vhp * thp, axis=1)))
        vz = s[:, 5]
        T = 0.5 * M * np.sum(s[:, 3:6] ** 2, axis=1) / QE
        print("{}\n  {:,d} transmitted; mean position r {:.1f} mm, azimuth {:+.2f} deg, z {:+.2f} mm; circulation {} seen from +z".format(
            path, len(s), 1e3 * r, np.degrees(np.arctan2(pm[1], pm[0])), 1e3 * s[:, 2].mean(), "counter-clockwise" if sense > 0 else "clockwise"))
        print("  mean momentum: direction {:+.2f} deg (x-y), alpha {:+.2f} deg from the tangent (+ = outward); per particle alpha mean {:+.2f}, "
              "rms {:.2f} deg; mean vz/v {:+.4f}; T {:.2f} keV".format(np.degrees(np.arctan2(vh[1], vh[0])), alpha, ap.mean(), ap.std(),
                                                                 vz.mean() / np.linalg.norm(s[:, 3:6], axis=1).mean(), 1e-3 * T.mean()))
        print("  outward = {} rotation of the momentum seen from +z".format("clockwise" if sense > 0 else "counter-clockwise"))


if __name__ == "__main__":
    main()
