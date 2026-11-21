"""phase_coupling.py - how much of the hand-off phase length is transverse-to-longitudinal coupling of the line (RFQ exit -> electrode
exit), the (l|x), (l|x'), (l|y), (l|y') terms that Barnard et al., PRAB 24, 023501 (2021) remove with shaped field gradients
(their inflector C3: rf spread +-39 deg -> +-22 deg). CPU only.

    python phase_coupling.py <bunch_*.npz> [...] --labels A B ...

Per every-particle run, over the transmitted CORE particles (the tail's arrival is set by its RFQ phase, not by the line):
linear least squares of the hand-off RF phase on the start values (phase phi0, relative energy delta, x, x', y, y' at the RFQ exit,
start plane = SI_Z_START) and the variance split: what phi0 + delta (the RFQ and the drift, R55/R56 like) explain, what the four
transverse start coordinates (R51..R54 like) add, and the phase length the run would have with those four terms at zero.
2026-10-21: the coefficients stay a fit to the core (they are properties of the line), but the phase lengths reported first are
over EVERY transmitted particle, core and tail (beam_report.py's wrapped phase), with and without the fitted transverse part
(Daniel: no core-only numbers); the core-sample lines are labelled as the fit sample.
"""
import argparse
import os

import numpy as np


def main(argv=None):
    """Transverse-to-longitudinal coupling of the line (command line: see the module docstring)."""
    p = argparse.ArgumentParser()
    p.add_argument("npz", nargs="+")
    p.add_argument("--labels", nargs="+", default=None)
    p.add_argument("--rf-mhz", type=float, default=32.8)
    p.add_argument("--z-start", type=float, default=float(os.environ.get("SI_Z_START", "-0.270")))
    a = p.parse_args(argv)
    F = a.rf_mhz * 1e6
    M_H2P, QE = 3.3462e-27, 1.602176634e-19
    for path, lab in zip(a.npz, a.labels or a.npz):
        if not os.path.exists(path):    # e.g. a no-space-charge run, no longer made by default (2026-11-07)
            print("(skipped {}: {} does not exist)".format(lab, path))
            continue
        d = np.load(path, allow_pickle=True)
        tail = d["is_tail"]
        ok = d["crossed"] & ~d["post_exit_hit"] & np.all(np.isfinite(d["handoff_state"]), axis=1) & np.isfinite(d["handoff_time"]) & ~tail
        r0, v0 = d["r0"][ok], d["v0"][ok]
        sp = np.linalg.norm(v0, axis=1)
        # when the particle crosses the start plane: its injection time (timed injection, --inject-core / the tail) plus, for a
        # particle placed in a spatial bunch, its drift to the plane (2026-10-23: the injection time was missing, so a run with
        # timed injection had a constant start phase)
        t_start = np.asarray(d["t_inject"], float)[ok] - (r0[:, 2] - a.z_start) / v0[:, 2]
        ph0 = 360.0 * F * t_start
        th = np.asarray(d["handoff_time"], float)[ok]
        ph = 360.0 * F * th
        ph = (ph - np.median(ph) + 180.0) % 360.0 - 180.0
        T0 = 0.5 * M_H2P * sp ** 2 / QE
        delta = 1e3 * (T0 / T0.mean() - 1.0)                                # per mille
        xp, yp = v0[:, 0] / v0[:, 2], v0[:, 1] / v0[:, 2]
        dz = r0[:, 2] - a.z_start
        X = np.column_stack([1e3 * (r0[:, 0] - xp * dz), 1e3 * xp, 1e3 * (r0[:, 1] - yp * dz), 1e3 * yp])   # mm, mrad, mm, mrad
        y = ph - ph.mean()
        L = np.column_stack([ph0 - ph0.mean(), delta - delta.mean()])
        T = X - X.mean(axis=0)
        cL = np.linalg.lstsq(L, y, rcond=None)[0]
        cA = np.linalg.lstsq(np.column_stack([L, T]), y, rcond=None)[0]
        r_L = y - L @ cL
        r_A = y - np.column_stack([L, T]) @ cA
        trans_part = T @ cA[2:]
        no_coupling = y - trans_part                                        # the phase with the four transverse terms set to zero
        # EVERY transmitted particle, core and tail: beam_report.py's wrapped hand-off phase (around the core's median), and the
        # same phase with the core-fitted transverse part taken off each particle (start coordinates at the start plane: the
        # tail is injected AT the plane at its arrival time t_inject, so its x, x', y, y' are its file values)
        okA = d["crossed"] & ~d["post_exit_hit"] & np.all(np.isfinite(d["handoff_state"]), axis=1) & np.isfinite(d["handoff_time"])
        r0A, v0A = d["r0"][okA], d["v0"][okA]
        xpA, ypA = v0A[:, 0] / v0A[:, 2], v0A[:, 1] / v0A[:, 2]
        dzA = r0A[:, 2] - a.z_start
        XA = np.column_stack([1e3 * (r0A[:, 0] - xpA * dzA), 1e3 * xpA, 1e3 * (r0A[:, 1] - ypA * dzA), 1e3 * ypA])
        phA = 360.0 * F * np.asarray(d["handoff_time"], float)[okA]
        phA = (phA - np.median(phA[~tail[okA]]) + 180.0) % 360.0 - 180.0
        transA = (XA - X.mean(axis=0)) @ cA[2:]
        phA_nc = (phA - transA + 180.0) % 360.0 - 180.0
        print("{}: every transmitted particle ({:,d}, {:,d} of them tail): hand-off phase {:.1f} deg rms; with the line's (l|x), (l|x'), "
              "(l|y), (l|y') taken off: {:.1f} deg rms (the transverse part alone {:.1f} deg rms)".format(
                  lab, int(okA.sum()), int((okA & tail).sum()), phA.std(), phA_nc.std(), transA.std()))
        print("  coefficients, fitted to the transmitted core ({:,d} particles): (l|x) {:+.2f} deg/mm, (l|x') {:+.2f} deg/mrad, (l|y) {:+.2f} deg/mm, "
              "(l|y') {:+.2f} deg/mrad; phi0 {:+.3f}, delta {:+.2f} deg/permille".format(int(ok.sum()), cA[2], cA[3], cA[4], cA[5], cA[0], cA[1]))
        print("  per plane, the rms phase the fitted terms put in: x-plane {:.1f} deg, y-plane {:.1f} deg (start beam rms x {:.2f} mm / x' {:.1f} mrad, "
              "y {:.2f} mm / y' {:.1f} mrad)".format(((XA[:, :2] - X.mean(axis=0)[:2]) @ cA[2:4]).std(), ((XA[:, 2:] - X.mean(axis=0)[2:]) @ cA[4:6]).std(),
                                                      *XA.std(axis=0)))
        print("  fit sample (transmitted core only, NOT a beam number): phase {:.1f} deg rms; RFQ start phase {:.1f} deg rms, energy spread "
              "{:.1f} per mille rms".format(y.std(), ph0.std(), delta.std()))
        print("  fit sample: explained by phi0 + delta alone {:.1f} deg (residual {:.1f}); the transverse terms carry {:.1f} deg rms; adding "
              "them leaves {:.1f} deg; with them at zero {:.1f} deg".format(np.sqrt(max(y.var() - r_L.var(), 0)), r_L.std(), trans_part.std(),
                                                                            r_A.std(), no_coupling.std()))


if __name__ == "__main__":
    main()
