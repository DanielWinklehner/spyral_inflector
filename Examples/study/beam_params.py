"""Transverse and longitudinal beam parameters of every-particle bunch runs, side by side, in the Baseline frame
(+z up: heights and vertical angles are the negatives of the deck-frame values).

    python beam_params.py <bunch_sc_all_level.npz> [<more> ...] --labels A B ...

All transmitted particles (core + tail; housing-exit hits excluded), as beam_report.py counts them.
Hand-off  : the plane through the design orbit's exit point, normal to the design velocity (handoff_state). Vertical z
            and z' = vz / v_horizontal; horizontal u along t = z x n (n the design exit direction) and u' = v.t / v.n;
            centroid u relative to the design orbit; rms emittances (geometric, and normalized with beta*gamma of each
            particle's speed), Twiss alpha/beta; energy and RF-phase spread.
55 mm     : asym_state, each particle's state 210 steps (21 ns) after its own crossing, ~54 mm at the design speed. NOT
            a plane: vertical quantities only (the midplane field bends the beam horizontally by ~1 rad over that
            distance, so a fixed-time horizontal snapshot mixes in the energy spread).
"""
import argparse
import json
import os
import pickle

import numpy as np

M_H2P, QE, C = 3.3462e-27, 1.602176634e-19, 299792458.0
p = argparse.ArgumentParser()
p.add_argument("npz", nargs="+")
p.add_argument("--labels", nargs="+", default=None)
p.add_argument("--rf-mhz", type=float, default=32.8)
a = p.parse_args()
labels = a.labels or [os.path.basename(f) for f in a.npz]


def emit(x, xp):
    x, xp = x - x.mean(), xp - xp.mean()
    s11, s22, s12 = np.mean(x * x), np.mean(xp * xp), np.mean(x * xp)
    e = np.sqrt(max(s11 * s22 - s12 ** 2, 0.0))
    return e, (s11 / e if e else np.nan), (-s12 / e if e else np.nan)


rows = []
for path, lab in zip(a.npz, labels):
    js = json.load(open(path[:-4] + ".json"))
    d = np.load(path, allow_pickle=True)
    out = d["crossed"] & ~d["post_exit_hit"]
    # the hand-off plane = the design exit of the basis the run tracked in (shifted like the run)
    bdir = (js.get("superpose") or {}).get("basis_dir")
    with open(os.path.join(bdir, "si_state_spiral.pickle"), "rb") as fh:
        st = pickle.load(fh)
    shift = float(js.get("shift_z_m") or 0.0)
    p0 = np.asarray(st["trj_design"], float)[-1] + np.array([0.0, 0.0, shift])
    n3 = np.asarray(st["v_design"], float)[-1]
    nh = np.array([n3[0], n3[1], 0.0]) / np.hypot(n3[0], n3[1])
    t = np.array([-nh[1], nh[0], 0.0])
    r = {"label": lab, "T": out.mean(), "R": js.get("assembly_rotation_deg"), "exit_bcs_mm": 1e3 * p0 * np.array([1, 1, -1])}
    hs = d["handoff_state"]
    ok = out & np.all(np.isfinite(hs), axis=1)
    x, v = hs[ok, :3], hs[ok, 3:]
    vh = np.hypot(v[:, 0], v[:, 1])
    z, zp = -1e3 * (x[:, 2]), -1e3 * v[:, 2] / vh               # BCS
    u, up = 1e3 * ((x - p0) @ t), 1e3 * (v @ t) / (v @ nh)
    speed = np.linalg.norm(v, axis=1)
    bg = speed / C / np.sqrt(1 - (speed / C) ** 2)
    ez, bz, az = emit(z, zp)
    eu, bu, au = emit(u, up)
    E = 0.5 * M_H2P * speed ** 2 / QE / 1e3
    th = np.asarray(d["handoff_time"], float)[ok]
    ph = 360.0 * a.rf_mhz * 1e6 * th
    tail = d["is_tail"][ok]
    ref = np.median(ph[~tail]) if (~tail).any() else np.median(ph)
    phw = (ph - ref + 180.0) % 360.0 - 180.0
    r["h"] = dict(n=int(ok.sum()), z_c=z.mean(), zp_c=zp.mean(), z_rms=z.std(), zp_rms=zp.std(),
                  zp_rob=0.5 * np.subtract(*np.percentile(zp, [84, 16])), ez=ez, ezn=ez * bg.mean(), bz=bz, az=az,
                  u_c=u.mean(), up_c=up.mean(), u_rms=u.std(), up_rms=up.std(), eu=eu, eun=eu * bg.mean(), bu=bu, au=au,
                  E=E.mean(), dE=E.std(), dphi=phw.std())
    sa = d["asym_state"]
    ok = out & np.all(np.isfinite(sa), axis=1)
    x, v = sa[ok, :3], sa[ok, 3:]
    vh = np.hypot(v[:, 0], v[:, 1])
    z, zp = -1e3 * x[:, 2], -1e3 * v[:, 2] / vh
    speed = np.linalg.norm(v, axis=1)
    bg = speed / C / np.sqrt(1 - (speed / C) ** 2)
    ez, bz, az = emit(z, zp)
    r["a"] = dict(n=int(ok.sum()), z_c=z.mean(), zp_c=zp.mean(), z_rms=z.std(), zp_rms=zp.std(),
                  zp_rob=0.5 * np.subtract(*np.percentile(zp, [84, 16])), ez=ez, ezn=ez * bg.mean(), bz=bz, az=az)
    rows.append(r)

w = 22
print("{:46s}".format("") + "".join("{:>{w}s}".format(r["label"][:w - 1], w=w) for r in rows))


def line(lab, f):
    print("{:46s}".format(lab) + "".join("{:>{w}s}".format(f(r), w=w) for r in rows))


line("transmitted (all particles)", lambda r: "{:.2f} %".format(100 * r["T"]))
line("R [deg]", lambda r: "{:+.4f}".format(r["R"]))
line("design exit point, BCS [mm]", lambda r: "{:.1f} {:.1f} {:+.2f}".format(*r["exit_bcs_mm"]))
print("HAND-OFF (electrode exit plane)")
line("  particles", lambda r: "{:,d}".format(r["h"]["n"]))
line("  vertical centroid z / z' [mm / mrad]", lambda r: "{:+.2f} / {:+.1f}".format(r["h"]["z_c"], r["h"]["zp_c"]))
line("  rms z / rms z' [mm / mrad]", lambda r: "{:.2f} / {:.1f}".format(r["h"]["z_rms"], r["h"]["zp_rms"]))
line("  z' robust (p84-p16)/2 [mrad]", lambda r: "{:.1f}".format(r["h"]["zp_rob"]))
line("  vertical emittance rms / norm [mm mrad]", lambda r: "{:.1f} / {:.3f}".format(r["h"]["ez"], r["h"]["ezn"]))
line("  vertical Twiss beta [m] / alpha", lambda r: "{:.3f} / {:+.2f}".format(r["h"]["bz"], r["h"]["az"]))
line("  horiz. centroid u / u' vs design [mm/mrad]", lambda r: "{:+.2f} / {:+.1f}".format(r["h"]["u_c"], r["h"]["up_c"]))
line("  rms u / rms u' [mm / mrad]", lambda r: "{:.2f} / {:.1f}".format(r["h"]["u_rms"], r["h"]["up_rms"]))
line("  horizontal emittance rms / norm [mm mrad]", lambda r: "{:.1f} / {:.3f}".format(r["h"]["eu"], r["h"]["eun"]))
line("  horizontal Twiss beta [m] / alpha", lambda r: "{:.3f} / {:+.2f}".format(r["h"]["bu"], r["h"]["au"]))
line("  mean energy / rms spread [keV]", lambda r: "{:.2f} / {:.2f}".format(r["h"]["E"], r["h"]["dE"]))
line("  rms RF phase [deg]", lambda r: "{:.1f}".format(r["h"]["dphi"]))
print("55 mm (21 ns after each particle's exit crossing; vertical only)")
line("  vertical centroid z / z' [mm / mrad]", lambda r: "{:+.2f} / {:+.1f}".format(r["a"]["z_c"], r["a"]["zp_c"]))
line("  rms z / rms z' [mm / mrad]", lambda r: "{:.2f} / {:.1f}".format(r["a"]["z_rms"], r["a"]["zp_rms"]))
line("  z' robust (p84-p16)/2 [mrad]", lambda r: "{:.1f}".format(r["a"]["zp_rob"]))
line("  vertical emittance rms / norm [mm mrad]", lambda r: "{:.1f} / {:.3f}".format(r["a"]["ez"], r["a"]["ezn"]))
line("  vertical Twiss beta [m] / alpha", lambda r: "{:.3f} / {:+.2f}".format(r["a"]["bz"], r["a"]["az"]))
