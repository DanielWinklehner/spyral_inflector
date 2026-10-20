"""Beam and loss report of an every-particle bunch run (core + RFQ tail), CPU only (numpy).

    python beam_report.py <bunch_sc_all_level.npz> [<more runs> ...] [--labels A B ...] [--current-ma 10]

Everything is over ALL particles of the RFQ file (Daniel, 2026-10-15: the tail is real beam):
  transmission (out of the housing, post-exit hits counted as lost), and how much of it is tail;
  at the hand-off plane (electrode exit) and 55 mm on: rms vertical angle and height, centroid;
  energy spread at the hand-off (rms, and the fraction more than 3.5 % below the mean);
  bunch length as the rms RF PHASE at the hand-off, wrapped to +-180 deg around the bunch (a tail
  particle that arrives N periods late sees the same phase as one on time);
  losses by electrode and the beam power they take, split into HV (spiral electrodes, quad poles)
  and grounded surfaces. Impact energy ~ start energy - q * electrode potential (the beam's own
  potential, < 0.3 kV, is ignored). The current is spread over every particle, as bunch.py does.
Pole potentials: the run's quad voltages (json 'superpose') times the sign pattern of the basis
(geometry/voltages.csv next to the steps: poles at +-3500 V).
"""
import argparse
import csv
import json
import os

import numpy as np

M_H2P, QE = 3.3462e-27, 1.602176634e-19
p = argparse.ArgumentParser()
p.add_argument("npz", nargs="+")
p.add_argument("--labels", nargs="+", default=None)
p.add_argument("--current-ma", type=float, default=10.0)
p.add_argument("--rf-mhz", type=float, default=32.8)
a = p.parse_args()
labels = a.labels or [os.path.basename(os.path.dirname(os.path.dirname(f))) for f in a.npz]


def potentials(js, npz_path):
    """Electrode potentials [V] of the run: the spiral from the json, the poles from 'superpose'."""
    V = {k: float(v) for k, v in js.get("voltages", {}).items()}
    sp = js.get("superpose") or {}
    vs = float(js.get("vscale") or 1.0)
    for k in ("SI_Anode", "SI_Cathode"):
        if k in V:
            V[k] *= vs
    quads = sp.get("quads") or []
    if quads:
        unit = float(sp.get("unit", 3500.0))
        # the basis sign pattern: geometry/voltages.csv of the run (two levels up from the npz's folder)
        here = os.path.dirname(os.path.abspath(npz_path))
        csvs = [os.path.join(here, "..", "geometry", "voltages.csv"), os.path.join(here, "geometry", "voltages.csv")]
        pattern = {}
        for c in csvs:
            if os.path.exists(c):
                with open(c) as fh:
                    for row in csv.DictReader(fh):
                        pattern[row["Electrode"]] = float(row["Voltage (V)"])
                break
        for i, (q, _alpha) in enumerate(quads):
            for j in range(4):
                name = "D{}".format(4 * i + j)
                sign = np.sign(pattern.get(name, unit if j < 2 else -unit))
                V[name] = sign * abs(float(q)) * np.sign(float(q))
    return V


rows = []
for path, lab in zip(a.npz, labels):
    js = json.load(open(path[:-4] + ".json"))
    d = np.load(path, allow_pickle=True)
    n = len(d["crossed"])
    i_p = 1e-3 * a.current_ma / n
    tail = d["is_tail"]
    out = d["crossed"] & ~d["post_exit_hit"]
    T0 = 0.5 * M_H2P * np.sum(d["v0"] ** 2, axis=1) / QE / 1e3
    V = potentials(js, path)

    def plane(state, sel):
        ok = sel & np.all(np.isfinite(state), axis=1)
        z, v = 1e3 * state[ok, 2], state[ok, 3:6]
        zp = 1e3 * v[:, 2] / np.hypot(v[:, 0], v[:, 1])
        return ok, z, zp, v

    okh, zh, zph, vh = plane(d["handoff_state"], out)
    oka, za, zpa, _ = plane(d["asym_state"], out)
    Th = 0.5 * M_H2P * np.sum(vh ** 2, axis=1) / QE / 1e3
    th = np.asarray(d["handoff_time"], float)[okh]
    ph = 360.0 * a.rf_mhz * 1e6 * th
    ref = np.median(ph[~tail[okh]]) if (~tail[okh]).any() else np.median(ph)
    phw = (ph - ref + 180.0) % 360.0 - 180.0
    names = list(d["electrode_names"])
    hit = d["hit_electrode"]
    lost_by = {}
    for k in range(n):
        if out[k]:
            continue
        # the npz 'crossed' already excludes the housing-exit hits (crossed_plane & ~post_exit), so those have
        # to be read from post_exit_hit; they used to land in "still inside"
        nm = ("housing exit opening" if d["post_exit_hit"][k] else "after the exit" if d["crossed"][k]
              else (names[hit[k]] if hit[k] >= 0 else "still inside"))
        lost_by.setdefault(nm, []).append(k)
    hv_names = [nm for nm in lost_by if nm in ("SI_Anode", "SI_Cathode") or (nm.startswith("D") and nm[1:].isdigit())]
    def power(idx, nm):
        e_imp = np.clip(T0[idx] - V.get(nm, 0.0) / 1e3, 0.0, None)
        return float(np.sum(e_imp) * 1e3 * i_p)
    rows.append({
        "label": lab, "n": n, "T": out.mean(), "T_tail_part": (out & tail).sum() / max(out.sum(), 1), "tail_through": (out & tail).sum(),
        "ang_h": zph.std(), "z_h": zh.std(), "cent_h": (zh.mean(), zph.mean()),
        "ang_55": zpa.std(), "z_55": za.std(), "cent_55": (za.mean(), zpa.mean()),
        "dE": Th.std(), "f_low": np.mean(Th < 0.965 * Th.mean()), "dphi": phw.std(),
        "lost": {nm: len(ix) / n for nm, ix in lost_by.items()},
        "P": {nm: power(np.array(ix), nm) for nm, ix in lost_by.items()},
        "hv": hv_names, "V": V})

w = 24
print("{:44s}".format("") + "".join("{:>{w}s}".format(r["label"][:w - 1], w=w) for r in rows))
def line(lab, f):
    print("{:44s}".format(lab) + "".join("{:>{w}s}".format(f(r), w=w) for r in rows))
line("particles (core + tail)", lambda r: "{:,d}".format(r["n"]))
line("TRANSMITTED, all particles", lambda r: "{:.2f} %".format(100 * r["T"]))
line("  of which tail particles", lambda r: "{} ({:.1f} %)".format(r["tail_through"], 100 * r["T_tail_part"]))
line("rms vertical angle, hand-off", lambda r: "{:.1f} mrad".format(r["ang_h"]))
line("rms vertical angle, 55 mm on", lambda r: "{:.1f} mrad".format(r["ang_55"]))
line("rms height, hand-off / 55 mm", lambda r: "{:.2f} / {:.2f} mm".format(r["z_h"], r["z_55"]))
line("centroid at the hand-off (height / angle)", lambda r: "{:+.2f} mm / {:+.1f} mrad".format(*r["cent_h"]))
line("centroid at 55 mm (height / angle)", lambda r: "{:+.2f} mm / {:+.1f} mrad".format(*r["cent_55"]))
line("energy spread, hand-off (rms)", lambda r: "{:.2f} keV".format(r["dE"]))
line("  more than 3.5 % below the mean", lambda r: "{:.1f} %".format(100 * r["f_low"]))
line("rms RF phase at the hand-off", lambda r: "{:.1f} deg".format(r["dphi"]))
allnames = sorted({k for r in rows for k in r["lost"]}, key=lambda k: -max(r["lost"].get(k, 0) for r in rows))
print("\nLOSSES (fraction of all particles | beam power at {:g} mA)".format(a.current_ma))
for nm in allnames:
    hv = any(nm in r["hv"] for r in rows)
    line("  {} {}".format(nm, "(HV)" if hv else "(ground)"),
         lambda r, nm=nm: "{:.2f} % | {:.1f} W".format(100 * r["lost"].get(nm, 0), r["P"].get(nm, 0)) if nm in r["lost"] else "-")
line("lost on HV electrodes", lambda r: "{:.2f} % | {:.1f} W".format(100 * sum(v for k, v in r["lost"].items() if k in r["hv"]),
                                                                   sum(v for k, v in r["P"].items() if k in r["hv"])))
line("  of which spiral electrodes", lambda r: "{:.2f} % | {:.1f} W".format(100 * sum(r["lost"].get(k, 0) for k in ("SI_Anode", "SI_Cathode")),
                                                                          sum(r["P"].get(k, 0) for k in ("SI_Anode", "SI_Cathode"))))
line("lost on grounded surfaces", lambda r: "{:.2f} % | {:.1f} W".format(100 * sum(v for k, v in r["lost"].items() if k not in r["hv"]),
                                                                       sum(v for k, v in r["P"].items() if k not in r["hv"])))
print("\npotentials used [kV]:")
for r in rows:
    print("  {}: ".format(r["label"]) + ", ".join("{} {:+.2f}".format(k, v / 1e3) for k, v in r["V"].items() if abs(v) > 0))
