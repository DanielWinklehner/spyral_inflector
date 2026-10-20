"""Particle files at a probe plane and at the start (Daniel, 2026-10-20: "add a probe plane halfway through the SI
entrance aperture and record all particles that cross the plane. Plane normal is the z axis. Recorded particles in the
format x, y, z (will be the same for all), vx, vy, vz in m and m/s, plus the particle phase in rad. Then also save the
initial particles in that same format in a separate file.")

Both files: Baseline (machine) frame, +z up, the beam travels towards -z (so vz < 0); the beam is turned about z by
the run's phi (90 deg + R) exactly as it was tracked, so the two files share one frame. Columns
    x [m]  y [m]  z [m]  vx [m/s]  vy [m/s]  vz [m/s]  phase [rad]
phase = 2 pi f_rf t, the TraceWin convention of the RFQ file: omega * time relative to the RFQ reference particle,
whose crossing of the start plane is t = 0 (positive = late). It is NOT wrapped: a tail particle that arrives N RF
periods late carries a phase N * 2 pi larger. Rows keep the RFQ file's order (unlost particles); the probe file lists
only the particles that reached the plane, and a companion *_index.txt gives each row's index in the initial file.

    python probe_files.py <bunch npz with probe_*> --out <dir> [--particles <.dst>]
"""
import argparse
import json
import os

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("npz")
p.add_argument("--out", required=True)
p.add_argument("--particles", default=None, help="default: the run json's particle file")
p.add_argument("--prefix", default=None, help="file-name prefix (default: the npz's run folder name)")
a = p.parse_args()

from spyral_inflector.tracking.deck import load_particles_with_tail, orient_beam, Z_START  # noqa: E402

js = json.load(open(a.npz[:-4] + ".json"))
d = np.load(a.npz, allow_pickle=True)
rf_hz = 1e6 * float(js["sc"]["rf_MHz"]) if js.get("sc") else 32.8e6
phi = float(js["phi_deg"])
particles = a.particles or js["particles"]
run = os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(a.npz)))) if os.path.basename(
    os.path.dirname(os.path.abspath(a.npz))) in ("probe", "sc_level") else os.path.basename(os.path.dirname(os.path.abspath(a.npz)))
prefix = a.prefix or run
os.makedirs(a.out, exist_ok=True)
MIR = np.array([1.0, 1.0, -1.0])
HEAD = ("Baseline (machine) frame: +z up, beam travels towards -z; beam turned about z by phi = {phi:+.4f} deg as tracked "
        "(run {run}, R {R:+.4f} deg).\nphase = 2 pi f t with f = {f:.4f} MHz: omega * time relative to the RFQ reference particle "
        "(TraceWin convention, t = 0 at its start-plane crossing; positive = late; not wrapped).\n{what}\n"
        "x [m]  y [m]  z [m]  vx [m/s]  vy [m/s]  vz [m/s]  phase [rad]")

# initial particles: the RFQ file at the start plane, as the run loaded and oriented them
r0, v0, ion, rows, t_inject, is_tail = load_particles_with_tail(particles, 10 ** 7, rf_mhz=rf_hz / 1e6)
r0, v0 = orient_beam(r0, v0, phi)
n = len(r0)
if n != len(d["crossed"]):
    raise SystemExit("the RFQ file gives {} particles but the run tracked {}".format(n, len(d["crossed"])))
ph0 = np.radians(rows[:, 6])
r_start = np.column_stack([r0[:, 0], r0[:, 1], np.full(n, float(Z_START))]) * MIR      # every particle AT the start plane
init = np.column_stack([r_start, v0 * MIR, ph0])
fn0 = os.path.join(a.out, "{}_initial_particles.txt".format(prefix))
np.savetxt(fn0, init, fmt="%+.9e", header=HEAD.format(phi=phi, run=run, R=float(js["assembly_rotation_deg"]), f=rf_hz / 1e6,
           what="INITIAL particles: every unlost particle of {} ({} core + {} RFQ tail) at the start plane z = {:+.4f} m "
                "(RFQ exit).".format(os.path.basename(particles), int((~is_tail).sum()), int(is_tail.sum()), -float(Z_START))))
print("wrote {} ({} particles)".format(fn0, n))

# probe planes
for k, zp in enumerate(np.atleast_1d(d["probe_z"])):
    st, tm = d["probe_state"][k], d["probe_time"][k]
    ok = np.all(np.isfinite(st), axis=1) & np.isfinite(tm)
    out = np.column_stack([st[ok, :3] * MIR, st[ok, 3:] * MIR, 2.0 * np.pi * rf_hz * tm[ok]])
    tag = "probe_z{:+.2f}mm".format(-1e3 * float(zp)).replace("+", "p").replace("-", "m").replace(".", "_")
    fn = os.path.join(a.out, "{}_{}.txt".format(prefix, tag))
    np.savetxt(fn, out, fmt="%+.9e", header=HEAD.format(phi=phi, run=run, R=float(js["assembly_rotation_deg"]), f=rf_hz / 1e6,
               what="PROBE plane z = {:+.5f} m (normal along z): every particle's first crossing, {} of {} ({} RFQ tail); "
                    "{:g} mA {}.".format(-float(zp), int(ok.sum()), n, int((ok & d["is_tail"]).sum()),
                                         float(js["sc"]["current_mA"]) if js.get("sc") else 0.0,
                                         "with space charge" if js.get("space_charge") else "no space charge")))
    np.savetxt(fn[:-4] + "_index.txt", np.where(ok)[0], fmt="%d",
               header="row index in {} of each row of {}".format(os.path.basename(fn0), os.path.basename(fn)))
    print("wrote {} ({} of {} particles crossed; phase {:.2f}..{:.2f} rad, median {:.2f})".format(
        fn, int(ok.sum()), n, out[:, 6].min(), out[:, 6].max(), np.median(out[:, 6])))
