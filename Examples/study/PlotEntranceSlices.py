"""Electrode cross-sections at several heights near the entrance (BCS frame, top view), with
the design orbit's point at that height, the axis and the cylinder bore.

    python PlotEntranceSlices.py --run Results/final/final1 [--out fig.png]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spyral_inflector.tracking.deck import load_step_assembly, mesh_assembly, load_state
from spyral_inflector.tracking.exit_plane import rotz

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
p = argparse.ArgumentParser()
p.add_argument("--run", required=True)
p.add_argument("--out", default=None)
p.add_argument("--heights", type=float, nargs="+", default=[104, 98, 90, 80, 70, 55], help="BCS z [mm]")
a = p.parse_args()
run = a.run if os.path.isabs(a.run) else os.path.join(DECK, a.run)

R = json.load(open(os.path.join(run, "exit_plane.json")))["rotation_deg"]
Rz = rotz(R)
st = load_state(os.path.join(run, "geometry", "state.pickle"))
trj = np.asarray(st["trj_design"]) @ Rz.T
assy = load_step_assembly(os.path.join(run, "steps"))
for u in [u for u, e in assy.electrodes.items() if e.name not in ("SI_Anode", "SI_Cathode", "Housing_Cylinder", "Entrance_Aperture", "Housing")]:
    assy.electrodes.pop(u)
mesh_assembly(assy)
V = {e.name: np.asarray(e._gmsh_msh["vertices"], float) @ Rz.T for e in assy.electrodes.values()}
zb_trj = -trj[:, 2]                                    # BCS heights of the orbit points
colors = {"SI_Anode": "tab:red", "SI_Cathode": "tab:blue", "Housing_Cylinder": "tab:green", "Housing": "0.3", "Entrance_Aperture": "tab:purple"}

n = len(a.heights)
fig, axes = plt.subplots(2, (n + 1) // 2, figsize=(5.2 * ((n + 1) // 2), 10.5))
axes = axes.ravel()
for ax, zb in zip(axes, a.heights):
    zb_m = zb * 1e-3
    for name, v in V.items():
        sel = np.abs(-v[:, 2] - zb_m) < 0.0025
        if sel.any():
            ax.plot(1e3 * v[sel, 0], 1e3 * v[sel, 1], ".", ms=3, color=colors[name], label=name)
    # orbit point at this height (interpolate along the orbit)
    k = int(np.argmin(np.abs(zb_trj - zb_m)))
    if abs(zb_trj[k] - zb_m) < 0.006:
        ax.plot(1e3 * trj[k, 0], 1e3 * trj[k, 1], "k*", ms=14, label="design orbit here")
        ax.annotate("orbit r = {:.1f} mm".format(1e3 * np.hypot(trj[k, 0], trj[k, 1])), (1e3 * trj[k, 0], 1e3 * trj[k, 1]),
                    textcoords="offset points", xytext=(8, 8), fontsize=9)
    ax.plot(1e3 * trj[:, 0], 1e3 * trj[:, 1], "-", color="0.6", lw=1, label="design orbit (projected)")
    ax.plot(0, 0, "k+", ms=12, mew=2, label="axis")
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(45 * np.cos(th), 45 * np.sin(th), ":", color="tab:green", lw=1, label="cylinder bore r = 45")
    ax.set_aspect("equal"); ax.set_xlim(-60, 60); ax.set_ylim(-60, 60); ax.grid(alpha=0.3)
    ax.set_title("BCS z = +{:.0f} mm  ({:.0f} mm below the entrance plate)".format(zb, 113 - zb))
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
axes[0].legend(fontsize=8, loc="upper left")
fig.suptitle("{}: electrode cross-sections near the entrance, R = {:+.2f} deg (BCS top view; the entrance plate's inner face is at z = +113 mm)".format(
    os.path.basename(run), R), fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = a.out or os.path.join(run, "fig_entrance_slices.png")
fig.savefig(out, dpi=110)
print("wrote", out)
