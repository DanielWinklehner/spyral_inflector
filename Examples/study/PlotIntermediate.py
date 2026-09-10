"""Intermediate picture of a final-push run: the geometry at its current rotation over the
field map (BCS frame, top view) at three heights, the design orbit, the exit plane and the
first gap; plus the on-axis field and the rotation history.

    python PlotIntermediate.py --run Results/final/final1 [--out fig.png]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spyral_inflector.tracking.deck import load_step_assembly, mesh_assembly, load_bfield, load_state
from spyral_inflector.tracking.exit_plane import rotz

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DECK = os.environ.get("SI_DECK") or os.path.dirname(SCRIPTS)
p = argparse.ArgumentParser()
p.add_argument("--run", required=True)
p.add_argument("--out", default=None)
a = p.parse_args()
run = a.run if os.path.isabs(a.run) else os.path.join(DECK, a.run)

spec_fn = os.path.join(run, "exit_plane.json") if os.path.exists(os.path.join(run, "exit_plane.json")) else os.path.join(run, "exit_plane_R0.json")
spec = json.load(open(spec_fn))
R = spec["rotation_deg"]
steps = os.path.join(run, "steps") if os.path.isdir(os.path.join(run, "steps")) and os.listdir(os.path.join(run, "steps")) else os.path.join(run, "steps0")
geo_dir = os.path.join(run, "geometry") if os.path.exists(os.path.join(run, "geometry", "state.pickle")) else os.path.join(run, "geometry0")
bf_path = json.load(open(os.path.join(run, "phase_1.done")))["info"]["pickle"]
bf = load_bfield(bf_path)
Rz = rotz(R)
st = load_state(os.path.join(geo_dir, "state.pickle"))
trj = np.asarray(st["trj_design"]) @ Rz.T
assy = load_step_assembly(steps)
mesh_assembly(assy)
V = {e.name: np.asarray(e._gmsh_msh["vertices"], float) @ Rz.T for e in assy.electrodes.values()}

hist = [json.load(open(os.path.join(run, "exit_plane_R0.json")))["rotation_deg"]] if os.path.exists(os.path.join(run, "exit_plane_R0.json")) else []
log = open(os.path.join(run, "log_run.txt"), encoding="utf-8").read().splitlines()
for line in log:
    if "rotation " in line and "->" in line and "deg (change" in line:
        try:
            hist.append(float(line.split("->")[1].split("deg")[0]))
        except ValueError:
            pass

fig = plt.figure(figsize=(18, 11))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.55])
ext = 0.15
xs = np.linspace(-ext, ext, 301)
X, Y = np.meshgrid(xs, xs, indexing="ij")
colors = {"SI_Anode": "tab:red", "SI_Cathode": "tab:blue", "Housing": "0.25", "Housing_Cylinder": "tab:green",
          "Entrance_Aperture": "tab:purple"}
for k, (zb, title) in enumerate(((0.0, "median plane, BCS z = 0"), (0.045, "BCS z = +45 mm (cutout level)"), (0.070, "BCS z = +70 mm (pole hole)"))):
    ax = fig.add_subplot(gs[0, k])
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, -zb)])          # deck z = -BCS z
    B = np.linalg.norm(bf(pts), axis=1).reshape(X.shape)
    im = ax.pcolormesh(1e3 * X, 1e3 * Y, B, cmap="viridis", vmin=0, vmax=2.4, shading="auto")
    ax.contour(1e3 * X, 1e3 * Y, B, levels=[1.4], colors="white", linewidths=0.8)
    for name, v in V.items():
        sel = np.abs(-v[:, 2] - zb) < 0.004                                        # vertices within 4 mm of this height
        if sel.any():
            ax.plot(1e3 * v[sel, 0], 1e3 * v[sel, 1], ".", ms=2, color=colors.get(name, "orange"), label=name if k == 0 else None)
    ax.plot(1e3 * trj[:, 0], 1e3 * trj[:, 1], "w-", lw=1.5, label="design orbit" if k == 0 else None)
    # exit plane and the first gap
    pb = np.array(spec["baseline"]["point_mm"]); nb = np.array(spec["baseline"]["normal"])
    tvec = np.array([-nb[1], nb[0]])
    ax.plot([pb[0] - 25 * tvec[0], pb[0] + 25 * tvec[0]], [pb[1] - 25 * tvec[1], pb[1] + 25 * tvec[1]], "y-", lw=2.5, label="exit plate outer face" if k == 0 else None)
    ax.plot([0, 150 * np.cos(np.radians(34))], [0, 150 * np.sin(np.radians(34))], "--", color="orange", lw=1.5, label="first gap, 34 deg" if k == 0 else None)
    ax.set_aspect("equal"); ax.set_xlim(-150, 150); ax.set_ylim(-150, 150)
    ax.set_title(title); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    if k == 0:
        ax.legend(loc="lower left", fontsize=8, framealpha=0.85)
    if k == 2:
        fig.colorbar(im, ax=ax, label="|B| [T]")

# on-axis field, new vs old map
ax = fig.add_subplot(gs[1, 0])
zz = np.linspace(-0.05, 0.40, 451)
Bn = bf(np.column_stack([np.zeros_like(zz), np.zeros_like(zz), -zz]))[:, 2]
ax.plot(1e3 * zz, Bn, label="new map (BCS)")
old = os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle")
if os.path.exists(old) and os.path.abspath(old) != os.path.abspath(bf_path):
    bo = load_bfield(old)
    Bo = bo(np.column_stack([np.zeros_like(zz), np.zeros_like(zz), -zz]))[:, 2]
    ax.plot(1e3 * zz, Bo, "--", label="previous map")
ax.axvline(275, color="0.5", lw=0.8); ax.text(277, ax.get_ylim()[0] if ax.get_ylim()[0] < 0 else -0.9, "beam start", fontsize=8, color="0.4")
ax.set_xlabel("BCS z [mm]"); ax.set_ylabel("Bz on axis [T]"); ax.set_title("on-axis field"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

# rotation history
ax = fig.add_subplot(gs[1, 1])
if hist:
    ax.plot(range(len(hist)), hist, "o-")
    for i, h in enumerate(hist):
        ax.annotate("{:+.2f}".format(h), (i, h), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=9)
ax.set_xlabel("iteration"); ax.set_ylabel("rotation R [deg]"); ax.set_title("rotation, re-optimized in the rotated field"); ax.grid(alpha=0.3)

# side view: r vs BCS z of the geometry and the orbit
ax = fig.add_subplot(gs[1, 2])
for name, v in V.items():
    if name in ("Housing", "Housing_Cylinder", "SI_Anode", "SI_Cathode", "Entrance_Aperture"):
        ax.plot(-1e3 * v[:, 2], 1e3 * np.hypot(v[:, 0], v[:, 1]), ".", ms=1.5, color=colors.get(name, "orange"), label=name)
ax.plot(-1e3 * trj[:, 2], 1e3 * np.hypot(trj[:, 0], trj[:, 1]), "k-", lw=1.5, label="design orbit")
ax.axhline(49, color="0.5", lw=0.8, ls=":"); ax.text(60, 50.5, "pole hole r = 49 mm", fontsize=8, color="0.4")
ax.set_xlabel("BCS z [mm]"); ax.set_ylabel("r [mm]"); ax.set_title("side view (radius vs height)"); ax.legend(fontsize=7, loc="upper right"); ax.grid(alpha=0.3)
ax.set_xlim(-30, 130)

opt = json.load(open(os.path.join(geo_dir, "summary.json")))["optimizer"] if os.path.exists(os.path.join(geo_dir, "summary.json")) else {}
fig.suptitle("{}: geometry of {} at R = {:+.3f} deg over map {}   |   optimizer: V {:.0f} V, dz {:+.2f} mm, exit angle {:+.3f} deg, z {:+.2f} mm".format(
    os.path.basename(run), os.path.basename(geo_dir), R, os.path.basename(bf_path), opt.get("voltage_V", float("nan")), opt.get("dz_mm", float("nan")),
    opt.get("residual_final", {}).get("angle_deg", float("nan")), opt.get("residual_final", {}).get("z_offset_mm", float("nan"))), fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = a.out or os.path.join(run, "fig_intermediate.png")
fig.savefig(out, dpi=110)
print("wrote", out)
