"""Render an inflector + quadrupole assembly with the tracked trajectories, in the Baseline (machine) frame.

Two figures in the style of the HCHC-60 audit renders (Results/audit/fig_doublet_triplet):
  <out>_render.png  3-D view: spiral anode red, cathode blue, quad poles orange, grounded plates grey, housing light grey;
                    trajectories green when transmitted, red when lost (the loss point appended)
  <out>_side.png    the x-z and y-z projections: each electrode as one translucent silhouette, all tracks on top

Everything is in BCS millimetres: +z up, the beam enters from +z and leaves at the median plane z = 0. The geometry is
put where the TRACKING had it - the levelled STEP folder rotated by the run's assembly rotation (not the re-solved export
rotation), each quad by its own total rotation on top, shifted by the run's levelling shift, mirrored and scaled - so the
tracks and the metal agree exactly. All recorded tracks are drawn (bunch.py records an evenly spaced subsample of every
particle, core and tail).

    python render_inflector.py <bunch_sc_all_level.npz> --steps <run>/sc_level/t<trunc>/steps
                               --quad-rotations 19.49 21.64 26.88 --out <png stem> [--title ...]
"""
import argparse
import glob
import json
import os
import tempfile

import numpy as np


def _cad_files(folder):
    """Electrode CAD files of a STEP folder: NNN_<Name>.stp (or older .step); assemblies and _negative solids
    (the gap volume) are not electrodes."""
    out = []
    for f in sorted(os.listdir(folder)):
        stem, ext = os.path.splitext(f)
        if ext.lower() in (".stp", ".step") and "assembly" not in stem.lower() and not stem.lower().endswith("_negative"):
            out.append(os.path.join(folder, f))
    return out


def _cad_name(path):
    """'NNN_<Name>.stp' -> '<Name>'."""
    return os.path.splitext(os.path.basename(path))[0].split("_", 1)[1]

p = argparse.ArgumentParser()
p.add_argument("npz")
p.add_argument("--steps", required=True, help="the levelled STEP folder the run tracked in (deck frame, metres)")
p.add_argument("--quad-rotations", type=float, nargs="+", required=True, help="TOTAL quad rotations on top of R [deg]")
p.add_argument("--rotation", type=float, default=None, help="assembly rotation [deg] (default: the run json's)")
p.add_argument("--shift-z", type=float, default=None, help="levelling shift [m, deck z] (default: the run json's)")
p.add_argument("--out", required=True, help="output stem; writes <stem>_render.png and <stem>_side.png")
p.add_argument("--title", default=None)
p.add_argument("--mesh-mm", type=float, default=1.5)
p.add_argument("--elev", type=float, default=22.0)
p.add_argument("--azim", type=float, default=-58.0)
p.add_argument("--current-ma", type=float, default=10.0)
p.add_argument("--zmin", type=float, default=-25.0, help="lowest height drawn [mm, BCS]; tracks continue past the exit")
p.add_argument("--n-tracks", type=int, default=0,
               help="draw this many evenly spaced recorded tracks (0 = all); 1000 on a full record picks the same particles "
                    "as bunch.py's default 1000-track recorder")
a = p.parse_args()
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

js = json.load(open(a.npz[:-4] + ".json"))
R = js.get("assembly_rotation_deg", 0.0) if a.rotation is None else a.rotation
shift = js.get("shift_z_m", 0.0) if a.shift_z is None else a.shift_z

# ------------------------------------------------------------------ geometry in BCS mm, as tracked
from spyral_inflector.tracking.export import steps_to_baseline_mm  # noqa: E402
import gmsh  # noqa: E402

tmp = tempfile.mkdtemp(prefix="render_bcs_")
steps_to_baseline_mm(a.steps, tmp, rotation_deg=R, quad_rotation=tuple(a.quad_rotations), shift_z=shift,
                     combined=None, log=lambda m: None)
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
meshes = []
for f in _cad_files(tmp):
    name = _cad_name(f)
    gmsh.model.add(name)
    gmsh.model.occ.importShapes(f)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", a.mesh_mm)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.3 * a.mesh_mm)
    gmsh.model.mesh.generate(2)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    lut = {int(t): i for i, t in enumerate(tags)}
    tris = []
    for etype, _etags, enodes in zip(*gmsh.model.mesh.getElements(2)):
        if etype == 2:                                   # 3-node triangles
            tris.append(np.array([lut[int(n)] for n in enodes]).reshape(-1, 3))
    gmsh.model.remove()
    if tris:
        meshes.append((name, xyz, np.vstack(tris)))
gmsh.finalize()


def style(name):
    n = name.lower()
    if name == "SI_Anode":
        return "#d6614d", 0.55
    if name == "SI_Cathode":
        return "#4f8fc0", 0.55
    if name.startswith("D") and name[1:].isdigit():
        return "#f08a24", 0.55
    if "housing" in n:
        return "#c8c8c8", 0.10
    return "#8c8c8c", 0.40                                # grounded plates / apertures


# ------------------------------------------------------------------ trajectories in BCS mm
d = np.load(a.npz, allow_pickle=True)
traj, idx = d["traj"], d["traj_idx"]
if a.n_tracks and len(idx) > a.n_tracks:
    _sel = np.unique(np.linspace(0, len(idx) - 1, a.n_tracks).astype(int))
    traj, idx = traj[:, _sel, :], idx[_sel]
crossed, hit, hit_point = d["crossed"], d["hit_electrode"], d["hit_point"]
post = d["post_exit_hit"] if "post_exit_hit" in d.files else np.zeros(len(crossed), bool)
post_point = d["post_exit_point"] if "post_exit_point" in d.files else None
tracks = []
for j, pidx in enumerate(idx):
    x = traj[:, j, :].astype(float)
    x = x[np.isfinite(x[:, 0])]
    if post[pidx]:
        fate = "lost"
        if post_point is not None and np.all(np.isfinite(post_point[pidx])):
            x = np.vstack([x, post_point[pidx]])
    elif crossed[pidx]:
        fate = "ok"
    elif hit[pidx] >= 0:
        fate = "lost"
        if np.all(np.isfinite(hit_point[pidx])):
            x = np.vstack([x, hit_point[pidx]])
    else:
        fate = "other"
    if len(x) > 1:
        x = 1e3 * x * np.array([1.0, 1.0, -1.0])         # deck frame (m) -> BCS (mm)
        tracks.append((x, fate))
n_ok = sum(f == "ok" for _, f in tracks)
n_lost = sum(f == "lost" for _, f in tracks)
OK, LOST, OTHER = "#2f9e44", "#e03131", "#868e96"
title = a.title or os.path.basename(os.path.dirname(os.path.dirname(a.npz)))
sub = "{:g} mA, all particles (core + RFQ tail), {:.1f} % transmitted; green transmitted, red lost ({} recorded tracks)".format(
    a.current_ma, 100 * js["transmission"], len(tracks))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

# ------------------------------------------------------------------ 3-D view
fig = plt.figure(figsize=(10, 11.5))
# the beam is the point of the picture: draw every track above the translucent metal instead of letting
# matplotlib's per-artist depth sort hide it inside the quadrupoles
ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
order = sorted(meshes, key=lambda m: ("housing" not in m[0].lower(), m[0]))   # housing first (behind)
for name, v, f in order:
    col, al = style(name)
    ax.add_collection3d(Poly3DCollection(v[f], facecolor=col, edgecolor="none", alpha=al, linewidths=0, zorder=1))
tracks.sort(key=lambda t: {"lost": 0, "other": 1, "ok": 2}[t[1]])       # lost first, transmitted drawn over them
ZMIN = a.zmin
for x, fate in tracks:
    x = x[x[:, 2] >= ZMIN]
    if len(x) < 2:
        continue
    col, al, lw = {"ok": (OK, 0.45, 0.4), "lost": (LOST, 0.55, 0.45)}.get(fate, (OTHER, 0.5, 0.4))
    ax.plot(x[:, 0], x[:, 1], x[:, 2], color=col, lw=lw, alpha=al, zorder=3 if fate == "ok" else 2)
# the frame is the metal plus the transmitted and lost tracks: a handful of tail particles that turn round and leave
# upstream ('other', still alive at the end) ran to z = 900 mm and squashed the views (2026-10-20); they are still
# drawn, clipped at the frame
allv = np.vstack([v for _n, v, _f in meshes] + [x[x[:, 2] >= ZMIN] for x, fate in tracks if fate != "other"])
lo, hi = allv.min(axis=0), allv.max(axis=0)
ax.set_xlim(lo[0], hi[0])
ax.set_ylim(lo[1], hi[1])
ax.set_zlim(lo[2], hi[2])
ax.set_box_aspect(hi - lo)
ax.view_init(elev=a.elev, azim=a.azim)
ax.set_xlabel("x [mm]", fontsize=13)
ax.set_ylabel("y [mm]", fontsize=13)
ax.set_zlabel("height above midplane [mm]", fontsize=13)
ax.set_title("{}\n{}".format(title, sub), fontsize=12)
for lab, col in (("transmitted ({})".format(n_ok), OK), ("lost ({})".format(n_lost), LOST)):
    ax.plot([], [], color=col, lw=1.5, label=lab)
for lab, col in (("anode (+V)", "#d6614d"), ("cathode (-V)", "#4f8fc0"), ("quadrupole poles", "#f08a24"),
                 ("grounded plates", "#8c8c8c")):
    ax.plot([], [], color=col, lw=7, alpha=0.6, label=lab)
ax.legend(loc="upper left", fontsize=9, frameon=False)
fig.tight_layout()
fig.savefig(a.out + "_render.png", dpi=110)
plt.close(fig)
print("wrote", a.out + "_render.png", "({} transmitted, {} lost tracks)".format(n_ok, n_lost))

# ------------------------------------------------------------------ side views: x-z and y-z silhouettes
from PIL import Image, ImageDraw  # noqa: E402

PX = 0.25                                                # mm per pixel of the silhouette masks
fig, axs = plt.subplots(1, 2, figsize=(10.5, 9.5), sharey=True, gridspec_kw={"wspace": 0.08})
zlo, zhi = max(lo[2] - 2, ZMIN), hi[2] + 2
for k, (axis, lab) in enumerate(((0, "x [mm]"), (1, "y [mm]"))):
    ax = axs[k]
    hlo, hhi = lo[axis] - 2, hi[axis] + 2
    W, H = int((hhi - hlo) / PX) + 1, int((zhi - zlo) / PX) + 1
    for name, v, f in order:
        col, al = style(name)
        if "housing" in name.lower():
            al = 0.18
        img = Image.new("L", (W, H), 0)
        dr = ImageDraw.Draw(img)
        P = np.column_stack([(v[:, axis] - hlo) / PX, (zhi - v[:, 2]) / PX])
        for tri in f:
            dr.polygon([tuple(P[t]) for t in tri], fill=255)
        mask = np.asarray(img, dtype=float) / 255.0
        rgba = np.zeros((H, W, 4))
        rgba[..., :3] = matplotlib.colors.to_rgb(col)
        rgba[..., 3] = al * mask
        ax.imshow(rgba, extent=(hlo, hhi, zlo, zhi), origin="upper", interpolation="nearest", aspect="equal", zorder=1)
    for x, fate in tracks:
        col, al, lw = {"ok": (OK, 0.5, 0.4), "lost": (LOST, 0.6, 0.4)}.get(fate, (OTHER, 0.5, 0.4))
        ax.plot(x[:, axis], x[:, 2], color=col, lw=lw, alpha=al, zorder=3 if fate == "ok" else 2)
    ax.axhline(0.0, color="k", lw=0.6, ls=":", zorder=4)
    ax.set_xlim(hlo, hhi)
    ax.set_ylim(zlo, zhi)
    ax.set_xlabel(lab, fontsize=12)
    ax.grid(True, alpha=0.25)
axs[0].set_ylabel("height above midplane [mm]", fontsize=12)
fig.suptitle("{}: side views\n{}".format(title, sub), fontsize=12)
fig.tight_layout()
fig.savefig(a.out + "_side.png", dpi=120)
plt.close(fig)
print("wrote", a.out + "_side.png")
