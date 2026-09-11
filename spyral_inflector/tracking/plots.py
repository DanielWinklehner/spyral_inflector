"""3-D plots of the solid electrode geometry with a subsample of particle trajectories:
an isometric view with -z upwards (the beam enters from the top and leaves at the median
plane) and a side view perpendicular to the design exit velocity with the median plane
drawn, to judge how level the exit beam is. Trajectories come from the `traj` array of a
bunch run's npz (an evenly spaced subsample, continued past the exit plane)."""
import argparse
import json
import os
import pickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

from .deck import load_step_assembly, QUAD_NAMES  # noqa: E402

COLORS = {"SI_Anode": ("#d62728", 0.30), "SI_Cathode": ("#1f77b4", 0.30)}
QUAD_COLOR = ("#ff7f0e", 0.10)
PLATE_COLOR = ("#7f7f7f", 0.10)
HOUSING_COLOR = ("#bbbbbb", 0.04)
OK_COLOR, LOST_COLOR, OTHER_COLOR = "#2ca02c", "#d62728", "#555555"


def electrode_meshes(step_dir, h=0.005, rotation_deg=0.0):
    """[(name, vertices (n,3), faces (m,3))] of every electrode in the STEP folder, the
    assembly rigidly rotated about z by rotation_deg first (the system rotation R of a
    run whose tracks were made with rotate_assembly; the STEP files stay unrotated)."""
    assembly = load_step_assembly(step_dir)
    if rotation_deg:
        from .deck import rotate_assembly
        rotate_assembly(assembly, float(rotation_deg))
    out = []
    for e in assembly.electrodes.values():
        if e.generate_mesh(brep_h=h) != 0:
            raise RuntimeError("failed to mesh " + e.name)
        m = e._gmsh_msh
        out.append((e.name, np.asarray(m["vertices"], dtype=float), np.asarray(m["elements"], dtype=int)))
    return out


def _style(name):
    if name in COLORS:
        return COLORS[name]
    if name in QUAD_NAMES:
        return QUAD_COLOR
    if "housing" in name.lower():
        return HOUSING_COLOR
    return PLATE_COLOR


def load_trajectories(npz_path, n_traj):
    """[(points (k,3), fate)] for an evenly spaced subsample of the recorded trajectories."""
    d = np.load(npz_path, allow_pickle=True)
    if "traj" not in d or d["traj"].shape[0] == 0:
        raise RuntimeError("{} has no recorded trajectories".format(npz_path))
    traj, idx = d["traj"], d["traj_idx"]
    crossed, hit, hit_point = d["crossed"], d["hit_electrode"], d["hit_point"]
    post = d["post_exit_hit"] if "post_exit_hit" in d.files else np.zeros(len(crossed), dtype=bool)
    post_point = d["post_exit_point"] if "post_exit_point" in d.files else None
    keep = np.unique(np.linspace(0, len(idx) - 1, min(n_traj, len(idx))).astype(int))
    out = []
    for j in keep:
        p = idx[j]
        x = traj[:, j, :].astype(float)
        x = x[np.isfinite(x[:, 0])]
        if post[p]:
            fate = "lost"
            if post_point is not None and np.all(np.isfinite(post_point[p])):
                x = np.vstack([x, post_point[p]])
        elif crossed[p]:
            fate = "ok"
        elif hit[p] >= 0:
            fate = "lost"
            if np.all(np.isfinite(hit_point[p])):
                x = np.vstack([x, hit_point[p]])
        else:
            fate = "other"
        if len(x) > 1:
            out.append((x, fate))
    return out, d


def _draw(ax, meshes, trajs, show_housing=False, zclip=None, lw=0.5):
    for name, verts, faces in meshes:
        if "housing" in name.lower() and not show_housing:
            continue
        if zclip is not None:
            faces = faces[np.all(verts[faces][:, :, 2] >= zclip, axis=1)]
            if len(faces) == 0:
                continue
        color, alpha = _style(name)
        ax.add_collection3d(Poly3DCollection(verts[faces], facecolor=color, edgecolor="none", alpha=alpha, linewidths=0))
    n_ok = n_lost = 0
    for x, fate in trajs:
        if zclip is not None:
            x = x[x[:, 2] >= zclip - 0.01]
            if len(x) < 2:
                continue
        col, al = {"ok": (OK_COLOR, 0.7), "lost": (LOST_COLOR, 0.6)}.get(fate, (OTHER_COLOR, 0.5))
        ax.plot(x[:, 0], x[:, 1], x[:, 2], color=col, lw=lw, alpha=al)
        n_ok += fate == "ok"
        n_lost += fate == "lost"
    return n_ok, n_lost


def _legend(ax, n_ok, n_lost):
    ax.plot([], [], color=OK_COLOR, label="transmitted ({})".format(n_ok))
    ax.plot([], [], color=LOST_COLOR, label="lost ({})".format(n_lost))
    ax.plot([], [], color=COLORS["SI_Anode"][0], lw=6, alpha=0.5, label="anode (+V)")
    ax.plot([], [], color=COLORS["SI_Cathode"][0], lw=6, alpha=0.5, label="cathode (-V)")
    ax.plot([], [], color=QUAD_COLOR[0], lw=6, alpha=0.5, label="quadrupoles")
    ax.plot([], [], color=PLATE_COLOR[0], lw=6, alpha=0.5, label="grounded apertures")
    ax.legend(loc="upper left", fontsize=9)


def plot_geometry_trajectories(npz_path, step_dir, out_png, title=None, n_traj=300, show_housing=False,
                               meshes=None, elev=35.264, azim=135):
    """Isometric view of the whole system, -z upwards. Returns (n_ok, n_lost) drawn."""
    trajs, _ = load_trajectories(npz_path, n_traj)
    if meshes is None:
        meshes = electrode_meshes(step_dir)
    fig = plt.figure(figsize=(11, 11))
    ax = fig.add_subplot(111, projection="3d")
    n_ok, n_lost = _draw(ax, meshes, trajs, show_housing)
    allv = np.vstack([v for name, v, _f in meshes if show_housing or "housing" not in name.lower()] + [x for x, _ in trajs])
    lo, hi = allv.min(axis=0), allv.max(axis=0)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect((hi - lo))
    ax.invert_zaxis()                      # -z upwards: the beam enters at the top
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m], -z up")
    _legend(ax, n_ok, n_lost)
    ax.set_title(title or os.path.basename(npz_path), fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)
    return n_ok, n_lost


def plot_side_view(npz_path, step_dir, out_png, title=None, n_traj=300, meshes=None, zmin=-0.13, zmax=0.05, azim=None):
    """Side view: horizontal line of sight perpendicular to the design exit velocity (from
    the si_state pickle next to the npz when available), orthographic, from z = zmin down
    past the median plane, with z = 0 drawn."""
    trajs, d = load_trajectories(npz_path, n_traj)
    if meshes is None:
        meshes = electrode_meshes(step_dir)
    if azim is None:
        azim = 90.0
        base = os.path.basename(npz_path)
        state_fn = os.path.join(os.path.dirname(npz_path), "si_state_{}.pickle".format(base[len("bunch_"):-len(".npz")]))
        try:
            with open(state_fn, "rb") as fh:
                vd = np.asarray(pickle.load(fh)["v_design"])[-1]
            azim = np.degrees(np.arctan2(vd[1], vd[0])) + 90.0     # look perpendicular to the exit velocity
        except Exception:  # noqa: BLE001
            pass
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")
    n_ok, n_lost = _draw(ax, meshes, trajs, show_housing=False, zclip=zmin, lw=0.6)
    allv = np.vstack([v for name, v, _f in meshes if "housing" not in name.lower()] + [x for x, _ in trajs])
    lo, hi = allv.min(axis=0), allv.max(axis=0)
    xx = [lo[0], hi[0], hi[0], lo[0], lo[0]]
    yy = [lo[1], lo[1], hi[1], hi[1], lo[1]]
    ax.plot(xx, yy, [0.0] * 5, color="k", lw=1.0, alpha=0.8)
    ax.plot([], [], color="k", lw=1.0, label="median plane z = 0")
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(zmin, zmax)
    try:
        ax.set_box_aspect((hi[0] - lo[0], hi[1] - lo[1], zmax - zmin), zoom=1.25)
    except TypeError:                       # matplotlib < 3.6
        ax.set_box_aspect((hi[0] - lo[0], hi[1] - lo[1], zmax - zmin))
    ax.invert_zaxis()
    ax.set_proj_type("ortho")
    ax.view_init(elev=0.0, azim=azim)
    along_y = abs(np.sin(np.radians(azim))) >= abs(np.cos(np.radians(azim)))
    if along_y:
        ax.set_yticks([])
        ax.set_xlabel("x [m]")
        ax.set_ylabel("")
    else:
        ax.set_xticks([])
        ax.set_xlabel("")
        ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m], -z up")
    _legend(ax, n_ok, n_lost)
    ax.set_title((title or os.path.basename(npz_path)) + "\nside view, line of sight perpendicular to the exit velocity",
                 fontsize=10, y=1.06)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)
    return n_ok, n_lost


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("npz")
    p.add_argument("--step-dir", required=True)
    p.add_argument("--out", required=True, help="isometric view png")
    p.add_argument("--side", default=None, help="side view png (optional)")
    p.add_argument("--title", default=None)
    p.add_argument("--n-traj", type=int, default=300)
    p.add_argument("--housing", action="store_true")
    p.add_argument("--elev", type=float, default=35.264)
    p.add_argument("--azim", type=float, default=135)
    p.add_argument("--side-azim", type=float, default=None, help="default: perpendicular to the design exit velocity")
    p.add_argument("--rotate", type=float, default=None,
                   help="rigid rotation of the STEP assembly about z [deg] to match tracks made with a rotated system; "
                        "default: the assembly_rotation_deg recorded in the npz's json summary, else 0")
    a = p.parse_args(argv)
    js = os.path.splitext(a.npz)[0] + ".json"
    title = a.title
    if title is None and os.path.exists(js):
        s = json.load(open(js))
        title = "{}: {:,d} particles, transmission {:.1f} %".format(os.path.basename(a.npz), s["n_particles"], 100 * s["transmission"])
    rotation = a.rotate
    if rotation is None and os.path.exists(js):
        _s = json.load(open(js))
        rotation = float(_s.get("assembly_rotation_deg", _s.get("rotation_deg", 0.0)) or 0.0)
    if rotation:
        print("assembly rotated by {:+.4f} deg about z to match the tracks".format(rotation))
    meshes = electrode_meshes(a.step_dir, rotation_deg=rotation or 0.0)
    n_ok, n_lost = plot_geometry_trajectories(a.npz, a.step_dir, a.out, title=title, n_traj=a.n_traj, show_housing=a.housing,
                                              meshes=meshes, elev=a.elev, azim=a.azim)
    print("wrote {} ({} transmitted, {} lost trajectories drawn)".format(a.out, n_ok, n_lost))
    if a.side:
        plot_side_view(a.npz, a.step_dir, a.side, title=title, n_traj=a.n_traj, meshes=meshes, azim=a.side_azim)
        print("wrote {}".format(a.side))


if __name__ == "__main__":
    main()
