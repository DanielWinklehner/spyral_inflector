"""Convert a COMSOL h5part field map of the HCHC-60 magnet into the field pickle the
tracking scripts load, in the deck frame: beam travels in +z from the bore (negative z)
to the median plane at z = 0, on-axis Bz negative (-0.891 T at the centre).

Two export conventions are supported (--mode):
  rotate  the map covers the bore side of the model (z > 0 there): rotate 180 deg about x,
          (x, y, z) -> (x, -y, -z), B -> (Bx, -By, -Bz). Proper rotation, no symmetry
          assumed. This is what flip_bfield_z.py did.
  negate  the map already covers z < 0 in the model frame (the other half of a
          median-plane-symmetric magnet): keep the axes, negate all three components.
          Identical to 'rotate' when the model is mirror-symmetric about the median
          plane and one vertical plane, which the 2026-09 map satisfies to < 10 uT.
  none    the export is already in the deck frame.

--fix-seams replaces the values on the grid planes x = 0, y = 0 and z = 0 by the mean
of the two neighbouring planes. A model built with symmetry planes and exported on a
grid that has points exactly on them carries interpolation glitches there (up to 2 T at
iron boundaries, 0.1-20 mT near the beam); the average across the seam is the
symmetric value (zero for the component that is odd across that plane).

--crop-xy limits the transverse extent (the housing needs about +-0.11 m).

    python convert_bfield.py SRC.h5part DST.pickle --mode negate --fix-seams --crop-xy 0.15
"""
import argparse
import os
import time

import h5py
import numpy as np

from PyPATools.field import Field

parser = argparse.ArgumentParser()
parser.add_argument("src")
parser.add_argument("dst")
parser.add_argument("--mode", choices=["rotate", "negate", "none", "auto"], default="auto",
                    help="auto: the map is taken to cover the bore at negative z already, and is negated "
                         "only if its on-axis Bz(0,0,0) is positive")
parser.add_argument("--fix-seams", action="store_true")
parser.add_argument("--crop-xy", type=float, default=None, help="keep |x|, |y| <= this [m]")
parser.add_argument("--group", default="Step#0/Block/Hfield")
args = parser.parse_args()

def read_h5part(path, group_name):
    with h5py.File(path, "r") as fh:
        group = fh[group_name]
        origin = np.asarray(group.attrs["__Origin__"], dtype=float)
        spacing = np.asarray(group.attrs["__Spacing__"], dtype=float)
        nz, ny, nx = group["0"].shape
        raw = {}
        for comp, dset in (("x", "0"), ("y", "1"), ("z", "2")):
            # file order (nz, ny, nx) -> (nx, ny, nz)
            raw[comp] = np.ascontiguousarray(np.transpose(np.asarray(group[dset]), (2, 1, 0)))
    return (origin[0] + np.arange(nx) * spacing[0], origin[1] + np.arange(ny) * spacing[1],
            origin[2] + np.arange(nz) * spacing[2], raw)


def read_comsol_text(path):
    """COMSOL spreadsheet export: '%' header lines, then one node per line
    'x y z Bx By Bz' (metres, Tesla) with x varying fastest, then y, then z."""
    import pandas as pd
    n_nodes = None
    with open(path, "r") as fh:
        for line in fh:
            if not line.startswith("%"):
                break
            if "Nodes:" in line:
                n_nodes = int(line.split(":")[1])
            if any(k in line for k in ("Length unit", "Description", "Dimension", "Model:")):
                print("  " + line.strip())
    df = pd.read_csv(path, sep=r"\s+", comment="%", header=None,
                     names=["x", "y", "z", "bx", "by", "bz"], dtype=np.float64, engine="c")
    print("  {:,d} nodes read ({:.0f} s)".format(len(df), time.time() - t0), flush=True)
    assert n_nodes is None or len(df) == n_nodes, "node count does not match the header"
    coords = {c: np.unique(np.round(df[c].to_numpy(), 6)) for c in "xyz"}
    nx, ny, nz = (len(coords[c]) for c in "xyz")
    assert nx * ny * nz == len(df), "not a full regular grid ({} x {} x {} != {})".format(nx, ny, nz, len(df))
    for c in "xyz":
        d = np.diff(coords[c])
        assert np.allclose(d, d[0], atol=1e-7), "{} spacing is not uniform".format(c)
    assert np.allclose(np.round(df["x"].to_numpy()[:nx], 6), coords["x"]), "x is not the fastest index"
    assert abs(df["y"].to_numpy()[nx] - coords["y"][1]) < 1e-6, "y is not the second index"
    assert abs(df["z"].to_numpy()[nx * ny] - coords["z"][1]) < 1e-6, "z is not the slowest index"
    raw = {}
    for comp, col in (("x", "bx"), ("y", "by"), ("z", "bz")):
        raw[comp] = np.ascontiguousarray(np.transpose(df[col].to_numpy().reshape(nz, ny, nx), (2, 1, 0)))
    del df
    return coords["x"], coords["y"], coords["z"], raw


t0 = time.time()
print("Reading {}".format(os.path.basename(args.src)))
if args.src.lower().endswith((".h5part", ".h5")):
    x, y, z, raw = read_h5part(args.src, args.group)
else:
    x, y, z, raw = read_comsol_text(args.src)
nx, ny, nz = len(x), len(y), len(z)
print("  source grid: x {:+.4f}..{:+.4f} ({}), y {:+.4f}..{:+.4f} ({}), z {:+.4f}..{:+.4f} ({}), spacing {:.3f} mm ({:.0f} s)".format(
    x[0], x[-1], nx, y[0], y[-1], ny, z[0], z[-1], nz, 1e3 * (x[1] - x[0]), time.time() - t0), flush=True)
bz_src = raw["z"][int(np.argmin(np.abs(x))), int(np.argmin(np.abs(y))), int(np.argmin(np.abs(z)))]
print("  source on-axis Bz(0,0,0) = {:+.5f} T (before any sign change)".format(bz_src))
if args.mode == "auto":
    assert z[0] < -0.1, "auto mode expects the bore at negative z in the source"
    args.mode = "negate" if bz_src > 0 else "none"
    print("  mode auto -> {}".format(args.mode))

if args.mode == "rotate":
    assert np.allclose(-y[::-1], y), "y axis is not symmetric; the rotation would move it"
    vals = {"x": np.ascontiguousarray(raw["x"][:, ::-1, ::-1]),
            "y": np.ascontiguousarray(-raw["y"][:, ::-1, ::-1]),
            "z": np.ascontiguousarray(-raw["z"][:, ::-1, ::-1])}
    z = -z[::-1]
elif args.mode == "negate":
    vals = {k: -raw[k] for k in "xyz"}
else:
    vals = raw
del raw

if args.crop_xy is not None:
    mx = np.abs(x) <= args.crop_xy + 1e-9
    my = np.abs(y) <= args.crop_xy + 1e-9
    vals = {k: np.ascontiguousarray(v[mx][:, my]) for k, v in vals.items()}
    x, y = x[mx], y[my]

if args.fix_seams:
    # Seam planes x = 0, y = 0, z = 0. A point on one seam gets the mean of its two
    # neighbours across that plane; a point on a seam LINE (two planes) the mean of the
    # four diagonal neighbours off both planes; the origin the mean of the eight corner
    # neighbours. Neighbours on another seam are never used, and every replacement is
    # computed from the original values, so glitches do not propagate.
    seams = []
    for axis, coord in ((0, x), (1, y), (2, z)):
        i0 = int(np.argmin(np.abs(coord)))
        if abs(coord[i0]) < 1e-6 and 0 < i0 < len(coord) - 1:
            seams.append((axis, i0))
        else:
            print("  no grid plane on {} = 0, nothing to fix there".format("xyz"[axis]))
    for comp in "xyz":
        a = vals[comp]
        src = a.copy()
        changed, worst = 0, 0.0
        n_seam = len(seams)
        # enumerate all non-empty subsets of the seams: a point lying on exactly that subset
        for mask in range(1, 1 << n_seam):
            on = [seams[b] for b in range(n_seam) if mask >> b & 1]
            off = [seams[b] for b in range(n_seam) if not mask >> b & 1]
            # index grid for the points on exactly these seams
            idx = [slice(None)] * 3
            for axis, i0 in on:
                idx[axis] = i0
            sub = src[tuple(idx)]
            new = np.zeros_like(sub)
            for signs in np.array(np.meshgrid(*[[-1, 1]] * len(on), indexing="ij")).reshape(len(on), -1).T:
                nb = [slice(None)] * 3
                for (axis, i0), s in zip(on, signs):
                    nb[axis] = i0 + int(s)
                new += src[tuple(nb)]
            new /= 2 ** len(on)
            if np.ndim(sub) == 0:
                # all seams at once: the origin
                d0 = abs(float(sub) - float(new))
                changed += int(d0 > 1e-3)
                worst = max(worst, d0)
                a[tuple(idx)] = float(new)
                continue
            # exclude the points that also lie on the other seams (handled by their own subset)
            keep = np.ones(sub.shape, dtype=bool)
            for axis, i0 in off:
                # position of that axis among the remaining (non-fixed) axes of `sub`
                rem = [ax for ax in range(3) if ax not in [o[0] for o in on]]
                sl = [slice(None)] * sub.ndim
                sl[rem.index(axis)] = i0
                keep[tuple(sl)] = False
            d = np.abs(sub - new)[keep]
            changed += int((d > 1e-3).sum())
            worst = max(worst, float(d.max()) if d.size else 0.0)
            target = a[tuple(idx)]
            target[keep] = new[keep]
            a[tuple(idx)] = target
        print("  seams, B{}: {} points moved by more than 1 mT, largest change {:.1f} mT".format(comp, changed, 1e3 * worst))
    i0, j0, k0 = (int(np.argmin(np.abs(c))) for c in (x, y, z))
    print("  Bz at the origin: {:+.5f} T (was {:+.5f} T)".format(vals["z"][i0, j0, k0], bz_src))

field = Field.from_arrays(grid={"x": x, "y": y, "z": z}, values=vals,
                          label="HCHC-60 B-field, {} ({}seams fixed)".format(args.mode, "" if args.fix_seams else "no "),
                          dim=3, units="m", interpolator_backend="scipy")
print("  deck grid: x {:+.4f}..{:+.4f}, y {:+.4f}..{:+.4f}, z {:+.4f}..{:+.4f}".format(x[0], x[-1], y[0], y[-1], z[0], z[-1]))
print("  on-axis check:")
for zz in (z[0], -0.30, -0.25, -0.20, -0.10, -0.05, 0.0, 0.05):
    if z[0] <= zz <= z[-1]:
        print("    z = {:+.3f} m -> Bz = {:+.5f} T".format(zz, field(np.array([[0.0, 0.0, zz]]))[0, 2]))
bz0 = field(np.array([[0.0, 0.0, 0.0]]))[0, 2]
if bz0 > 0:
    print("  WARNING: Bz(0,0,0) is positive; the deck expects a negative field. Check --mode.")
os.makedirs(os.path.dirname(os.path.abspath(args.dst)), exist_ok=True)
field.save(args.dst)
print("  wrote {} ({:.2f} GB, {:.0f} s)".format(args.dst, os.path.getsize(args.dst) / 1e9, time.time() - t0))
