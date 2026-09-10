"""Coordinate frames of the B-field, and rotating it about the axis.

The deck integrates in the mirror image of the machine through the median plane: the
beam is tracked upward in +z (see export.py). Field maps are now delivered in the
machine ("Baseline") frame: right-handed, +z up, the beam entering from +z. B is a
pseudovector, so the mirror z -> -z takes (Bx, By, Bz)(x, y, z) to (-Bx, -By, +Bz)(x, y, -z).
For a median-plane-symmetric magnet that is the same field -- but the map's z-range lies
on the other side of the plane, so the relabeling is still needed to cover the bore.

Which frame a map is in is read off its z-range: the bore (beam) side is at positive z
in the Baseline frame and at negative z in the deck frame. The file on disk is never
changed; the mirror happens in memory at load.
"""
import numpy as np

from PyPATools.field import Field


def bfield_frame(field):
    """'baseline' if the map extends further towards +z than -z, 'deck' for the reverse,
    None for a field without a grid (a constant) or a z-range too symmetric to tell."""
    grid = getattr(field, "grid", None)
    if grid is None:
        return None
    z = np.asarray(grid["z"], dtype=float)
    lo, hi = float(z.min()), float(z.max())
    if abs(hi + lo) < 0.02 * (hi - lo):
        return None
    return "baseline" if hi > -lo else "deck"


def mirror_bfield_z(field, label=None):
    """The mirror image of a gridded B-field through the plane z = 0 (pseudovector rule)."""
    g, v = field.grid, field.grid_values
    z = np.asarray(g["z"], dtype=float)
    order = np.argsort(-z)
    grid = {"x": np.asarray(g["x"], dtype=float), "y": np.asarray(g["y"], dtype=float), "z": -z[order]}
    values = {"x": -np.asarray(v["x"])[:, :, order], "y": -np.asarray(v["y"])[:, :, order],
              "z": np.asarray(v["z"])[:, :, order]}
    return Field.from_arrays(grid=grid, values=values, dim=3, units="m",
                             label=label or "{} (z-mirrored)".format(getattr(field, "label", "B-field")))


def to_deck_frame(field, frame="auto", log=print):
    """The field in the deck frame. frame='auto' reads the frame off the z-range,
    'baseline' forces the mirror, 'deck' returns the field unchanged. Under 'auto' the
    call is idempotent: a field already in the deck frame is left alone."""
    if frame == "deck" or getattr(field, "grid", None) is None:
        return field
    if frame == "auto":
        detected = bfield_frame(field)
        if detected is None:
            log("B-field frame: symmetric z-range, taken as deck frame (no mirror)")
            return field
        frame = detected
    if frame == "baseline":
        z = np.asarray(field.grid["z"], dtype=float)
        log("B-field frame: Baseline (z {:+.3f}..{:+.3f} m) -> mirrored into the deck frame".format(z.min(), z.max()))
        return mirror_bfield_z(field)
    return field


def rotate_bfield_z(field, deg, label=None, log=print):
    """The field a system rotated by +deg about z sees, on the field's own grid:
    B'(x) = Rz(-deg) B(Rz(deg) x). A rotation about z never mixes z-slices, so each slice
    is resampled in 2-D (bilinear); positions that fall outside the grid take the nearest
    edge value, which only affects the corners beyond the tracked region."""
    from scipy.ndimage import map_coordinates

    g, v = field.grid, field.grid_values
    x, y, z = (np.asarray(g[k], dtype=float) for k in "xyz")
    bx, by, bz = (np.asarray(v[k], dtype=float) for k in "xyz")
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    X, Y = np.meshgrid(x, y, indexing="ij")
    xr, yr = c * X - s * Y, s * X + c * Y                       # Rz(+deg) applied to the grid point
    ix = (xr - x[0]) / (x[1] - x[0])
    iy = (yr - y[0]) / (y[1] - y[0])
    coords = np.stack([ix.ravel(), iy.ravel()])
    out = {k: np.empty_like(bx) for k in "xyz"}
    for k in range(len(z)):
        sx = map_coordinates(bx[:, :, k], coords, order=1, mode="nearest").reshape(X.shape)
        sy = map_coordinates(by[:, :, k], coords, order=1, mode="nearest").reshape(X.shape)
        out["x"][:, :, k] = c * sx + s * sy                     # Rz(-deg) applied to the vector
        out["y"][:, :, k] = -s * sx + c * sy
        out["z"][:, :, k] = map_coordinates(bz[:, :, k], coords, order=1, mode="nearest").reshape(X.shape)
    log("B-field rotated by {:+.3f} deg about z on its {}x{}x{} grid".format(deg, len(x), len(y), len(z)))
    return Field.from_arrays(grid={"x": x, "y": y, "z": z}, values=out, dim=3, units="m",
                             label=label or "{} (rotated {:+.2f} deg)".format(getattr(field, "label", "B-field"), deg))
