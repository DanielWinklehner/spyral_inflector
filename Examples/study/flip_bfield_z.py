"""Re-orient the HCHC-60 cyclotron field map for the spiral inflector.

The map from the magnet model has its median plane at z = 0 and the centre bore
running to z = +0.25 m, so the beam would travel in -z. The spiral inflector code
expects the opposite: the beam comes down the bore from negative z to the median
plane at z = 0. The map therefore has to be turned over.

The operation is a 180 degree rotation about the x-axis, (x, y, z) -> (x, -y, -z),
which for a vector field means

    B'_x(x, y, z) = +B_x(x, -y, -z)
    B'_y(x, y, z) = -B_y(x, -y, -z)
    B'_z(x, y, z) = -B_z(x, -y, -z)

This is a proper rotation, not a mirror, so it is physically just "mount the magnet
the other way up" and it preserves div B = 0. A bare mirror in z would NOT do: the
map's Bz is symmetric about z = 0, so mirroring alone leaves Bz positive, whereas
the inflector expects Bz negative near the median plane (the commented-out uniform
field in the generation script is -0.901 T).

Result: z runs -0.25 .. +0.05 m, and on-axis Bz is -0.891 T at z = 0.

Note the flipped map still stops at z = -0.25 m, while the bunch starts at
z = -0.275 m and the first quadrupole begins at z = -0.27 m. Those sit outside the
map, where the interpolator returns 0.

Run from anywhere; paths below are absolute.
"""

import os

import h5py
import numpy as np

from PyPATools.field import Field

SRC = (r"D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\Shared_With_Jarrett_and_Sam"
       r"\Fields\HCHC-60"
       r"\20260903_HCHC-60_3D_CenterBore_Res1mm_z-5to25cm_xy20cm_IBA_BH_HiCoil_CutPole_Tesla_m.h5part")

DST = (r"D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron"
       r"\Spiral_inflector\Fields\HCHC-60_CentralBField_zflipped.pickle")

# (x, y, z) component sign after the rotation
SIGNS = {"x": +1.0, "y": -1.0, "z": -1.0}
DATASETS = {"x": "0", "y": "1", "z": "2"}


def main():
    print("Reading {}".format(os.path.basename(SRC)))

    values = {}

    with h5py.File(SRC, "r") as fh:
        group = fh["Step#0/Block/Hfield"]

        origin = np.asarray(group.attrs["__Origin__"], dtype=float)
        spacing = np.asarray(group.attrs["__Spacing__"], dtype=float)
        nz, ny, nx = group["0"].shape

        for comp, dset in DATASETS.items():
            # File order is (nz, ny, nx). Reversing the z and y axes evaluates the
            # source at (x, -y, -z); the transpose then puts it in the (nx, ny, nz)
            # order the interpolators want. ascontiguousarray drops the reversed
            # view so the original array is not kept alive by the result.
            raw = np.asarray(group[dset])
            flipped = np.transpose(raw[::-1, ::-1, :], (2, 1, 0))
            values[comp] = np.ascontiguousarray(SIGNS[comp] * flipped)
            del raw, flipped
            print("  {} component done".format(comp))

    x = origin[0] + np.arange(nx) * spacing[0]
    y = origin[1] + np.arange(ny) * spacing[1]
    z = origin[2] + np.arange(nz) * spacing[2]

    # New axes: y' = -y reversed, z' = -z reversed. y is symmetric so it is unchanged.
    y_new = -y[::-1]
    z_new = -z[::-1]

    assert np.allclose(y_new, y), "y axis is not symmetric; the rotation would move it"

    print("")
    print("  x: {:+.3f} .. {:+.3f} m".format(x[0], x[-1]))
    print("  y: {:+.3f} .. {:+.3f} m".format(y_new[0], y_new[-1]))
    print("  z: {:+.3f} .. {:+.3f} m   (was {:+.3f} .. {:+.3f})".format(
        z_new[0], z_new[-1], z[0], z[-1]))

    field = Field.from_arrays(
        grid={"x": x, "y": y_new, "z": z_new},
        values=values,
        label="HCHC-60 central B-field, rotated 180 deg about x",
        dim=3,
        units="m",
        # scipy rather than 'auto': these interpolators are pickled, and a GPU
        # backend would put cupy arrays in the file.
        interpolator_backend="scipy",
    )

    print("")
    print("  on-axis check (x = y = 0):")
    probe = np.array([[0.0, 0.0, zz] for zz in (-0.25, -0.20, -0.10, -0.05, 0.0, 0.05)])
    for point, value in zip(probe, field(probe)):
        print("    z = {:+.3f} m -> Bz = {:+.5f} T".format(point[2], value[2]))

    os.makedirs(os.path.dirname(DST), exist_ok=True)
    field.save(DST)

    print("")
    print("  wrote {} ({:.2f} GB)".format(DST, os.path.getsize(DST) / 1e9))


if __name__ == "__main__":
    main()
