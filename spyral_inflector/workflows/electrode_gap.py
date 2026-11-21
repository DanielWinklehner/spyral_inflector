"""Minimum distance between the two spiral electrodes (anode / cathode) of a geometry, and where
it is, from the STEP files (exact solid-to-solid distance from OpenCASCADE via gmsh), plus a
profile of the gap along the electrodes from a fine surface mesh. CPU only.

    python -m spyral_inflector.workflows.electrode_gap <steps folder> [<steps folder> ...] [--labels ...] [--voltage-kv 23.2] [--h 0.0008]

Frame and units are the STEP files' (deck frame, metres, for the run folders' steps/).
The profile bins the anode's surface nodes by their distance to the cathode surface along z.
"""
import argparse
import glob
import os

import numpy as np


def main(argv=None):
    """Minimum anode-cathode distance of STEP geometries (command line: see the module docstring)."""
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
    p.add_argument("steps", nargs="+")
    p.add_argument("--labels", nargs="+", default=None)
    p.add_argument("--voltage-kv", type=float, nargs="+", default=None, help="anode-cathode potential difference [kV] per geometry (for kV/mm)")
    p.add_argument("--h", type=float, default=0.0008, help="surface mesh size for the profile [m]")
    p.add_argument("--plot", default=None)
    a = p.parse_args(argv)
    labels = a.labels or [os.path.basename(os.path.dirname(os.path.abspath(s))) for s in a.steps]

    import gmsh  # noqa: E402
    from scipy.spatial import cKDTree  # noqa: E402

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.Verbosity", 1)
    results = []
    for k, (sdir, lab) in enumerate(zip(a.steps, labels)):
        fa = [f for f in _cad_files(sdir) if _cad_name(f) == "SI_Anode"][0]
        fc = [f for f in _cad_files(sdir) if _cad_name(f) == "SI_Cathode"][0]
        gmsh.model.add(lab)
        va = gmsh.model.occ.importShapes(fa)
        vc = gmsh.model.occ.importShapes(fc)
        gmsh.model.occ.synchronize()
        dmin, pa, pc = np.inf, None, None
        for da, ta in va:
            for dc, tc in vc:
                d, x1, y1, z1, x2, y2, z2 = gmsh.model.occ.getDistance(da, ta, dc, tc)
                if d < dmin:
                    dmin, pa, pc = d, np.array([x1, y1, z1]), np.array([x2, y2, z2])
        # profile: mesh both surfaces finely, distance of each anode node to the cathode nodes
        gmsh.option.setNumber("Mesh.MeshSizeMax", a.h)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.25 * a.h)
        gmsh.model.mesh.generate(2)
        def nodes_of(vols):
            pts = []
            for dim, tag in vols:
                for _, s in gmsh.model.getBoundary([(dim, tag)], oriented=False):
                    nt, xyz, _ = gmsh.model.mesh.getNodes(2, s, includeBoundary=True)
                    pts.append(np.asarray(xyz).reshape(-1, 3))
            return np.vstack(pts)
        A, C = nodes_of(va), nodes_of(vc)
        dA, _ = cKDTree(C).query(A)
        gmsh.model.remove()
        mid = 0.5 * (pa + pc)
        res = {"label": lab, "dmin": dmin, "at": mid, "r": np.hypot(mid[0], mid[1]), "az": np.degrees(np.arctan2(mid[1], mid[0])),
               "A": A, "dA": dA}
        results.append(res)
        V = a.voltage_kv[k] if a.voltage_kv and k < len(a.voltage_kv) else None
        print("{}: minimum anode-cathode distance {:.2f} mm at z {:+.1f} mm, r {:.1f} mm, azimuth {:+.1f} deg (deck frame){}".format(
            lab, 1e3 * dmin, 1e3 * mid[2], 1e3 * res["r"], res["az"],
            "; {:.1f} kV across it = {:.1f} kV/mm mean field".format(V, V / (1e3 * dmin)) if V else ""))
        zs = A[:, 2]
        edges = np.linspace(zs.min(), zs.max(), 13)
        line = []
        for z0, z1 in zip(edges[:-1], edges[1:]):
            m = (zs >= z0) & (zs < z1)
            if m.any():
                line.append("{:+.0f}..{:+.0f}: {:.1f}".format(1e3 * z0, 1e3 * z1, 1e3 * dA[m].min()))
        print("   smallest gap per z band [mm]: " + "; ".join(line))
    gmsh.finalize()

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 5))
        for res in results:
            zs, d = 1e3 * res["A"][:, 2], 1e3 * res["dA"]
            order = np.argsort(zs)
            bins = np.linspace(zs.min(), zs.max(), 60)
            idx = np.digitize(zs, bins)
            zb = [0.5 * (bins[i - 1] + bins[i]) for i in range(1, len(bins)) if np.any(idx == i)]
            db = [d[idx == i].min() for i in range(1, len(bins)) if np.any(idx == i)]
            ax.plot(zb, db, "-o", ms=3, label="{} (min {:.1f} mm)".format(res["label"], 1e3 * res["dmin"]))
        ax.set_xlabel("z along the inflector [mm, deck frame: entrance at -110, exit near 0]")
        ax.set_ylabel("smallest anode-cathode distance in the band [mm]")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=120)
        print("-> " + a.plot)


if __name__ == "__main__":
    main()
