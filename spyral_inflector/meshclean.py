"""Clean a surface mesh before it goes to the BEM solver.

The electrode solids come out of gmsh with unmerged duplicate vertices at surface seams
and, where a truncated electrode end leaves a thin stub, sliver triangles (edges of
20 um, aspect ratios above 10,000). Both wreck the conditioning of the boundary-element
system: a pair of SI electrodes that should converge in a handful of GMRES iterations
took 400, and with the housing around them the solve hit the iteration cap with the
residual an order of magnitude above tolerance. Merging coincident vertices and dropping
degenerate triangles removes a negligible amount of surface and restores the conditioning.

Mesh layout as py_electrodes delivers it: verts (3, N), elems (3, M), domns (M,).
"""
import numpy as np


def clean_surface_mesh(mesh, merge_tol=1e-7, min_area=1e-8, max_aspect=100.0, min_edge=1e-4, log=None):
    """Return a cleaned copy of {"verts", "elems", "domns"} and a dict describing what changed.

    merge_tol: vertices closer than this [m] are one vertex. min_area [m^2], max_aspect
    (longest edge squared over area; an equilateral triangle has 2.3) and min_edge [m] drop
    degenerate triangles. On a 5 mm mesh a typical element has 10 mm2; the 0.4 mm vee edge
    of the electrodes meshes at ~1 mm2 and aspect ~25, so 0.01 mm2 / 100 / 0.1 mm keep every
    real feature and remove only needles and specks -- and it is the specks, elements
    thousands of times smaller than their neighbours, that stall GMRES, not the aspect
    ratio alone. Unreferenced vertices are removed and the indices compacted."""
    V = np.asarray(mesh["verts"], dtype=float)
    T = np.asarray(mesh["elems"], dtype=np.int64)
    D = np.asarray(mesh["domns"])
    transposed = V.shape[0] == 3 and V.shape[1] != 3
    if transposed:
        V, T = V.T, T.T
    n_v0, n_t0 = len(V), len(T)

    # 1. merge coincident vertices
    key = np.round(V / merge_tol).astype(np.int64)
    _, first, inverse = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    V = V[first]
    T = inverse[T]
    n_merged = n_v0 - len(V)

    # 2. drop collapsed and degenerate triangles
    ok = (T[:, 0] != T[:, 1]) & (T[:, 1] != T[:, 2]) & (T[:, 0] != T[:, 2])
    A, B, C = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(B - A, C - A), axis=1)
    emax2 = np.maximum(np.maximum(((B - A) ** 2).sum(1), ((C - B) ** 2).sum(1)), ((A - C) ** 2).sum(1))
    aspect = emax2 / np.maximum(area, 1e-300)
    n_collapsed = int((~ok).sum())
    emin2 = np.minimum(np.minimum(((B - A) ** 2).sum(1), ((C - B) ** 2).sum(1)), ((A - C) ** 2).sum(1))
    ok &= (area >= min_area) & (aspect <= max_aspect) & (emin2 >= min_edge ** 2)
    n_sliver = n_t0 - n_collapsed - int(ok.sum())
    area_dropped = float(area[~ok & (area > 0)].sum())
    T, D = T[ok], D[ok]

    # 3. compact the vertex list
    used, T = np.unique(T.ravel(), return_inverse=True)
    T = T.reshape(-1, 3)
    V = V[used]

    info = {"vertices_merged": int(n_merged), "triangles_collapsed": n_collapsed, "triangles_sliver": int(n_sliver),
            "area_dropped_mm2": 1e6 * area_dropped, "n_vertices": int(len(V)), "n_triangles": int(len(T)),
            "n_vertices_before": int(n_v0), "n_triangles_before": int(n_t0)}
    if log is not None and (n_merged or n_collapsed or n_sliver):
        log("   mesh clean: merged {} duplicate vertices, dropped {} collapsed and {} sliver triangles "
            "({:.4f} mm2 of surface); {} triangles remain".format(n_merged, n_collapsed, n_sliver, 1e6 * area_dropped, len(T)))
    out = {"verts": V.T if transposed else V, "elems": T.T if transposed else T, "domns": D}
    for k, v in mesh.items():
        if k not in out:
            out[k] = v
    return out, info
