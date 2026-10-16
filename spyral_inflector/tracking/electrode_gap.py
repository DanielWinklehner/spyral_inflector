"""The smallest distance between the two spiral electrodes, from their STEP solids.

Daniel, 2026-10-15: "The two spiral electrodes almost touch toward the end. That will spark 100%.
We need to have a minimum gap requirement." The design gap (the knob `gap`) is the spacing along
the beam; the electrodes' ends were never checked, and on the e071 / e016e shapes they converge to
1.9-2.4 mm just past the beam exit (23 kV across: ~12 kV/mm). This measures the exact solid-to-
solid distance with OpenCASCADE (gmsh), about a second per geometry, so every build can report it
and an optimizer can be held to a minimum.
"""
import glob
import os

import numpy as np

ANODE, CATHODE = "SI_Anode", "SI_Cathode"


def _volumes(gmsh, path):
    return gmsh.model.occ.importShapes(path)


def electrode_gap(steps_dir=None, anode_step=None, cathode_step=None):
    """{'min_mm', 'point_mm' (midpoint, deck frame), 'z_mm', 'r_mm', 'azimuth_deg'} of the closest
    approach between the anode and the cathode solids. Either a STEP folder (files named
    NNN_SI_Anode.stp / NNN_SI_Cathode.stp, as build_geometry writes them) or the two files."""
    import gmsh
    if steps_dir is not None:
        from .deck import find_cad
        anode_step = anode_step or find_cad(steps_dir, ANODE)
        cathode_step = cathode_step or find_cad(steps_dir, CATHODE)
        if anode_step is None or cathode_step is None:
            raise FileNotFoundError("no {} / {} CAD file in {}".format(ANODE, CATHODE, steps_dir))
    started = not gmsh.isInitialized()
    if started:
        gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    previous = gmsh.model.getCurrent() if not started else None
    gmsh.model.add("electrode_gap_check")
    try:
        va = _volumes(gmsh, anode_step)
        vc = _volumes(gmsh, cathode_step)
        gmsh.model.occ.synchronize()
        best = (np.inf, None, None)
        for da, ta in va:
            for dc, tc in vc:
                d, x1, y1, z1, x2, y2, z2 = gmsh.model.occ.getDistance(da, ta, dc, tc)
                if d < best[0]:
                    best = (d, np.array([x1, y1, z1]), np.array([x2, y2, z2]))
    finally:
        gmsh.model.remove()
        if started:
            gmsh.finalize()
        elif previous:
            try:
                gmsh.model.setCurrent(previous)
            except Exception:  # noqa: BLE001
                pass
    d, pa, pc = best
    mid = 0.5 * (pa + pc)
    return {"min_mm": 1e3 * float(d), "point_mm": (1e3 * mid).tolist(), "z_mm": 1e3 * float(mid[2]),
            "r_mm": 1e3 * float(np.hypot(mid[0], mid[1])), "azimuth_deg": float(np.degrees(np.arctan2(mid[1], mid[0])))}
