"""Build a spiral inflector with quadrupoles from its knobs, let the trajectory optimizer
set the spiral voltage (and the axial shift dz, unless fixed) for the design particle,
and export the electrodes as STEP files with the design-orbit state the tracking needs.

This is the geometry stage of the HCHC-60 study (GeometryPoint.py); the knob names and
defaults are the deck's."""
import csv
import json
import os
import pickle
import re
import shutil
import time

import numpy as np

DEFAULT_KNOBS = dict(
    volt=12000.0, gap=0.019, tilt=31.0, dx=0.01, sigma=0.0022, aspect=2.4, gamma=5.0, angling=11.0,
    quad_bore=0.013, quad_z1=-0.27, quad_z2=-0.19, quad_len=0.045, quad_len2=None, shared_plates=False,
    n_quads=2, quad_z3=None, quad_len3=None,
    aper_hole=0.0125, entrance_hole=None, entrance_rect=None, exit_rects=None, plate_gap=0.001, plate_thickness=0.005,
    slot_width=0.015, slot_length=0.040, exit_opening=None,
    housing_gap=0.006, housing_thickness=0.004, top_distance=0.005, bottom_distance=0.010, rotation=0.0,
    q2_exit_plate=True,
)
KNOB_HELP = {
    "volt": "nominal spiral electrode voltage [V] (the optimizer scales it)", "gap": "electrode gap [m]",
    "tilt": "k' tilt [deg]", "dx": "electrode offset [m]", "sigma": "V-depth [m]", "aspect": "electrode width / gap",
    "gamma": "gammaAng, exit tilt [deg]", "angling": "anglingAng, inner face angling [deg]",
    "quad_bore": "quadrupole hyperbola vertex radius a = b [m]", "quad_z1": "quad 1 start z [m]", "quad_z2": "quad 2 start z [m]",
    "quad_len": "quad 1 length [m]", "quad_len2": "quad 2 length [m] (default: quad_len)",
    "n_quads": "2 (doublet) or 3 (triplet)", "quad_z3": "quad 3 start z [m] (shared_plates: derived)", "quad_len3": "quad 3 length [m] (default: quad_len)",
    "shared_plates": "one grounded plate between the quads (quad 2 then starts at z1 + len1 + 2 plate_gap + plate_thickness)",
    "aper_hole": "quad aperture hole radius [m]", "entrance_hole": "hole radius of the first plate [m] (default: aper_hole)",
    "entrance_rect": "[half-width u, half-height w [m], angle [deg]]: a RECTANGULAR first-plate hole instead of the round one; "
                     "u is turned by angle from the deck x axis (the plate turns with the assembly, not with quad 1)",
    "exit_rects": "per quad, None or [half-width u, half-height w [m], angle [deg]]: a RECTANGULAR hole in the plate AFTER "
                  "that quad (index 0 = the plate between quads 1 and 2) instead of the round aper_hole",
    "plate_gap": "axial gap between a quad's aperture plates and its pole ends [m]", "plate_thickness": "quad aperture plate thickness [m]",
    "slot_width": "inflector entrance slot width (gap direction) [m]", "slot_length": "inflector entrance slot length [m]",
    "exit_opening": "housing exit opening (across, along the tilted gap) [m] (default: the entrance slot)",
    "housing_gap": "clearance between the electrodes' convex hull and the housing's inner wall [m]",
    "housing_thickness": "housing wall thickness [m]", "top_distance": "entrance plate height above the electrodes [m]",
    "bottom_distance": "distance of the housing exit opening past the electrode end [m]",
    "rotation": "rotation of the whole inflector (electrodes, housing, apertures) about the axis [deg]; rotate the quads and the beam by the same angle",
    "q2_exit_plate": "False: no exit plate on quad 2, the inflector entrance aperture terminates it (set quad_len2 to reach it)",
}


def _sanitize(name):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(name)).strip("_")


def rmtree_retry(path, log=None, tries=30, wait=10.0):
    """shutil.rmtree that survives a folder held open for a moment (Dropbox indexing a just-written
    folder on the deck gives WinError 32 on rmdir): retries every `wait` s up to `tries` times."""
    for attempt in range(tries):
        if not os.path.isdir(path):
            return
        try:
            shutil.rmtree(path)
            return
        except PermissionError as exc:
            if attempt == tries - 1:
                raise
            if log:
                log("  folder busy ({}); retrying in {:.0f} s".format(exc, wait))
            time.sleep(wait)


def build_geometry(out_dir, steps_dir, bfield, energy_mev, knobs=None, fix_truncations=(0.34, 0.77), fix_dz=None,
                   maxiter=8, res=0.005, h=0.005, log=None):
    """Generate the geometry, optimize the design particle, export the STEP files.

    energy_mev: design kinetic energy (the mean of the RFQ core file). knobs: overrides
    of DEFAULT_KNOBS. fix_truncations: entrance/exit truncation angles held fixed [deg];
    None frees that truncation as a knob of the optimizer (e.g. (0.34, None)).
    fix_dz: hold the axial shift of the whole system at this value [m] (then only the
    voltage is optimized). Returns a dict with steps_dir, voltages_csv, state_pickle,
    the optimized voltage/dz and the residuals; writes out_dir/summary.json.
    """
    from PyPATools.particles import ParticleDistribution
    from PyPATools.species import IonSpecies
    from ..spyral_inflector import SpiralInflector
    from .deck import SPECIES, write_voltages
    from .export import HOUSING_OUTLINE_NAME, write_housing_outline

    log = log or (lambda m: print(m, flush=True))
    k = dict(DEFAULT_KNOBS)
    k.update(knobs or {})
    if k["quad_len2"] is None:
        k["quad_len2"] = k["quad_len"]
    if k.get("quad_len3") is None:
        k["quad_len3"] = k["quad_len"]
    n_quads = int(k.get("n_quads", 2) or 2)
    if k["shared_plates"]:
        k["quad_z2"] = k["quad_z1"] + k["quad_len"] + 2.0 * k["plate_gap"] + k["plate_thickness"]
        if n_quads >= 3:
            k["quad_z3"] = k["quad_z2"] + k["quad_len2"] + 2.0 * k["plate_gap"] + k["plate_thickness"]
    z_starts = [k["quad_z1"], k["quad_z2"]] + ([k["quad_z3"]] if n_quads >= 3 else [])
    lengths = [k["quad_len"], k["quad_len2"]] + ([k["quad_len3"]] if n_quads >= 3 else [])
    os.makedirs(out_dir, exist_ok=True)
    t_start = time.time()
    log("GEOMETRY: {}".format(k))
    log("  design energy {:.5f} MeV, truncations {} deg, dz {}".format(
        energy_mev, fix_truncations, "fixed at {:+.2f} mm".format(1e3 * fix_dz) if fix_dz is not None else "optimized"))

    ion = ParticleDistribution(species=IonSpecies(SPECIES))
    ion.set_mean_energy_z_mev(energy_mev)
    si = SpiralInflector(ion=ion, method="numerical", solver="bempp",
                         volt=k["volt"], gap=k["gap"], tilt=k["tilt"], dx=k["dx"], sigma=k["sigma"],
                         vee_shape="parabolic", ns=100, aspect_ratio=k["aspect"], rotation=k["rotation"],
                         debug=False, gammaAng=k["gamma"], anglingAng=k["angling"])
    si.load_bfield(bfield=bfield)
    si.initialize()
    si.set_parameter(key="h", value=h)
    si.set_parameter(key="make_aperture", value=True)
    si.set_parameter(key="aperture_params", value={"thickness": 4e-3, "radius": 47e-3, "length": k["slot_length"], "width": k["slot_width"],
                                                   "top_distance": k["top_distance"], "bottom_distance": k["bottom_distance"],
                                                   "hole_type": "rectangle", "voltage": 0.0,
                                                   **({"exit_width": k["exit_opening"][0], "exit_length": k["exit_opening"][1]}
                                                      if k["exit_opening"] else {})})
    si.set_parameter(key="make_housing", value=True)
    si.set_parameter(key="housing_params", value={"zmin": -0.12, "zmax": 0.03, "span": True, "gap": k["housing_gap"],
                                                  "thickness": k["housing_thickness"], "voltage": 0.0, "experimental": True})
    si.set_parameter(key="make_quadrupoles", value=True)
    si.set_parameter(key="quadrupole_params", value={"a": k["quad_bore"], "b": k["quad_bore"], "radius": 0.04,
                                                     "z_starts": z_starts, "lengths": lengths,
                                                     "voltages": [3500] * len(z_starts), "aper_rad": 2.0 * k["aper_hole"],
                                                     "plate_gap": k["plate_gap"], "plate_thickness": k["plate_thickness"],
                                                     "shared_plates": bool(k["shared_plates"]),
                                                     "last_exit_plate": bool(k["q2_exit_plate"]),
                                                     "entrance_aper_rad": 2.0 * (k["entrance_hole"] if k["entrance_hole"] else k["aper_hole"]),
                                                     **({"entrance_hole_rect": (2.0 * float(k["entrance_rect"][0]), 2.0 * float(k["entrance_rect"][1]),
                                                                                float(np.radians(float(k["entrance_rect"][2]))))}
                                                        if k.get("entrance_rect") else {}),
                                                     **({"exit_hole_rects": [(2.0 * float(x[0]), 2.0 * float(x[1]), float(np.radians(float(x[2]))))
                                                                             if x else None for x in k["exit_rects"]]}
                                                        if k.get("exit_rects") else {})})
    si.generate_geometry()
    log("  geometry generated ({:.0f} s)".format(time.time() - t_start))

    # a truncation given as None is a free knob of the optimizer (starting from the default)
    t_ent = 0.34 if fix_truncations[0] is None else float(fix_truncations[0])
    t_exit = 0.77 if fix_truncations[1] is None else float(fix_truncations[1])
    fixed = {i: v for i, v in ((0, fix_truncations[0]), (1, fix_truncations[1])) if v is not None}
    if fix_dz is not None:
        fixed[2] = fix_dz
    result = si.optimize_trajectory(maxiter=maxiter, solver="dfols", res=res,
                                    initial_guess=[t_ent, t_exit, 1.85e-3 if fix_dz is None else fix_dz, 0.97],
                                    fixed=fixed, bounds=((0.0, 15.0), (0.0, 15.0), (-15.0e-3, 15.0e-3), (0.85, 1.3)),
                                    exclude_quadrupoles=True)
    m = result["measurements"]
    optimizer = {
        "status": result["status"], "converged": result["converged"], "n_evaluations": result["n_evaluations"],
        "truncations_deg": [float(result["db_entrance"]), float(result["db_exit"])],
        "free_truncations": [fix_truncations[0] is None, fix_truncations[1] is None],
        "dz_mm": 1e3 * result["dz"], "volt_scale": result["volt_scale"], "voltage_V": result["voltage"],
        "residual_final": {"angle_deg": float(result["residual_final"][0]), "centering_mm": 1e3 * float(result["residual_final"][1]),
                           "z_offset_mm": 1e3 * float(result["residual_final"][2]), "width_mm": 1e3 * float(result["residual_final"][3])},
        "clearance_exit_mm": {"anode": 1e3 * m["clearance_anode_exit"], "cathode": 1e3 * m["clearance_cathode_exit"]},
        "min_clearance_mm": 1e3 * m["min_clearance"], "exit_point_mm": (1e3 * m["exit_point"]).tolist()}
    log("  optimizer: voltage {:.1f} V, dz {:+.2f} mm, residuals angle {:+.3f} deg, z {:+.2f} mm ({}, {:.0f} s)".format(
        result["voltage"], 1e3 * result["dz"], result["residual_final"][0], 1e3 * result["residual_final"][2], result["status"],
        time.time() - t_start))

    # export: one STEP file per electrode, the voltages and the design-orbit state
    assembly = si.numerical_variables["objects"]
    names, used = {}, set()
    for i, electrode in enumerate(assembly.electrodes.values()):
        name = _sanitize(electrode.name) or "electrode_{:03d}".format(i)
        unique, j = name, 2
        while unique.lower() in used:
            unique = "{}_{}".format(name, j)
            j += 1
        used.add(unique.lower())
        names[i] = unique
    electrode_voltages = {names[i]: float(e.voltage) for i, e in enumerate(assembly.electrodes.values())}
    mesh = si.numerical_variables["full mesh"]
    ctx = result["context"]
    state = {"trj_design": si.analytic_variables["trj_design"] + result["shift"],
             "v_design": si.analytic_variables["v_design"],
             "trj_design_full": si.analytic_variables["trj_design_full"] + result["shift"],
             "v_design_full": si.analytic_variables["v_design_full"],
             "b_lim_deg": result["b_lim_deg"], "shift": result["shift"], "shift_lab": result["shift_lab"],
             "volt_scale": result["volt_scale"], "voltage": result["voltage"], "electrode_voltages": electrode_voltages,
             "z_start": ctx["z_start"], "v0": ctx["v0"], "dt": ctx["dt"], "nsteps": ctx["nsteps"],
             "track_r": m["r"], "track_v": m["v"],
             "mesh": {"verts": np.asarray(mesh["verts"]), "elems": np.asarray(mesh["elems"]), "domns": np.asarray(mesh["domns"]),
                      "names": {int(e.bempp_domain): names[i] for i, e in enumerate(assembly.electrodes.values())}},
             "residual_final": np.asarray(result["residual_final"]), "res_m": res, "knobs": k, "energy_mev": energy_mev}
    state_fn = os.path.join(out_dir, "state.pickle")
    with open(state_fn, "wb") as fh:
        pickle.dump(state, fh)
    volt_fn = os.path.join(out_dir, "voltages.csv")
    write_voltages(volt_fn, electrode_voltages)
    rmtree_retry(steps_dir, log=log)
    os.makedirs(steps_dir)
    for i, electrode in enumerate(assembly.electrodes.values()):
        if electrode.export(os.path.join(steps_dir, "{:03d}_{}.stp".format(i, names[i]))) != 0:
            raise RuntimeError("STEP export failed for {}".format(electrode.name))
    log("  exported {} STEP files -> {}".format(len(names), steps_dir))
    # the negative space between the two spiral electrodes as an extra CAD object (Daniel, 2026-10-16); its
    # "_negative" suffix keeps every STEP-folder loader from taking it for an electrode. A failure here must not
    # stop a build.
    try:
        from ..geometry import generate_gap_volume
        from .export import GAP_VOLUME_NAME
        import gmsh as _gmsh
        _gap_fn = os.path.join(steps_dir, GAP_VOLUME_NAME)
        if generate_gap_volume(si).export(_gap_fn) != 0:
            raise RuntimeError("export failed")
        _gmsh.initialize()
        try:
            _gmsh.option.setNumber("General.Terminal", 0)
            _gmsh.model.occ.importShapes(_gap_fn)
            _gmsh.model.occ.synchronize()
            _vol = sum(_gmsh.model.occ.getMass(3, t) for _d, t in _gmsh.model.getEntities(3))
        finally:
            _gmsh.finalize()
        log("  gap volume between the spiral electrodes: {:.1f} cm3 -> {}".format(1e6 * _vol, GAP_VOLUME_NAME))
    except Exception as exc:  # noqa: BLE001
        log("  gap volume export failed: {}".format(exc))
    # the closest approach of the two spiral electrodes (Daniel, 2026-10-15: they nearly touched at the
    # exit end on the e071 / e016e shapes -- a guaranteed spark); every build reports it
    try:
        from .electrode_gap import electrode_gap
        gap_info = electrode_gap(steps_dir)
        v_ac = abs(electrode_voltages.get("SI_Anode", 0.0) - electrode_voltages.get("SI_Cathode", 0.0))
        log("  electrode gap: anode-cathode minimum {:.2f} mm at z {:+.1f} mm, r {:.1f} mm, azimuth {:+.1f} deg "
            "({:.1f} kV across it: {:.1f} kV/mm)".format(gap_info["min_mm"], gap_info["z_mm"], gap_info["r_mm"],
                                                        gap_info["azimuth_deg"], 1e-3 * v_ac,
                                                        1e-3 * v_ac / max(gap_info["min_mm"], 1e-6)))
    except Exception as exc:  # noqa: BLE001
        log("  electrode gap check failed: {}".format(exc))
        gap_info = None
    # The central region wants the housing as a shape rather than as this solid, so write what
    # the solid was made from: the polygon the wall is extruded from and the four corners of the
    # exit opening. Same frame as the STEP files just written - the electrode's own rotation and
    # translation (the centering shift) applied exactly as PyElectrode.export applies them.
    _hs = assembly.get_electrode_by_name("Housing")
    _ol = getattr(_hs, "outline", None) if _hs is not None else None
    if _ol:
        _xy = _hs._transformation.apply_to_points(
            np.column_stack([_ol["outer"], np.zeros(len(_ol["outer"]))]))[:, :2]
        _ap = _hs._transformation.apply_to_points(np.array(_ol["aperture"], dtype=float))
        write_housing_outline(os.path.join(steps_dir, HOUSING_OUTLINE_NAME), _xy, _ap,
                              curves=_ol.get("outer_curves"),
                              comments=("geometry: {}".format(out_dir),
                                        "the housing solid is this polygon extruded over z = {:+.4f} .. {:+.4f} m "
                                        "with a {:.1f} mm wall inwards and {:.1f} mm clearance to the electrodes."
                                        .format(_ol["z_min"], _ol["z_max"], 1e3 * _ol["thickness"], 1e3 * _ol["gap"]),
                                        "the exit opening is {:.1f} x {:.1f} mm (along the gap x across it)."
                                        .format(1e3 * _ol["aperture_size"][0], 1e3 * _ol["aperture_size"][1])))
        log("  housing outline: {} outer vertices + 4 aperture corners -> {}".format(
            len(_xy), os.path.join(steps_dir, HOUSING_OUTLINE_NAME)))
    # surface field of the last BEM solve (the optimized voltage), per electrode
    from .bem_reload import surface_field_stats
    try:
        sf = surface_field_stats(si.numerical_variables["n_fun_coeff"], mesh, {int(e.bempp_domain): names[i] for i, e in enumerate(assembly.electrodes.values())})
        for name, st in sorted(sf.items(), key=lambda kv: -kv[1]["p99_kv_cm"])[:4]:
            log("  surface field {:<20s} p99 {:.1f} kV/cm, top-1%-area mean {:.1f}, max {:.1f} (mesh h {:.0f} mm)".format(
                name, st["p99_kv_cm"], st["top1pct_mean_kv_cm"], st["max_kv_cm"], 1e3 * h))
    except Exception as exc:  # noqa: BLE001
        log("  surface field statistics failed: {}".format(exc))
        sf = {}
    out = {"steps_dir": steps_dir, "voltages_csv": volt_fn, "state_pickle": state_fn, "knobs": k, "energy_mev": energy_mev,
           "fixed_truncations_deg": list(fix_truncations), "fix_dz_m": fix_dz, "res_m": res, "optimizer": optimizer,
           "voltage": result["voltage"], "dz_mm": 1e3 * result["dz"], "surface_field": sf, "surface_mesh_h_m": h,
           "electrode_gap": gap_info, "wall_s": time.time() - t_start}
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    del si
    return out
