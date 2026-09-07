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
    aper_hole=0.0125, entrance_hole=None, plate_gap=0.001, plate_thickness=0.005,
    slot_width=0.015, slot_length=0.040, exit_opening=None,
)
KNOB_HELP = {
    "volt": "nominal spiral electrode voltage [V] (the optimizer scales it)", "gap": "electrode gap [m]",
    "tilt": "k' tilt [deg]", "dx": "electrode offset [m]", "sigma": "V-depth [m]", "aspect": "electrode width / gap",
    "gamma": "gammaAng, exit tilt [deg]", "angling": "anglingAng, inner face angling [deg]",
    "quad_bore": "quadrupole hyperbola vertex radius a = b [m]", "quad_z1": "quad 1 start z [m]", "quad_z2": "quad 2 start z [m]",
    "quad_len": "quad 1 length [m]", "quad_len2": "quad 2 length [m] (default: quad_len)",
    "shared_plates": "one grounded plate between the quads (quad 2 then starts at z1 + len1 + 2 plate_gap + plate_thickness)",
    "aper_hole": "quad aperture hole radius [m]", "entrance_hole": "hole radius of the first plate [m] (default: aper_hole)",
    "plate_gap": "axial gap between a quad's aperture plates and its pole ends [m]", "plate_thickness": "quad aperture plate thickness [m]",
    "slot_width": "inflector entrance slot width (gap direction) [m]", "slot_length": "inflector entrance slot length [m]",
    "exit_opening": "housing exit opening (across, along the tilted gap) [m] (default: the entrance slot)",
}


def _sanitize(name):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(name)).strip("_")


def build_geometry(out_dir, steps_dir, bfield, energy_mev, knobs=None, fix_truncations=(0.34, 0.77), fix_dz=None,
                   maxiter=8, res=0.005, h=0.005, log=None):
    """Generate the geometry, optimize the design particle, export the STEP files.

    energy_mev: design kinetic energy (the mean of the RFQ core file). knobs: overrides
    of DEFAULT_KNOBS. fix_truncations: entrance/exit truncation angles held fixed [deg].
    fix_dz: hold the axial shift of the whole system at this value [m] (then only the
    voltage is optimized). Returns a dict with steps_dir, voltages_csv, state_pickle,
    the optimized voltage/dz and the residuals; writes out_dir/summary.json.
    """
    from PyPATools.particles import ParticleDistribution
    from PyPATools.species import IonSpecies
    from ..spyral_inflector import SpiralInflector
    from .deck import SPECIES, write_voltages

    log = log or (lambda m: print(m, flush=True))
    k = dict(DEFAULT_KNOBS)
    k.update(knobs or {})
    if k["quad_len2"] is None:
        k["quad_len2"] = k["quad_len"]
    if k["shared_plates"]:
        k["quad_z2"] = k["quad_z1"] + k["quad_len"] + 2.0 * k["plate_gap"] + k["plate_thickness"]
    os.makedirs(out_dir, exist_ok=True)
    t_start = time.time()
    log("GEOMETRY: {}".format(k))
    log("  design energy {:.5f} MeV, truncations {} deg, dz {}".format(
        energy_mev, fix_truncations, "fixed at {:+.2f} mm".format(1e3 * fix_dz) if fix_dz is not None else "optimized"))

    ion = ParticleDistribution(species=IonSpecies(SPECIES))
    ion.set_mean_energy_z_mev(energy_mev)
    si = SpiralInflector(ion=ion, method="numerical", solver="bempp",
                         volt=k["volt"], gap=k["gap"], tilt=k["tilt"], dx=k["dx"], sigma=k["sigma"],
                         vee_shape="parabolic", ns=100, aspect_ratio=k["aspect"], rotation=0.0,
                         debug=False, gammaAng=k["gamma"], anglingAng=k["angling"])
    si.load_bfield(bfield=bfield)
    si.initialize()
    si.set_parameter(key="h", value=h)
    si.set_parameter(key="make_aperture", value=True)
    si.set_parameter(key="aperture_params", value={"thickness": 4e-3, "radius": 50e-3, "length": k["slot_length"], "width": k["slot_width"],
                                                   "top_distance": 5e-3, "bottom_distance": 10e-3, "hole_type": "rectangle", "voltage": 0.0,
                                                   **({"exit_width": k["exit_opening"][0], "exit_length": k["exit_opening"][1]}
                                                      if k["exit_opening"] else {})})
    si.set_parameter(key="make_housing", value=True)
    si.set_parameter(key="housing_params", value={"zmin": -0.12, "zmax": 0.03, "span": True, "gap": 6e-3, "thickness": 4e-3,
                                                  "voltage": 0.0, "experimental": True})
    si.set_parameter(key="make_quadrupoles", value=True)
    si.set_parameter(key="quadrupole_params", value={"a": k["quad_bore"], "b": k["quad_bore"], "radius": 0.04,
                                                     "z_starts": [k["quad_z1"], k["quad_z2"]], "lengths": [k["quad_len"], k["quad_len2"]],
                                                     "voltages": [3500, 3500], "aper_rad": 2.0 * k["aper_hole"],
                                                     "plate_gap": k["plate_gap"], "plate_thickness": k["plate_thickness"],
                                                     "shared_plates": bool(k["shared_plates"]),
                                                     "entrance_aper_rad": 2.0 * (k["entrance_hole"] if k["entrance_hole"] else k["aper_hole"])})
    si.generate_geometry()
    log("  geometry generated ({:.0f} s)".format(time.time() - t_start))

    fixed = {0: fix_truncations[0], 1: fix_truncations[1]}
    if fix_dz is not None:
        fixed[2] = fix_dz
    result = si.optimize_trajectory(maxiter=maxiter, solver="dfols", res=res,
                                    initial_guess=[fix_truncations[0], fix_truncations[1], 1.85e-3 if fix_dz is None else fix_dz, 0.97],
                                    fixed=fixed, bounds=((0.0, 15.0), (0.0, 15.0), (-15.0e-3, 15.0e-3), (0.85, 1.3)),
                                    exclude_quadrupoles=True)
    m = result["measurements"]
    optimizer = {
        "status": result["status"], "converged": result["converged"], "n_evaluations": result["n_evaluations"],
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
    if os.path.isdir(steps_dir):
        shutil.rmtree(steps_dir)
    os.makedirs(steps_dir)
    for i, electrode in enumerate(assembly.electrodes.values()):
        if electrode.export(os.path.join(steps_dir, "{:03d}_{}.step".format(i, names[i]))) != 0:
            raise RuntimeError("STEP export failed for {}".format(electrode.name))
    log("  exported {} STEP files -> {}".format(len(names), steps_dir))
    out = {"steps_dir": steps_dir, "voltages_csv": volt_fn, "state_pickle": state_fn, "knobs": k, "energy_mev": energy_mev,
           "fixed_truncations_deg": list(fix_truncations), "fix_dz_m": fix_dz, "res_m": res, "optimizer": optimizer,
           "voltage": result["voltage"], "dz_mm": 1e3 * result["dz"], "wall_s": time.time() - t_start}
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    del si
    return out
