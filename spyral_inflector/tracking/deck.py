"""Inputs of a spiral-inflector tracking run: the RFQ particle file, the STEP electrode
assembly, the field maps, and the electrode naming conventions of the exported geometry.

The beam species, the RF and the start plane are machine settings (spyral_inflector.machine), not
constants of this module.

Electrode names follow the deck's STEP export: ``SI_Anode``/``SI_Cathode`` (spiral
electrodes), ``Housing``, ``Entrance_Aperture``, the grounded quad plates ``ent0``/``ext0``/
``ent4``/``ext4`` and the quadrupole poles ``D0..D3`` (quad 1) and ``D4..D7`` (quad 2).
"""
import csv
import json
import os
import pickle

import numpy as np

from py_electrodes.py_electrodes import PyElectrode, PyElectrodeAssembly
from PyPATools.field import Field
from PyPATools.species import IonSpecies

from .. import machine

CLIGHT = 299792458.0

MAX_QUADS = 4               # the generator names the poles of quad i D<4i>..D<4i+3>


def quad_names(i):
    """The four pole electrodes of quadrupole i (0-based): D<4i>, D<4i+1> at +V, D<4i+2>, D<4i+3> at -V."""
    return tuple("D{}".format(4 * i + j) for j in range(4))


QUADS = [quad_names(i) for i in range(MAX_QUADS)]
QUAD1, QUAD2 = QUADS[0], QUADS[1]
QUAD_SIGN = {name: (+1 if j < 2 else -1) for q in QUADS for j, name in enumerate(q)}
SPIRAL_NAMES = {"SI_Anode", "SI_Cathode"}
QUAD_NAMES = set(QUAD_SIGN)


def quad_index(name):
    """0-based quadrupole index of a pole electrode name, or None."""
    for i, q in enumerate(QUADS):
        if name in q:
            return i
    return None


# Every quadrupole is bounded by a grounded plate at each end, and the plates are SHARED: with
# shared_plates=True the exit plate of quad i is the entrance plate of quad i+1, and with
# last_exit_plate=False the LAST quad has no plate of its own -- the inflector's own
# Entrance_Aperture terminates it (geometry.py runs that quad up to one plate_gap before the
# entrance plate's face). A plate must therefore appear in the solve of BOTH quads it sits
# between, and Entrance_Aperture must appear in the last quad's solve as well as the spiral's;
# leaving it out leaves that quad's field unterminated.
#
# ent<4i>/ext<4i> are discs with a circular hole (hole_type="ellipse", a = b = aper_rad), so they
# are surfaces of revolution about z and a quad turning inside them turns its field with it.
# Entrance_Aperture is NOT -- it carries a rectangular slot -- so once it is in the solve, turning
# one quad by itself is only approximate; a RIGID rotation of the whole sub-assembly stays exact.
QUAD_PLATE_PREFIXES = ("ent", "ext")
ENTRANCE_APERTURE = "Entrance_Aperture"


def is_quad_plate(name):
    """True for a grounded plate that bounds a quadrupole: the circular ent<4i>/ext<4i> discs of
    the quad line, and the inflector entrance aperture, which is the last quad's exit plate."""
    return name == ENTRANCE_APERTURE or any(
        name.startswith(p) and name[len(p):].isdigit() for p in QUAD_PLATE_PREFIXES)


def isolated_quad_drop(names, index):
    """Electrode names to drop so that only quadrupole `index` and the grounded plates that bound
    the quad line are left -- the sub-assembly an isolated quad basis is solved in.

    What goes is the spiral electrodes, the housing and the poles of the other quadrupoles. Every
    bounding plate stays, including the shared ones and the entrance aperture, because those are
    what terminate the field; dropping them is not a small perturbation.
    """
    keep = set(quad_names(index))
    return sorted(n for n in names if n not in keep and not is_quad_plate(n))


# ---------------------------------------------------------------- voltages
def read_voltages(voltages):
    """Electrode voltages as {name: V} from a dict or the deck's voltages.csv
    (columns 'Electrode', 'Voltage (V)')."""
    if voltages is None:
        return {}
    if isinstance(voltages, dict):
        return dict(voltages)
    out = {}
    with open(voltages, newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["Electrode"]] = float(row["Voltage (V)"])
    return out


def write_voltages(path, volts):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Electrode", "Voltage (V)"])
        for name, volt in volts.items():
            w.writerow([name, "{:.3f}".format(volt)])


def set_quad_voltages(volts, *qs):
    """Quad i: D<4i>,D<4i+1> = +q_i, D<4i+2>,D<4i+3> = -q_i for every voltage given (in place,
    returned); quads beyond the given voltages are left as they are."""
    for i, q in enumerate(qs):
        for j, name in enumerate(quad_names(i)):
            volts[name] = (+1 if j < 2 else -1) * float(q)
    return volts


def set_spiral_voltage(volts, v):
    volts["SI_Anode"] = float(v)
    volts["SI_Cathode"] = -float(v)
    return volts


def load_state(state):
    """The design-orbit state written by the geometry build (dict or pickle path)."""
    if isinstance(state, dict):
        return state
    with open(state, "rb") as fh:
        return pickle.load(fh)


# ---------------------------------------------------------------- particles
from PyPATools.particles_src.particle_io import read_tracewin_dst as read_dst  # noqa: E402  (the reader lives in PyPATools)


def dst_rows(path, core=True):
    """A TraceWin .dst as rows of the text format (x(mm) x'(mrad) y(mm) y'(mrad) z(mm)
    z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss). TraceWin's phase is omega*t relative to
    the reference particle, so a late particle gets a negative z = -phase * beta*lambda / 2pi.
    core=True keeps the bunch core only: |phase| <= 180 deg and energy >= half the median
    (drops the unaccelerated stragglers), which is how the Bevatech core file was made."""
    b = read_dst(path)
    e, m = b["energy_MeV"], b["mass_MeV"]
    gamma = 1.0 + e / m
    beta = np.sqrt(1.0 - 1.0 / gamma ** 2)
    bl = beta * CLIGHT / (b["freq_MHz"] * 1e6)
    z = -b["phase_rad"] / (2.0 * np.pi) * bl
    p = gamma * beta
    dpp = p / np.median(p) - 1.0
    rows = np.column_stack([1e3 * b["x"], 1e3 * b["xp"], 1e3 * b["y"], 1e3 * b["yp"], 1e3 * z, 1e3 * dpp,
                            np.degrees(b["phase_rad"]), np.zeros(b["n"]), e, np.zeros(b["n"])])
    if core:
        keep = (np.abs(np.degrees(b["phase_rad"])) <= 180.0) & (e >= 0.5 * np.median(e))
        rows = rows[keep]
    return rows


def core_selection(rows):
    """The bunch core of text-format rows: |phase| <= 180 deg and energy at least half the
    median (drops the RFQ's unaccelerated stragglers)."""
    phase, ekin = rows[:, 6], rows[:, 8]
    return (np.abs(phase) <= 180.0) & (ekin >= 0.5 * np.median(ekin))


def particle_rows(path, core=True):
    """Unlost particles of an RFQ file as text-format rows: a .dst is read directly,
    anything else is the 10-column text file. core=True keeps the bunch core only (the
    same selection for both formats; a text file that is already the core is unchanged)."""
    if path.lower().endswith(".dst"):
        return dst_rows(path, core=core)
    data = np.loadtxt(path, skiprows=1)
    data = data[data[:, 9] == 0]
    # 2026-10-23: every particle is placed with its OWN beta*lambda, as dst_rows does. The text exports carry
    # z = -phase/360 * beta*lambda at the mean energy of ALL particles (tail included), which put a 67 keV core
    # particle 2.6 % too close to the start plane (WARP 10 mA V1.10 export: 1.9 deg rms of arrival-time error)
    gamma = 1.0 + data[:, 8] / IonSpecies(machine.species()).mass_mev
    data[:, 4] = -data[:, 6] / 360.0 * np.sqrt(1.0 - 1.0 / gamma ** 2) * CLIGHT / machine.rf_hz() * 1e3
    return data[core_selection(data)] if core else data


def load_particles(path, n, seed=20260905, z_start=None, species=None, core=True):
    """Read an RFQ output file (TraceWin .dst, or the text format x(mm) x'(mrad) y(mm)
    y'(mrad) z(mm) z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss) and draw n particles from
    the unlost ones without replacement; n at or above the file size takes every particle once.

    Returns r (n, 3) [m], v (n, 3) [m/s], the IonSpecies and the raw rows drawn. The
    bunch starts at z_start with the file's z as the longitudinal spread. z_start and
    species default to the machine settings."""
    z_start = machine.z_start() if z_start is None else z_start
    species = machine.species() if species is None else species
    data = particle_rows(path, core=core)
    rng = np.random.default_rng(seed)
    if n >= len(data):
        idx = rng.permutation(len(data))          # every particle exactly once (n is capped at the file size)
    else:
        idx = rng.choice(len(data), size=n, replace=False)
    data = data[idx]
    ion = IonSpecies(species)
    x, xp, y, yp, z = (data[:, i] * 1e-3 for i in range(5))      # m, rad
    ekin = data[:, 8]                                             # MeV
    gamma = 1.0 + ekin / ion.mass_mev
    speed = np.sqrt(1.0 - 1.0 / gamma ** 2) * CLIGHT
    vz = speed / np.sqrt(1.0 + xp ** 2 + yp ** 2)
    r = np.column_stack([x, y, z_start + z])
    v = np.column_stack([vz * xp, vz * yp, vz])
    return r, v, ion, data


def load_particles_with_tail(path, n, seed=20260905, z_start=None, species=None, rf_mhz=None, inject_core=False):
    """Every particle of an RFQ file: the core (|phase| <= 180 deg and energy at least half
    the median) placed as load_particles places it -- a spatial bunch at z_start with the
    phase as longitudinal spread -- and the unaccelerated tail as particles to be injected
    AT z_start at their own arrival time t = phase / (2 pi f_rf), with their own energy.

    inject_core=True puts EVERY particle at the plane z_start and gives the core the same
    arrival-time treatment as the tail (the convention of the inflector -> central-region
    hand-off): the file's phase is the arrival time, t = phase_deg / (360 f_rf), which is
    exactly the time the spatial placement encodes as z = -phase/(2 pi) * beta*lambda, and
    the clock is shifted so the earliest particle starts at t = 0. No particle is then born
    inside an electrode bore or downstream of a field it never traversed. The price is that
    the charge still upstream of the plane is missing from the space-charge solve during the
    injection ramp; see the TODO in tracking/bunch.py.

    Returns r, v, ion, rows, t_inject [s] (0 for the core unless inject_core) and is_tail.
    n below the file size draws n particles at random from the whole file. z_start, species
    and rf_mhz default to the machine settings."""
    z_start = machine.z_start() if z_start is None else z_start
    species = machine.species() if species is None else species
    rf_mhz = machine.rf_mhz() if rf_mhz is None else rf_mhz
    data = particle_rows(path, core=False)
    rng = np.random.default_rng(seed)
    if n < len(data):
        data = data[rng.choice(len(data), size=n, replace=False)]
    is_tail = ~core_selection(data)
    phase, ekin = data[:, 6], data[:, 8]
    ion = IonSpecies(species)
    x, xp, y, yp, z = (data[:, i] * 1e-3 for i in range(5))
    gamma = 1.0 + ekin / ion.mass_mev
    speed = np.sqrt(1.0 - 1.0 / gamma ** 2) * CLIGHT
    vz = speed / np.sqrt(1.0 + xp ** 2 + yp ** 2)
    v = np.column_stack([vz * xp, vz * yp, vz])
    t_all = phase / 360.0 / (rf_mhz * 1e6)
    # 2026-10-23 (Daniel): a WARP export is a ONE-PERIOD snapshot of a DC beam -- the crossings of every recorded period
    # folded into [-180, 180) deg -- in which the tail overlaps the main bunch, as the tail of the preceding bunches does in
    # reality. Such a file (no particle beyond +-180 deg) is placed as one spatial snapshot, tail included: every particle
    # at z = -phase/360 * its own beta*lambda (particle_rows), nothing injected later. TraceWin's tail, trailing up to
    # thousands of degrees, is still injected at the plane at its arrival time.
    snapshot = bool(np.all(np.abs(phase) <= 180.0))
    if inject_core:
        r = np.column_stack([x, y, np.full(len(x), float(z_start))])
        t_inject = t_all - t_all.min()            # the earliest particle starts the clock
    elif snapshot:
        r = np.column_stack([x, y, z_start + z])
        t_inject = np.zeros(len(x))
    else:
        r = np.column_stack([x, y, z_start + np.where(is_tail, 0.0, z)])
        t_inject = np.where(is_tail, t_all, 0.0)
    return r, v, ion, data, t_inject, is_tail


def orient_beam(r, v, phi_deg=0.0, swap_xy=False):
    """Exchange x and y (first) and/or rotate the beam about z by phi_deg."""
    r, v = np.array(r, dtype=float), np.array(v, dtype=float)
    if swap_xy:
        r[:, [0, 1]] = r[:, [1, 0]]
        v[:, [0, 1]] = v[:, [1, 0]]
    if phi_deg:
        c, s = np.cos(np.radians(phi_deg)), np.sin(np.radians(phi_deg))
        rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        r, v = r @ rot.T, v @ rot.T
    return r, v


def mean_energy_mev(path, core=True):
    """Mean kinetic energy of the unlost particles of an RFQ file (.dst or text) [MeV]: the design energy."""
    return float(np.mean(particle_rows(path, core=core)[:, 8]))


# ---------------------------------------------------------------- geometry
# CAD files of a STEP folder (Daniel, 2026-10-16): every export is now ".stp"; runs from before have ".step".
# Solids that are NOT electrodes carry the stem suffix "_negative" (the gap volume between the spiral
# electrodes, which fills the beam channel and must never enter a field solve or a collision check); a
# combined *assembly* file is skipped as before.
CAD_EXTENSIONS = (".stp", ".step")
NEGATIVE_SUFFIX = "_negative"


def is_electrode_cad(fn):
    """True for an electrode's CAD file NNN_<Name>.stp / .step: not an assembly, not a _negative solid."""
    stem, ext = os.path.splitext(os.path.basename(fn))
    return (ext.lower() in CAD_EXTENSIONS and "assembly" not in stem.lower()
            and not stem.lower().endswith(NEGATIVE_SUFFIX))


def electrode_cad_files(step_dir):
    """Sorted electrode CAD file names (not paths) of a STEP folder."""
    return sorted(f for f in os.listdir(step_dir) if is_electrode_cad(f))


def cad_name(fn):
    """'NNN_<Name>.stp' -> '<Name>'."""
    return os.path.splitext(os.path.basename(fn))[0].split("_", 1)[1]


def find_cad(step_dir, name):
    """Path of the electrode file NNN_<name>.stp (or .step) in step_dir, None if there is none."""
    for f in electrode_cad_files(step_dir):
        if cad_name(f) == name:
            return os.path.join(step_dir, f)
    return None


# 2026-11-17: per-electrode mesh refinement of a STEP folder, {electrode name: {"mesh_per_2pi", "mesh_min" [m]}}, written
# by a build with a grooved anode (geometry_point). STEP files carry no mesh sizes, so without it the groove would be
# meshed (BEM and collisions) at the flat size h. Folders without the file mesh exactly as before.
MESH_REFINE_NAME = "mesh_refine.json"


def load_step_assembly(step_dir, voltages=None, name="spiral inflector from STEP"):
    """PyElectrodeAssembly from a folder of NNN_<Name>.stp (or older .step) files, one electrode each
    (a combined *assembly* export and _negative solids in the same folder are skipped). An electrode
    named in the folder's mesh_refine.json carries that refinement to mesh_assembly()."""
    volts = read_voltages(voltages)
    refine = {}
    refine_fn = os.path.join(step_dir, MESH_REFINE_NAME)
    if os.path.isfile(refine_fn):
        with open(refine_fn) as fh:
            refine = json.load(fh)
    assembly = PyElectrodeAssembly(name)
    for fn in electrode_cad_files(step_dir):
        ename = cad_name(fn)
        e = PyElectrode(name=ename, voltage=volts.get(ename, 0.0))
        if e.generate_from_file(os.path.join(step_dir, fn)) != 0:
            raise RuntimeError("could not load {}".format(fn))
        if ename in refine:
            e.mesh_refine = {k: float(v) for k, v in refine[ename].items()}
        assembly.add_electrode(e)
    return assembly


def _generate_mesh_refined(e, h, refine):
    """e.generate_mesh() with gmsh's mesh sizes from curvature switched on (py_electrodes opens a fresh
    gmsh session per electrode, so the options are set just before its mesh.generate call)."""
    import gmsh
    original = gmsh.model.mesh.generate

    def generate(dim=3):
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", float(refine["mesh_per_2pi"]))
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(refine["mesh_min"]))
        return original(dim)

    gmsh.model.mesh.generate = generate
    try:
        return e.generate_mesh(brep_h=h) if h else e.generate_mesh()
    finally:
        gmsh.model.mesh.generate = original


def rotate_quads(assembly, *alphas_deg):
    """Rotate quad i about the z axis by alphas_deg[i] (45 deg = skew quad); quads beyond the
    given angles are not rotated."""
    for e in assembly.electrodes.values():
        i = quad_index(e.name)
        if i is not None and i < len(alphas_deg) and alphas_deg[i] != 0.0:
            e.set_rotation_angle_axis(float(np.radians(alphas_deg[i])), np.array([0.0, 0.0, 1.0]), absolute=False)
    return assembly


def shift_assembly(assembly, dz):
    """Translate every electrode along z by dz [m] (the whole inflector system moved axially).
    Must be called before the meshes are generated."""
    if dz:
        for e in assembly.electrodes.values():
            e.set_translation(np.array([0.0, 0.0, float(dz)]), absolute=False)
    return assembly


def rotate_assembly(assembly, deg):
    """Rigidly rotate every electrode about the z axis [deg] -- the whole inflector system
    turned in the magnet. Must be called before the meshes are generated."""
    if deg:
        for e in assembly.electrodes.values():
            e.set_rotation_angle_axis(float(np.radians(deg)), np.array([0.0, 0.0, 1.0]), absolute=False)
    return assembly


def mesh_assembly(assembly, h=None):
    """Force the surface meshes (collisions and the Poisson conductor lookup build from
    them lazily); returns the triangle count."""
    n_tri = 0
    for e in assembly.electrodes.values():
        refine = getattr(e, "mesh_refine", None)
        if refine:
            ok = _generate_mesh_refined(e, h, refine)
        else:
            ok = e.generate_mesh(brep_h=h) if h else e.generate_mesh()
        if ok != 0 or e._gmsh_msh is None:
            raise RuntimeError("failed to mesh {}".format(e.name))
        n_tri += len(e._gmsh_msh["elements"])
    return n_tri


def drop_electrodes(assembly, names):
    """Remove electrodes by name (diagnostic: they then neither stop particles nor
    bound the Poisson solve)."""
    drop = set(names)
    for uuid in [u for u, e in assembly.electrodes.items() if e.name in drop]:
        assembly.electrodes.pop(uuid)
    return assembly


# ---------------------------------------------------------------- fields
def with_fast_interpolator(field, label):
    """Re-wrap a Field on the numba backend (the pickled maps use scipy's, which is far
    slower for the ~1e8 evaluations of a run)."""
    grid, values = field.grid, field.grid_values
    if grid is None or values is None:
        return field
    return Field.from_arrays(grid=grid, values=values, label=label, dim=3, units="m", interpolator_backend="numba")


def load_bfield(path, frame="auto"):
    """The B-field map in the deck frame on the fast interpolator; a map delivered in the
    machine (Baseline) frame is mirrored in memory (tracking/frames.py)."""
    from .frames import to_deck_frame
    return with_fast_interpolator(to_deck_frame(Field.from_file(path), frame=frame), "B-field")


def basis_is_rotatable(basis_dir):
    """True when the quad bases of basis_dir were solved in isolation, so a quad angle is applied
    by rotating the field about z rather than by the cos/sin skew identity. Written by the basis
    driver as basis_meta.json {"rotatable": true}."""
    try:
        with open(os.path.join(basis_dir, "basis_meta.json")) as fh:
            return bool(json.load(fh).get("rotatable"))
    except (OSError, ValueError):
        return False


def superpose_basis(basis_dir, *quads, vscale=1.0, unit=3500.0, log=None):
    """Vacuum field summed on the grid from the basis solves of basis_dir:
    quads = (q1, alpha1, q2, alpha2[, q3, alpha3]), voltages scale linearly.

    How alpha is applied depends on how the bases were solved.

    ROTATABLE bases (basis_meta.json "rotatable": true, written after solving each quad in
    isolation with --isolate-quad) carry the whole field of one quad inside surfaces of
    revolution, so turning the quad by alpha turns the field by alpha -- exact up to the bilinear
    resampling. Only ef_itp_q<i>.pickle is needed.

    Otherwise the old identity is used: a quad rotated by alpha is cos(2 alpha) normal +
    sin(2 alpha) skew, which needs ef_itp_q<i>skew.pickle and is exact only for the pure
    quadrupole term of a field that also carries the fringe and the plate images.
    """
    if len(quads) % 2:
        raise ValueError("superpose_basis needs (voltage, angle) pairs, got {}".format(quads))
    pairs = [(float(quads[2 * i]), float(quads[2 * i + 1])) for i in range(len(quads) // 2)]
    rotatable = basis_is_rotatable(basis_dir)
    need = ["spiral"] + ["q{}".format(i + 1) for i in range(len(pairs))] + \
           ([] if rotatable else ["q{}skew".format(i + 1) for i, (_, a) in enumerate(pairs) if a != 0.0])
    basis = {k: Field.from_file(os.path.join(basis_dir, "ef_itp_{}.pickle".format(k))) for k in need}
    for k in need[1:]:
        for c in "xyz":
            if basis[k].grid_values[c].shape != basis["spiral"].grid_values[c].shape:
                raise RuntimeError("basis field {} is not on the spiral field's grid".format(k))
    if rotatable:
        from .frames import rotate_source_z
        quiet = log or (lambda m: None)
        for i, (_, a) in enumerate(pairs):
            if a != 0.0:
                key = "q{}".format(i + 1)
                basis[key] = rotate_source_z(basis[key], a, log=quiet)

    def quad(c, key, volt, alpha):
        if rotatable:
            return volt / unit * basis[key].grid_values[c]          # already turned, above
        a = np.radians(alpha)
        v = np.cos(2 * a) * basis[key].grid_values[c]
        if alpha != 0.0:
            v = v + np.sin(2 * a) * basis[key + "skew"].grid_values[c]
        return volt / unit * v

    values = {}
    for c in "xyz":
        v = vscale * basis["spiral"].grid_values[c]
        for i, (q, a) in enumerate(pairs):
            v = v + quad(c, "q{}".format(i + 1), q, a)
        values[c] = v
    label = "superposed spiral x{:.4f} ".format(vscale) + " ".join("q{}={:.0f}@{:.0f}".format(i + 1, q, a) for i, (q, a) in enumerate(pairs))
    return Field.from_arrays(grid=basis["spiral"].grid, values=values, dim=3, units="m", label=label)
