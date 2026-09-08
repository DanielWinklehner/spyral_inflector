"""Inputs of a spiral-inflector tracking run: the RFQ particle file, the STEP electrode
assembly, the field maps, and the naming conventions of the HCHC-60 deck.

Electrode names follow the deck's STEP export: ``SI_Anode``/``SI_Cathode`` (spiral
electrodes), ``Housing``, ``Entrance_Aperture``, the grounded quad plates ``ent0``/``ext0``/
``ent4``/``ext4`` and the quadrupole poles ``D0..D3`` (quad 1) and ``D4..D7`` (quad 2).
"""
import csv
import os
import pickle

import numpy as np

from py_electrodes.py_electrodes import PyElectrode, PyElectrodeAssembly
from PyPATools.field import Field
from PyPATools.species import IonSpecies

Z_START = -0.275            # m, where the RFQ bunch starts (upstream of the first quad plate)
SPECIES = "H2_1+"
CLIGHT = 299792458.0
RF_FREQ_HZ = 32.8e6         # cyclotron RF: one bunch carries I / f_RF of charge

QUAD1 = ("D0", "D1", "D2", "D3")
QUAD2 = ("D4", "D5", "D6", "D7")
QUAD_SIGN = {"D0": +1, "D1": +1, "D2": -1, "D3": -1, "D4": +1, "D5": +1, "D6": -1, "D7": -1}
SPIRAL_NAMES = {"SI_Anode", "SI_Cathode"}
QUAD_NAMES = set(QUAD1) | set(QUAD2)


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


def set_quad_voltages(volts, q1, q2):
    """D0,D1 = +q1, D2,D3 = -q1, D4,D5 = +q2, D6,D7 = -q2 (in place, returned)."""
    for name, sgn in QUAD_SIGN.items():
        volts[name] = sgn * (q1 if name in QUAD1 else q2)
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


def particle_rows(path, core=True):
    """Unlost particles of an RFQ file as text-format rows: a .dst is read directly
    (core selection as above), anything else is the 10-column text file."""
    if path.lower().endswith(".dst"):
        return dst_rows(path, core=core)
    data = np.loadtxt(path, skiprows=1)
    return data[data[:, 9] == 0]


def load_particles(path, n, seed=20260905, z_start=Z_START, species=SPECIES, core=True):
    """Read an RFQ output file (TraceWin .dst, or the text format x(mm) x'(mrad) y(mm)
    y'(mrad) z(mm) z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss) and draw n particles from
    the unlost ones without replacement; n at or above the file size takes every particle once.

    Returns r (n, 3) [m], v (n, 3) [m/s], the IonSpecies and the raw rows drawn. The
    bunch starts at z_start with the file's z as the longitudinal spread."""
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
def load_step_assembly(step_dir, voltages=None, name="spiral inflector from STEP"):
    """PyElectrodeAssembly from a folder of NNN_<Name>.step files, one electrode each
    (a combined *assembly* export in the same folder is skipped)."""
    volts = read_voltages(voltages)
    assembly = PyElectrodeAssembly(name)
    for fn in sorted(os.listdir(step_dir)):
        if not fn.lower().endswith(".step") or "assembly" in fn.lower():
            continue
        ename = os.path.splitext(fn)[0].split("_", 1)[1]
        e = PyElectrode(name=ename, voltage=volts.get(ename, 0.0))
        if e.generate_from_file(os.path.join(step_dir, fn)) != 0:
            raise RuntimeError("could not load {}".format(fn))
        assembly.add_electrode(e)
    return assembly


def rotate_quads(assembly, alpha1_deg, alpha2_deg):
    """Rotate quad 1 (D0-D3) and quad 2 (D4-D7) about the z axis (45 deg = skew quad)."""
    for e in assembly.electrodes.values():
        if e.name in QUAD1 and alpha1_deg != 0.0:
            e.set_rotation_angle_axis(float(np.radians(alpha1_deg)), np.array([0.0, 0.0, 1.0]), absolute=False)
        elif e.name in QUAD2 and alpha2_deg != 0.0:
            e.set_rotation_angle_axis(float(np.radians(alpha2_deg)), np.array([0.0, 0.0, 1.0]), absolute=False)
    return assembly


def mesh_assembly(assembly, h=None):
    """Force the surface meshes (collisions and the Poisson conductor lookup build from
    them lazily); returns the triangle count."""
    n_tri = 0
    for e in assembly.electrodes.values():
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


def load_bfield(path):
    return with_fast_interpolator(Field.from_file(path), "B-field")


def superpose_basis(basis_dir, q1, alpha1, q2, alpha2, vscale=1.0, unit=3500.0):
    """Vacuum field summed on the grid from the basis solves ef_itp_spiral/q1/q1skew/q2/
    q2skew.pickle of basis_dir: a quad rotated by alpha is cos(2 alpha) normal +
    sin(2 alpha) skew (exact for the quadrupole term), voltages scale linearly."""
    need = ["spiral", "q1", "q2"] + (["q1skew"] if alpha1 != 0.0 else []) + (["q2skew"] if alpha2 != 0.0 else [])
    basis = {k: Field.from_file(os.path.join(basis_dir, "ef_itp_{}.pickle".format(k))) for k in need}
    for k in need[1:]:
        for c in "xyz":
            if basis[k].grid_values[c].shape != basis["spiral"].grid_values[c].shape:
                raise RuntimeError("basis field {} is not on the spiral field's grid".format(k))

    def quad(c, key, volt, alpha):
        a = np.radians(alpha)
        v = np.cos(2 * a) * basis[key].grid_values[c]
        if alpha != 0.0:
            v = v + np.sin(2 * a) * basis[key + "skew"].grid_values[c]
        return volt / unit * v

    values = {c: vscale * basis["spiral"].grid_values[c] + quad(c, "q1", q1, alpha1) + quad(c, "q2", q2, alpha2) for c in "xyz"}
    label = "superposed spiral x{:.4f} q1={:.0f}@{:.0f} q2={:.0f}@{:.0f}".format(vscale, q1, alpha1, q2, alpha2)
    return Field.from_arrays(grid=basis["spiral"].grid, values=values, dim=3, units="m", label=label)
