"""Workflow configs: one TOML file per design (and one per machine), read by the chain and the CLI.

A design config names everything a levelling chain needs; disk locations are written with the roots of the per-machine
paths file (``${work}/Results``, see ``workflows.paths``), never with drive letters::

    # design.toml
    include = ["machine.toml"]          # other configs merged in first, relative to this file; later wins

    [machine]                           # see spyral_inflector.machine (usually in the included machine.toml)
    species = "H2_1+"
    rf_mhz = 32.8
    z_start_m = -0.270
    inject_core = true
    tag = "M1"                          # optional: exported file names
    name = "Machine 1"                  # optional: headers of exported files

    [run]
    results = "${work}/Results"         # results root; runs go to <results>/<family>/<prefix><pass>
    family = "triplet"
    prefix = "d1o_c"
    source = "${work}/Results/triplet/d1"            # run whose levelled setting starts pass 1
    best_knobs = "${work}/Results/triplet/d1/prof/best_knobs.json"   # optional: knobs of an optimizer winner
    quads_json = "${work}/Results/triplet/d1/prof/best_quads.json"   # quad voltages for the level step
    # vscale = 1.05                     # optional (default: best_knobs' _vscale, else the source's levelling)
    # quad_rotation_deg = [18.5, 22.7]  # optional basis turns (default: the source's levelling)
    # r0_deg = -97.0                    # optional first-pass rotation guess
    max_passes = 4
    dr_tol_deg = 0.2

    [inputs]
    bfield = "${work}/Fields/map.pickle"             # the UNROTATED B map
    particles = "${inputs}/Particles/beam.txt"
    current_ma = 10.0

    [rotation]
    rule = "corner"                     # offset | parallel | corner
    gap_azimuth_deg = 26
    corner = "outer"
    corner_offset_mm = 3.25
    half_gap_m = 0.00325

    [reports]
    reference_run = "${work}/Results/triplet/d1"     # optional (default: [run] source)
    reference_label = "d1 sc"                         # optional (default: "<reference folder> sc")
    gap_voltage_kv = 24.8
    # mass_kg = 3.3462e-27              # optional: ion mass of the report tools (default: the species' mass)
    render = true

    [steps]                             # optional extra arguments appended to one step's command line
    rotate = []
    level = []
    level_final = []
    export = []

Numbers reach the step command lines as repr() of the TOML value, which reads back as the same float.
"""
import copy
import json
import os
import tomllib

from .paths import Paths, load_paths


class ConfigError(ValueError):
    pass


SECTIONS = ("machine", "run", "inputs", "rotation", "reports", "steps")
STEP_NAMES = ("rotate", "level", "level_final", "export", "export_final", "gap", "report", "phase", "render")
REQUIRED = {
    "machine": ("species", "rf_mhz", "z_start_m"),
    "run": ("results", "prefix", "source", "quads_json"),
    "inputs": ("bfield", "particles", "current_ma"),
    "rotation": ("rule", "gap_azimuth_deg", "half_gap_m"),
}
DEFAULTS = {
    "run": {"family": "triplet", "max_passes": 4, "dr_tol_deg": 0.2},
    "rotation": {"corner": "outer"},
    "reports": {"render": True},
}
PATH_KEYS = {
    "run": ("results", "source", "best_knobs", "quads_json", "logs"),
    "inputs": ("bfield", "particles"),
    "reports": ("reference_run",),
}


def _merge(base, over):
    out = copy.deepcopy(base)
    for key, value in over.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def read_toml(path, _seen=None):
    """One config file with its includes merged in (includes first, the file itself last)."""
    path = os.path.abspath(path)
    seen = set(_seen or ())
    if path in seen:
        raise ConfigError("include loop at {}".format(path))
    seen.add(path)
    try:
        with open(path, "rb") as fh:
            data = tomllib.load(fh)
    except FileNotFoundError:
        raise ConfigError("config {} does not exist".format(path)) from None
    except tomllib.TOMLDecodeError as err:
        raise ConfigError("{}: {}".format(path, err)) from None
    includes = data.pop("include", [])
    if isinstance(includes, str):
        includes = [includes]
    merged = {}
    for inc in includes:
        inc_path = inc if os.path.isabs(inc) else os.path.join(os.path.dirname(path), inc)
        merged = _merge(merged, read_toml(inc_path, seen))
    return _merge(merged, data)


class Config:
    """A design config: the merged tables with every ${root} resolved, plus where it came from."""

    def __init__(self, data, source, paths):
        self.data = data
        self.source = source
        self.paths = paths

    def __getitem__(self, section):
        return self.data[section]

    def get(self, section, key, default=None):
        return self.data.get(section, {}).get(key, default)

    def step_args(self, step):
        return list(self.data.get("steps", {}).get(step, []))

    def as_json(self):
        return json.dumps({"source": self.source, "paths_file": self.paths.source, "config": self.data}, indent=1)


def check_required(raw, where, required=REQUIRED):
    """Raise ConfigError for a missing key of `required` ({table: keys}); the chain needs all of REQUIRED."""
    for section, keys in required.items():
        missing = [k for k in keys if k not in raw.get(section, {})]
        if missing:
            raise ConfigError("{}: [{}] needs {}".format(where, section, ", ".join(missing)))
    rot = raw.get("rotation", {})
    if "rotation" in required and rot.get("rule") == "corner" and "corner_offset_mm" not in rot:
        raise ConfigError("{}: [rotation] rule = \"corner\" needs corner_offset_mm".format(where))


def load_config(path, paths=None, required=None):
    """Read a config (see the module docstring). paths: a Paths, or the --paths file (None: the search).
    required: {table: keys} that must be present (the chain passes REQUIRED); None checks only the table names."""
    raw = read_toml(path)
    unknown = [k for k in raw if k not in SECTIONS]
    if unknown:
        raise ConfigError("{}: unknown table(s) {} (known: {})".format(path, ", ".join(unknown), ", ".join(SECTIONS)))
    for section, defaults in DEFAULTS.items():
        raw[section] = _merge(defaults, raw.get(section, {}))
    if required:
        check_required(raw, path, required)
    for step in raw.get("steps", {}):
        if step not in STEP_NAMES:
            raise ConfigError("{}: [steps] {} is not a step (known: {})".format(path, step, ", ".join(STEP_NAMES)))
    if not isinstance(paths, Paths):
        paths = load_paths(paths)
    data = {}
    for section, table in raw.items():
        table = dict(table)
        for key in PATH_KEYS.get(section, ()):
            if key in table:
                table[key] = paths.resolve(table[key])
        if section == "steps":
            table = paths.resolve_tree(table)
        data[section] = table
    return Config(data, os.path.abspath(path), paths)
