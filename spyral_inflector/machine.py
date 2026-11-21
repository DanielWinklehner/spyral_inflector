"""Machine settings: the few numbers that belong to the cyclotron and its injector, not to the inflector design.

The package has no built-in machine. A run gets these from the ``[machine]`` table of its workflow config
(``spyral_inflector.workflows``), which exports them as environment variables, so every subprocess of a run sees the
same values:

=================  ==================  =====================================================================
setting            variable            meaning
=================  ==================  =====================================================================
``species``        ``SI_SPECIES``      PyPATools species name of the beam, e.g. ``H2_1+``
``rf_mhz``         ``SI_RF_MHZ``       cyclotron RF [MHz]: one bunch carries I / f_RF of charge, and a
                                       particle file's phase is converted to arrival time with it
``z_start_m``      ``SI_Z_START``      where the injected bunch starts, deck frame [m] (negative: upstream)
``inject_core``    ``SI_INJECT_CORE``  "1": every particle is injected at the start plane at its own arrival
                                       time (tracking.bunch --inject-core is then the default)
``tag``            ``SI_MACHINE``      short name used in exported file names (optional)
``name``           ``SI_MACHINE_NAME`` display name written into headers of exported files (optional)
=================  ==================  =====================================================================

A setting that a computation needs and nobody gave is an error naming the variable, never a silent default.
"""
import os
from decimal import Decimal

ENV = {
    "species": "SI_SPECIES",
    "rf_mhz": "SI_RF_MHZ",
    "z_start_m": "SI_Z_START",
    "inject_core": "SI_INJECT_CORE",
    "tag": "SI_MACHINE",
    "name": "SI_MACHINE_NAME",
}


class MachineSettingMissing(RuntimeError):
    pass


def _require(key):
    var = ENV[key]
    value = os.environ.get(var)
    if value is None or value == "":
        raise MachineSettingMissing(
            "machine setting '{}' is not set: give it in the [machine] table of the workflow config, or set the "
            "environment variable {} (see spyral_inflector.machine)".format(key, var))
    return value


def species():
    """PyPATools species name of the beam (SI_SPECIES)."""
    return _require("species")


def rf_mhz():
    """Cyclotron RF frequency [MHz] (SI_RF_MHZ)."""
    return float(_require("rf_mhz"))


def rf_hz():
    """Cyclotron RF frequency [Hz], converted from the decimal MHz value without binary rounding
    (32.8 MHz -> 32800000.0 exactly; 32.8 * 1e6 in floating point is 32799999.999999996)."""
    return float(Decimal(_require("rf_mhz")) * 1000000)


def z_start():
    """Start plane of the injected bunch, deck frame [m] (SI_Z_START)."""
    return float(_require("z_start_m"))


def inject_core():
    """True when every particle is injected at the start plane at its own arrival time (SI_INJECT_CORE=1)."""
    return os.environ.get(ENV["inject_core"], "0") == "1"


def tag(default=""):
    """Short machine name for exported file names (SI_MACHINE); empty when not set."""
    return os.environ.get(ENV["tag"]) or default


def name(default=""):
    """Display name for headers of exported files (SI_MACHINE_NAME, else the tag)."""
    return os.environ.get(ENV["name"]) or tag(default)


def to_env(settings):
    """The environment variables for a [machine] table (dict). Numbers are written with repr(), which reads back
    as the same float; strings are taken as they are."""
    out = {}
    for key, value in settings.items():
        if key not in ENV:
            raise KeyError("unknown machine setting '{}' (known: {})".format(key, ", ".join(sorted(ENV))))
        if isinstance(value, bool):
            out[ENV[key]] = "1" if value else "0"
        elif isinstance(value, float):
            out[ENV[key]] = repr(value)
        else:
            out[ENV[key]] = str(value)
    return out
