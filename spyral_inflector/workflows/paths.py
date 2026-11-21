"""The per-machine paths file: where the disk roots of the work are on this computer.

One TOML file per computer, shared by every package of the chain. A config names a root, never a drive letter::

    # paths.toml
    [meta]                      # optional
    name = "desktop"

    [roots]
    work   = 'F:\\Work'          # TOML literal strings keep Windows backslashes as they are
    inputs = '/data/inputs'

    [spyral_inflector]          # optional: settings only this package reads (other packages use their own table)

A string in a config may contain ``${name}``; it is replaced by that root. ``$$`` is a literal ``$``. A path with no
``${...}`` (absolute or relative) is used unchanged. An unknown name is an error that names the paths file read.

The file is found from, first match wins: the ``--paths`` argument (``explicit``), the shared environment variable
``MACHINE_PATHS``, the package variable ``SPYRAL_INFLECTOR_PATHS``, ``~/.config/spyral_inflector/paths.toml``.
No file is fine until a ``${name}`` has to be resolved. Roots are not checked for existence when the file is read
(a root may be an offline archive disk).
"""
import os
import re
import tomllib

PACKAGE = "spyral_inflector"
ENV_SHARED = "MACHINE_PATHS"
ENV_PACKAGE = "SPYRAL_INFLECTOR_PATHS"
_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
_REF = re.compile(r"\$\$|\$\{([^}]*)\}")


class PathsError(ValueError):
    pass


def default_locations():
    """The places searched after the explicit argument, in order, as (label, path or None)."""
    return [(ENV_SHARED, os.environ.get(ENV_SHARED) or None),
            (ENV_PACKAGE, os.environ.get(ENV_PACKAGE) or None),
            ("~/.config/{}/paths.toml".format(PACKAGE), os.path.join(os.path.expanduser("~"), ".config", PACKAGE, "paths.toml"))]


def find_paths_file(explicit=None):
    """The paths file to read, or None. An explicit file or one named by an environment variable must exist."""
    if explicit:
        if not os.path.isfile(explicit):
            raise PathsError("paths file {} (--paths) does not exist".format(explicit))
        return os.path.abspath(explicit)
    for label, path in default_locations():
        if path is None:
            continue
        if os.path.isfile(path):
            return os.path.abspath(path)
        if label != default_locations()[-1][0]:
            raise PathsError("paths file {} (from {}) does not exist".format(path, label))
    return None


class Paths:
    """The roots of one paths file (or none) and the substitution of ${name} in config strings."""

    def __init__(self, roots=None, source=None, tables=None, searched=None):
        self.roots = dict(roots or {})
        self.source = source
        self.tables = dict(tables or {})
        self.searched = list(searched or [])

    def package_settings(self, package=PACKAGE):
        """This package's own table of the paths file ({} when absent)."""
        return dict(self.tables.get(package, {}))

    def _where(self):
        if self.source:
            return "paths file {} (roots: {})".format(self.source, ", ".join(sorted(self.roots)) or "none")
        return "no paths file found (searched: --paths, {})".format(", ".join(
            "{} ({})".format(label, path) if path else label for label, path in self.searched))

    def resolve(self, value):
        """Replace ${name} in a string; anything else is returned unchanged. A resolved path is normalised."""
        if not isinstance(value, str) or "$" not in value:
            return value
        used = []

        def sub(m):
            if m.group(0) == "$$":
                return "$"
            name = m.group(1)
            if name not in self.roots:
                raise PathsError("unknown root ${{{}}} in '{}': {}".format(name, value, self._where()))
            used.append(name)
            return self.roots[name]

        out = _REF.sub(sub, value)
        return os.path.normpath(out) if used else out

    def resolve_tree(self, obj):
        """resolve() applied to every string of a nested dict / list."""
        if isinstance(obj, dict):
            return {k: self.resolve_tree(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [self.resolve_tree(v) for v in obj]
        return self.resolve(obj)


def load_paths(explicit=None):
    """Read the paths file (see the module docstring). Returns a Paths, empty when no file was found."""
    source = find_paths_file(explicit)
    searched = default_locations()
    if source is None:
        return Paths(searched=searched)
    with open(source, "rb") as fh:
        data = tomllib.load(fh)
    roots = data.get("roots", {})
    if not isinstance(roots, dict):
        raise PathsError("{}: [roots] must be a table".format(source))
    out = {}
    for name, path in roots.items():
        if not _NAME.match(name):
            raise PathsError("{}: root name '{}' is not an identifier".format(source, name))
        if not isinstance(path, str):
            raise PathsError("{}: root '{}' must be a string".format(source, name))
        path = os.path.expanduser(path)
        if not (os.path.isabs(path) or re.match(r"[A-Za-z]:[\\/]", path)):
            raise PathsError("{}: root '{}' = '{}' is not an absolute path".format(source, name, path))
        out[name] = path
    for key, value in data.items():
        if key not in ("meta", "roots") and not isinstance(value, dict):
            raise PathsError("{}: top-level key '{}': package settings go in a table named after the package, "
                             "e.g. [{}]".format(source, key, PACKAGE))
    tables = {k: v for k, v in data.items() if k not in ("meta", "roots")}
    return Paths(out, source=source, tables=tables, searched=searched)
