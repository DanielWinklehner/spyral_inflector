"""The per-machine paths file: lookup order, ${root} substitution and the errors."""
import os

import pytest

from spyral_inflector.workflows.paths import PathsError, load_paths, find_paths_file

from conftest import write


def _paths_file(folder, name="paths.toml", work=None):
    work = work or os.path.join(os.fspath(folder), "work")
    return write(os.path.join(os.fspath(folder), name),
                 "[meta]\nname = 'test'\n\n[roots]\nwork = '{}'\ninputs = '{}'\n\n[spyral_inflector]\nnote = 'x'\n".format(
                     work, os.path.join(os.fspath(folder), "inputs")))


def test_no_file_is_fine_until_a_root_is_used(clean_env):
    p = load_paths()
    assert p.source is None
    assert p.resolve("relative/path") == "relative/path"
    with pytest.raises(PathsError, match="no paths file found"):
        p.resolve("${work}/x")


def test_resolve_and_normalise(clean_env):
    f = _paths_file(clean_env)
    p = load_paths(f)
    assert p.source == os.path.abspath(f)
    assert p.resolve("${work}/Results/run") == os.path.normpath(os.path.join(os.fspath(clean_env), "work", "Results", "run"))
    assert p.resolve("cost $$5") == "cost $5"
    absolute = os.path.join(os.fspath(clean_env), "elsewhere")
    assert p.resolve(absolute) == absolute                 # absolute paths keep working
    assert p.package_settings() == {"note": "x"}


def test_unknown_root_names_the_file(clean_env):
    f = _paths_file(clean_env)
    with pytest.raises(PathsError) as err:
        load_paths(f).resolve("${archive}/x")
    assert "archive" in str(err.value) and os.path.abspath(f) in str(err.value)


def test_lookup_order(clean_env, monkeypatch):
    shared = _paths_file(clean_env / "shared")
    package = _paths_file(clean_env / "package")
    explicit = _paths_file(clean_env / "explicit")
    home = _paths_file(clean_env / "home" / ".config" / "spyral_inflector")
    assert find_paths_file() == os.path.abspath(home)
    monkeypatch.setenv("SPYRAL_INFLECTOR_PATHS", package)
    assert find_paths_file() == os.path.abspath(package)
    monkeypatch.setenv("MACHINE_PATHS", shared)
    assert find_paths_file() == os.path.abspath(shared)
    assert find_paths_file(explicit) == os.path.abspath(explicit)


def test_named_file_must_exist(clean_env, monkeypatch):
    with pytest.raises(PathsError, match="--paths"):
        find_paths_file(os.path.join(os.fspath(clean_env), "missing.toml"))
    monkeypatch.setenv("MACHINE_PATHS", os.path.join(os.fspath(clean_env), "missing.toml"))
    with pytest.raises(PathsError, match="MACHINE_PATHS"):
        find_paths_file()


@pytest.mark.parametrize("body, message", [
    ("[roots]\nwork = 'relative/dir'\n", "not an absolute path"),
    ("[roots]\n'1bad' = '/abs'\n", "not an identifier"),
    ("[roots]\nwork = 3\n", "must be a string"),
    ("results = '/abs'\n[roots]\nwork = '/abs'\n", "package settings go in a table"),
])
def test_invalid_files(clean_env, body, message):
    f = write(os.path.join(os.fspath(clean_env), "bad.toml"), body)
    with pytest.raises(PathsError, match=message):
        load_paths(f)
