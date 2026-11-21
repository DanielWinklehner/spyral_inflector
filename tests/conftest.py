"""Shared fixtures: a clean environment (no machine settings, no paths file from the user's setup)."""
import os

import pytest

from spyral_inflector import machine
from spyral_inflector.workflows import paths as paths_mod


@pytest.fixture
def clean_env(monkeypatch, tmp_path):
    """No SI_* machine settings, no paths-file variables, and a home folder without ~/.config/spyral_inflector."""
    for var in list(machine.ENV.values()) + [paths_mod.ENV_SHARED, paths_mod.ENV_PACKAGE, "SI_RESULTS"]:
        monkeypatch.setenv(var, "placeholder")      # records the original state, so a value the code under test
        monkeypatch.delenv(var)                     # sets is removed again after the test
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    return tmp_path


def write(path, text):
    path = os.fspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path
