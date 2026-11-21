"""Machine settings: no built-in machine, the environment round trip, exact RF in Hz."""
import pytest

from spyral_inflector import machine


def test_missing_setting_is_an_error(clean_env):
    for reader in (machine.species, machine.rf_mhz, machine.rf_hz, machine.z_start):
        with pytest.raises(machine.MachineSettingMissing):
            reader()
    assert machine.inject_core() is False
    assert machine.tag() == "" and machine.name() == ""


def test_round_trip(clean_env, monkeypatch):
    env = machine.to_env({"species": "H2_1+", "rf_mhz": 32.8, "z_start_m": -0.27, "inject_core": True,
                          "tag": "M1", "name": "Machine 1"})
    assert env == {"SI_SPECIES": "H2_1+", "SI_RF_MHZ": "32.8", "SI_Z_START": "-0.27", "SI_INJECT_CORE": "1",
                   "SI_MACHINE": "M1", "SI_MACHINE_NAME": "Machine 1"}
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    assert machine.species() == "H2_1+"
    assert machine.rf_mhz() == 32.8
    assert machine.z_start() == -0.27
    assert machine.inject_core() is True
    assert machine.tag() == "M1" and machine.name() == "Machine 1"


def test_rf_hz_is_exact(clean_env, monkeypatch):
    monkeypatch.setenv("SI_RF_MHZ", "32.8")
    assert machine.rf_hz() == 32.8e6                 # not 32.8 * 1e6 = 32799999.999999996
    assert 32.8 * 1e6 != 32.8e6


def test_name_falls_back_to_tag(clean_env, monkeypatch):
    monkeypatch.setenv("SI_MACHINE", "M2")
    assert machine.name() == "M2"


def test_unknown_key(clean_env):
    with pytest.raises(KeyError, match="unknown machine setting"):
        machine.to_env({"harmonic": 4})
