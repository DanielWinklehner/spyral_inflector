"""Design configs and the chain: includes, ${root} resolution, and the command lines of a dry run.

The expected command lines are the ones the PowerShell levelling runner (Invoke-Level) built for the same inputs,
with each deck script replaced by its workflows module; the rotate step also gets the corner-rule arguments the
deck script took as defaults, and --results/--family explicitly."""
import json
import os
import sys

import pytest

from spyral_inflector.workflows.chain import Chain, fmt
from spyral_inflector.workflows.config import REQUIRED, ConfigError, load_config

from conftest import write


def _setup(root):
    root = os.fspath(root)
    work = os.path.join(root, "work")
    src = os.path.join(work, "Results", "triplet", "d1")
    write(os.path.join(src, "sc_level", "sc_level.json"),
          json.dumps({"quad_rotation_deg": [18.51, 22.71, 28.46], "quads": [9302.0, -4.55921, -9917.0, -3.85248, -6858.0, -6.31506],
                      "vscale": 1.25, "t_exit_deg": 0.851555}))
    write(os.path.join(src, "prof", "best_knobs.json"), json.dumps({"_vscale": 1.0506682648059003, "tilt": 39.0}))
    write(os.path.join(src, "prof", "best_quads.json"), "{}")
    paths = write(os.path.join(root, "paths.toml"), "[roots]\nwork = '{}'\ninputs = '{}'\n".format(work, os.path.join(root, "inputs")))
    write(os.path.join(root, "cfg", "machine.toml"),
          "[machine]\nspecies = 'H2_1+'\nrf_mhz = 32.8\nz_start_m = -0.270\ninject_core = true\ntag = 'M1'\n")
    design = write(os.path.join(root, "cfg", "d1o.toml"), """include = ["machine.toml"]
[run]
results = "${work}/Results"
prefix = "d1o_c"
source = "${work}/Results/triplet/d1"
best_knobs = "${work}/Results/triplet/d1/prof/best_knobs.json"
quads_json = "${work}/Results/triplet/d1/prof/best_quads.json"
vscale = 1.0506682648059005
[inputs]
bfield = "${work}/Fields/map.pickle"
particles = "${inputs}/Particles/beam.txt"
current_ma = 9.434
[rotation]
rule = "corner"
gap_azimuth_deg = 26
corner_offset_mm = 3.25
half_gap_m = 0.00325
[reports]
gap_voltage_kv = 24.8
mass_kg = 3.3462e-27
render = false
""")
    return work, src, paths, design


def test_include_and_resolution(clean_env):
    work, src, paths, design = _setup(clean_env)
    cfg = load_config(design, paths, required=REQUIRED)
    assert cfg["machine"]["z_start_m"] == -0.27 and cfg["machine"]["tag"] == "M1"
    assert cfg["run"]["source"] == os.path.normpath(src)
    assert cfg["run"]["family"] == "triplet" and cfg["run"]["max_passes"] == 4      # defaults
    assert cfg["inputs"]["particles"] == os.path.normpath(os.path.join(os.fspath(clean_env), "inputs", "Particles", "beam.txt"))


def test_config_errors(clean_env):
    work, src, paths, design = _setup(clean_env)
    bad = write(os.path.join(os.fspath(clean_env), "cfg", "bad.toml"), "include = ['machine.toml']\n[nonsense]\nx = 1\n")
    with pytest.raises(ConfigError, match="unknown table"):
        load_config(bad, paths)
    partial = write(os.path.join(os.fspath(clean_env), "cfg", "partial.toml"), "include = ['machine.toml']\n[run]\nprefix = 'x'\n")
    load_config(partial, paths)                                   # fine without the chain's requirements
    with pytest.raises(ConfigError, match=r"\[run\] needs"):
        load_config(partial, paths, required=REQUIRED)
    loop = write(os.path.join(os.fspath(clean_env), "cfg", "loop.toml"), "include = ['loop.toml']\n")
    with pytest.raises(ConfigError, match="include loop"):
        load_config(loop, paths)


def test_fmt_round_trips():
    for v in (1.0506682648059005, 0.00325, 9.434, 18.51 + -4.55921, -0.27):
        assert float(fmt(v)) == v
    assert fmt(26) == "26"
    with pytest.raises(TypeError):
        fmt(True)


def test_dry_run_matches_the_legacy_runner(clean_env):
    work, src, paths, design = _setup(clean_env)
    cfg = load_config(design, paths, required=REQUIRED)
    chain = Chain(cfg, dry_run=True, out=lambda m: None)
    lrun = chain.run()
    res = os.path.normpath(os.path.join(work, "Results"))
    run = os.path.join(res, "triplet", "d1o_c1")
    bf = os.path.normpath(os.path.join(work, "Fields", "map.pickle"))
    parts = os.path.normpath(os.path.join(os.fspath(clean_env), "inputs", "Particles", "beam.txt"))
    bkf = os.path.join(src, "prof", "best_knobs.json")
    bqf = os.path.join(src, "prof", "best_quads.json")
    qrot = ["18.51", "22.71", "28.46"]
    qtot = [repr(18.51 + -4.55921), repr(22.71 + -3.85248), repr(28.46 + -6.31506)]
    rule = ["--rotation-rule", "corner", "--gap-azimuth", "26", "--corner", "outer", "--corner-offset-mm", "3.25", "--half-gap", "0.00325"]
    lvl = ["--run", run, "--best-json", bqf, "--quads"] + qrot + ["--base-run", run, "--bfield", bf, "--particles", parts,
                                                                 "--current-ma", "9.434", "--vscale", "1.0506682648059005"]
    exp = ["--run", run, "--quads"] + qtot + rule
    ref_npz = os.path.join(src, "sc_level", "bunch_sc_all_level.npz")
    run_npz = os.path.join(run, "sc_level", "bunch_sc_all_level.npz")
    expected = [
        ("rotate", ["--prev", src, "--name", "d1o_c1", "--results", res, "--family", "triplet", "--bfield", bf, "--particles", parts,
                    "--gap-azimuth", "26", "--corner", "outer", "--corner-offset-mm", "3.25", "--half-gap", "0.00325",
                    "--best-knobs", bkf], False),
        ("level", lvl, False),
        ("export", exp, True),
        ("level", lvl + ["--final"], False),
        ("export", exp, True),
        ("electrode_gap", [os.path.join(run, "steps"), "--labels", "d1o_c1 levelled", "--voltage-kv", "24.8"], True),
        ("beam_report", [ref_npz, run_npz, "--labels", "d1 sc", "d1o_c1 9.434mA", "--current-ma", "9.434", "--mass-kg", "3.3462e-27"], True),
        ("phase_coupling", [ref_npz, run_npz, "--labels", "d1 sc", "d1o_c1 9.434mA", "--mass-kg", "3.3462e-27"], True),
    ]
    got = [(argv[3].rsplit(".", 1)[1], argv[4:], cpu) for _, argv, cpu, _ in chain.commands]
    assert lrun == run
    assert [c[1] for c in chain.commands][0][:3] == [sys.executable, "-u", "-m"]
    assert got == expected


def test_machine_env_of_a_step(clean_env):
    work, src, paths, design = _setup(clean_env)
    chain = Chain(load_config(design, paths, required=REQUIRED), dry_run=True, out=lambda m: None)
    env = chain.env(cpu_only=True)
    assert env["SI_SPECIES"] == "H2_1+" and env["SI_RF_MHZ"] == "32.8" and env["SI_Z_START"] == "-0.27"
    assert env["SI_INJECT_CORE"] == "1" and env["SI_MACHINE"] == "M1" and env["CUDA_VISIBLE_DEVICES"] == ""
    assert "CUDA_VISIBLE_DEVICES" not in chain.env(cpu_only=False)


def test_vscale_default_is_the_json_value(clean_env):
    work, src, paths, design = _setup(clean_env)
    cfg = load_config(design, paths, required=REQUIRED)
    del cfg.data["run"]["vscale"]
    assert Chain(cfg, dry_run=True).vscale() == "1.0506682648059003"
