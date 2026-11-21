"""The run-folder comparison, the command line, and the exit-plane rotation rules."""
import json
import os

import numpy as np
import pytest

from spyral_inflector.tracking.exit_plane import rotation_for_corner, rotation_for_gap, rotation_parallel_to_gap
from spyral_inflector.workflows.cli import main as cli
from spyral_inflector.workflows.compare import compare_dirs

from conftest import write


def _run(folder, value=1.0, path_root=None):
    folder = os.fspath(folder)
    os.makedirs(folder, exist_ok=True)
    np.savez(os.path.join(folder, "bunch.npz"), r=np.array([[0.0, 1.0, value]]), alive=np.array([True]))
    write(os.path.join(folder, "sc_level", "sc_level.json"),
          json.dumps({"dz_m": value, "run": os.path.join(path_root or folder, "steps"), "setup_time_s": np.random.rand()}))
    write(os.path.join(folder, "steps", "001_SI_Anode.stp"),
          "ISO-10303-21;\nHEADER;\nFILE_NAME('x','{}');\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,{}));\nENDSEC;\n".format(
              np.random.rand(), value))
    write(os.path.join(folder, "log.txt"), "time {}\n".format(np.random.rand()))


def test_identical_runs(tmp_path):
    _run(tmp_path / "a")
    _run(tmp_path / "b")
    diffs, compared, skipped, only_a, only_b = compare_dirs(tmp_path / "a", tmp_path / "b")
    assert diffs == [] and only_a == [] and only_b == []
    assert "log.txt" in skipped and len(compared) == 3


def test_one_bit_differs(tmp_path):
    _run(tmp_path / "a", 1.0)
    _run(tmp_path / "b", np.nextafter(1.0, 2.0))
    diffs = compare_dirs(tmp_path / "a", tmp_path / "b")[0]
    assert any("bunch.npz" in d for d in diffs)
    assert any("dz_m" in d for d in diffs)
    assert any("STEP DATA" in d for d in diffs)


def test_cli_paths_and_config(clean_env, capsys):
    root = os.fspath(clean_env)
    paths = write(os.path.join(root, "paths.toml"), "[roots]\nwork = '{}'\n".format(os.path.join(root, "work")))
    cfg = write(os.path.join(root, "c.toml"), "[machine]\nspecies = 'H2_1+'\n[inputs]\nbfield = '${work}/map.pickle'\n")
    assert cli(["--paths", paths, "paths"]) == 0
    assert "${work}" in capsys.readouterr().out
    assert cli(["--paths", paths, "config", "--config", cfg]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["config"]["inputs"]["bfield"] == os.path.normpath(os.path.join(root, "work", "map.pickle"))
    with pytest.raises(SystemExit):
        cli(["no-such-command"])


def test_rotation_for_corner_places_the_corner():
    corner = np.array([0.060, 0.020])
    for gap, offset in ((26.0, 3.25), (24.0, 0.0), (34.0, 5.0)):
        R = rotation_for_corner(corner, gap, offset, r_near=-90.0)
        c, s = np.cos(np.radians(R)), np.sin(np.radians(R))
        p = np.array([c * corner[0] - s * corner[1], s * corner[0] + c * corner[1]])
        n_ray = np.array([-np.sin(np.radians(gap)), np.cos(np.radians(gap))])
        assert np.isclose(p @ n_ray, -1e-3 * offset, atol=1e-12)


def test_rotation_for_gap_and_parallel():
    point, normal = np.array([0.05, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    R = rotation_for_gap(point, normal, 30.0, 0.0)
    assert np.isclose(R, 30.0)
    Rp = rotation_parallel_to_gap(np.array([1.0, 0.0, 0.0]), 30.0)
    c, s = np.cos(np.radians(Rp)), np.sin(np.radians(Rp))
    n = np.array([c, s])
    ray = np.array([np.cos(np.radians(30.0)), np.sin(np.radians(30.0))])
    assert np.isclose(abs(n @ ray), 0.0, atol=1e-12)          # the plate normal is perpendicular to the ray
