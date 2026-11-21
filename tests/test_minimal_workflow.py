"""The neutral worked example (Examples/minimal_workflow) end to end: config + paths file + machine settings,
generated field and beam, the beam read through the machine settings, the analytic design."""
import importlib.util
import json
import os

import numpy as np

from conftest import write

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(os.path.dirname(HERE), "Examples", "minimal_workflow")


def test_minimal_workflow(clean_env):
    work = os.path.join(os.fspath(clean_env), "work")
    paths = write(os.path.join(os.fspath(clean_env), "paths.toml"), "[roots]\nwork = '{}'\n".format(work))
    spec = importlib.util.spec_from_file_location("run_minimal", os.path.join(EXAMPLE, "run_minimal.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    summary = mod.main(["--config", os.path.join(EXAMPLE, "design.toml"), "--paths", paths])
    out = os.path.join(work, "minimal_workflow")
    assert os.path.isfile(os.path.join(out, "summary.json"))
    with open(os.path.join(out, "summary.json"), encoding="utf-8") as fh:
        assert json.load(fh)["n_particles"] == summary["n_particles"] == 400
    # the bunch was placed at the machine's start plane (z_start_m = -0.10 in machine.toml)
    assert np.isclose(np.mean(summary["bunch_z_range_m"]), -0.10, atol=0.02)
    assert summary["bz_centre_T"] < 0 and summary["height_m"] > 0
