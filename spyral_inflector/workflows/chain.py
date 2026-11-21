"""The levelling chain of one design, driven by a workflow config (workflows.config).

Pass k (run <prefix><k> under <results>/<family>):

    rotate   the rotation fixed point of the plate rule and the basis solves, from the previous pass's levelled
             setting (pass 1: the source run, optionally with an optimizer winner's knobs)
    level    exit truncation + dz levelling on the space-charge centroid at the operating current
    export   the Baseline-frame export; the plate rule re-solves R on the levelled STEP set

The chain repeats while the export's R differs from the run's by dr_tol_deg or more (at most max_passes). The
converged pass then runs the final steps: level --final (every particle), export, the electrode-gap check, the
beam report and the phase coupling against the reference run, and the renders.

Every step is a subprocess ``python -u -m spyral_inflector.workflows.<step> ...``. The machine settings of the config
go into the environment of every step; the CPU-only steps (export, gap, reports, renders) run with
CUDA_VISIBLE_DEVICES="". Step output is appended to <logs>/level.log; chain.log records each command line.

    python -m spyral_inflector.workflows.chain --config design.toml [--paths paths.toml] [--dry-run] [--results <dir>]
"""
import argparse
import datetime
import json
import os
import subprocess
import sys
from decimal import Decimal, ROUND_HALF_UP

from .. import machine
from .config import REQUIRED, load_config

GPU_STEPS = ("rotate", "level")


def fmt(value):
    """A config value as a command-line string: floats as repr() (reads back as the same float), ints as they are."""
    if isinstance(value, bool):
        raise TypeError("a boolean is not a command-line number: {!r}".format(value))
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _jload(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


class Chain:
    """Builds and runs the steps of one config. dry_run prints the command lines instead (pass 1 counts as
    converged, as the legacy runners did)."""

    def __init__(self, cfg, dry_run=False, python=None, out=None, results=None):
        self.cfg = cfg
        self.dry_run = dry_run
        self.python = python or sys.executable
        self.out = out or (lambda m: print(m, flush=True))
        run = cfg["run"]
        self.results = os.path.abspath(results) if results else run["results"]
        self.family_dir = os.path.join(self.results, run["family"])
        self.logs = (run.get("logs") if not results else None) or os.path.join(self.results, "chain_logs", run["prefix"])
        self.commands = []          # (label, argv, cpu_only, log name), in order

    # ----------------------------------------------------------------------------------------------- settings
    def machine_env(self):
        return machine.to_env(self.cfg["machine"])

    def env(self, cpu_only):
        env = dict(os.environ)
        env.update(self.machine_env())
        env["SI_RESULTS"] = self.results
        env["PYTHONUNBUFFERED"] = "1"
        env["MPLBACKEND"] = "Agg"
        env.pop("MPI4PY_RC_INITIALIZE", None)        # set, it stops mpi4py's import in the BEM steps
        if cpu_only:
            env["CUDA_VISIBLE_DEVICES"] = ""
        else:
            env.pop("CUDA_VISIBLE_DEVICES", None)
        return env

    def source_levelling(self):
        """Basis turns and total turns of the quads as strings (from the config, else the source's levelling)."""
        run = self.cfg["run"]
        lv_path = os.path.join(run["source"], "sc_level", "sc_level.json")
        lv = _jload(lv_path) if os.path.isfile(lv_path) else None
        if "quad_rotation_deg" in run:
            qrot = [float(v) for v in run["quad_rotation_deg"]]
        elif lv is not None:
            qrot = [float(v) for v in lv["quad_rotation_deg"]]
        else:
            raise FileNotFoundError("{} missing: give [run] quad_rotation_deg".format(lv_path))
        if "quad_total_deg" in run:
            qtot = [float(v) for v in run["quad_total_deg"]]
        elif lv is not None:
            qtot = [qrot[i] + float(lv["quads"][2 * i + 1]) for i in range(len(qrot))]
        else:
            raise FileNotFoundError("{} missing: give [run] quad_total_deg".format(lv_path))
        return [fmt(v) for v in qrot], [fmt(v) for v in qtot]

    def vscale(self):
        run = self.cfg["run"]
        if "vscale" in run:
            return fmt(float(run["vscale"]))
        if run.get("best_knobs"):
            return fmt(float(_jload(run["best_knobs"])["_vscale"]))
        lv = _jload(os.path.join(run["source"], "sc_level", "sc_level.json"))
        return fmt(float(lv["vscale"]))

    def rule_args(self, with_rule=True):
        """The corner-rule arguments (the rotate step implements the corner rule only; export takes --rotation-rule)."""
        rot = self.cfg["rotation"]
        if rot["rule"] != "corner":
            raise ValueError("the chain's rotate step implements the corner rule only, not [rotation] rule = {!r}".format(rot["rule"]))
        out = ["--rotation-rule", "corner"] if with_rule else []
        return out + ["--gap-azimuth", fmt(rot["gap_azimuth_deg"]), "--corner", rot["corner"],
                      "--corner-offset-mm", fmt(rot["corner_offset_mm"]), "--half-gap", fmt(rot["half_gap_m"])]

    def current(self):
        return fmt(self.cfg["inputs"]["current_ma"])

    # --------------------------------------------------------------------------------------------- the steps
    def rotate_cmd(self, prev, name, first):
        run, inp = self.cfg["run"], self.cfg["inputs"]
        cmd = ["--prev", prev, "--name", name, "--results", self.results, "--family", run["family"],
               "--bfield", inp["bfield"], "--particles", inp["particles"]] + self.rule_args(with_rule=False)
        if first:
            if run.get("best_knobs"):
                cmd += ["--best-knobs", run["best_knobs"]]
            if run.get("r0_deg") is not None:
                cmd += ["--r0={}".format(fmt(float(run["r0_deg"])))]
        return cmd + self.cfg.step_args("rotate")

    def level_cmd(self, lrun, qrot, vscale, final=False):
        inp = self.cfg["inputs"]
        cmd = (["--run", lrun, "--best-json", self.cfg["run"]["quads_json"], "--quads"] + qrot +
               ["--base-run", lrun, "--bfield", inp["bfield"], "--particles", inp["particles"],
                "--current-ma", self.current(), "--vscale", vscale] + self.cfg.step_args("level"))
        return cmd + (["--final"] + self.cfg.step_args("level_final") if final else [])

    def export_cmd(self, lrun, qtot, final=False):
        return (["--run", lrun, "--quads"] + qtot + self.rule_args() + self.cfg.step_args("export") +
                (self.cfg.step_args("export_final") if final else []))

    def levelled_steps_dir(self, lrun):
        """The STEP folder the final level run tracked in: sc_level/t<t_exit:.3f>/steps, else <run>/steps."""
        default = os.path.join(lrun, "steps")
        lv_path = os.path.join(lrun, "sc_level", "sc_level.json")
        if self.dry_run or not os.path.isfile(lv_path):
            return default
        t = float(_jload(lv_path)["t_exit_deg"])
        names = ["t{:.3f}".format(t), "t{}".format(Decimal(repr(t)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))]
        for name in names:
            cand = os.path.join(lrun, "sc_level", name, "steps")
            if os.path.isdir(cand):
                return cand
        return default

    def final_cmds(self, lname, lrun, qtot):
        rep = self.cfg.data.get("reports", {})
        ref_run = rep.get("reference_run") or self.cfg["run"]["source"]
        ref_label = rep.get("reference_label") or "{} sc".format(os.path.basename(os.path.normpath(ref_run)))
        run_label = "{} {}mA".format(lname, self.current())
        ref_npz = os.path.join(ref_run, "sc_level", "bunch_sc_all_level.npz")
        run_npz = os.path.join(lrun, "sc_level", "bunch_sc_all_level.npz")
        mass = ["--mass-kg", fmt(float(rep["mass_kg"]))] if "mass_kg" in rep else []
        tdir = self.levelled_steps_dir(lrun)
        out = []
        gap = [tdir, "--labels", "{} levelled".format(lname)]
        if "gap_voltage_kv" in rep:
            gap += ["--voltage-kv", fmt(rep["gap_voltage_kv"])]
        out.append(("electrode gap " + lname, "electrode_gap", gap + self.cfg.step_args("gap"), "gap.log"))
        out.append(("beam report " + lname, "beam_report",
                    [ref_npz, run_npz, "--labels", ref_label, run_label, "--current-ma", self.current()] + mass +
                    self.cfg.step_args("report"), "report.log"))
        out.append(("phase coupling " + lname, "phase_coupling",
                    [ref_npz, run_npz, "--labels", ref_label, run_label] + mass + self.cfg.step_args("phase"), "report.log"))
        if rep.get("render", True):
            out.append(("render " + lname, "render",
                        [run_npz, "--steps", tdir, "--quad-rotations"] + qtot +
                        ["--out", os.path.join(lrun, "renders", lname), "--current-ma", self.current(),
                         "--title", "{} levelled, every particle {} mA".format(lname, self.current())] +
                        self.cfg.step_args("render"), "level.log"))
        return out

    # ---------------------------------------------------------------------------------------------- running
    def say(self, msg):
        line = "[{}] {}".format(datetime.datetime.now().strftime("%m-%d %H:%M:%S"), msg)
        if not self.dry_run:
            with open(os.path.join(self.logs, "chain.log"), "a", encoding="utf-8") as fh:
                fh.write(line + "\n")
        self.out(msg)

    def step(self, label, module, args, log="level.log", cpu_only=False):
        argv = [self.python, "-u", "-m", "spyral_inflector.workflows." + module] + list(args)
        self.commands.append((label, argv, cpu_only, log))
        if self.dry_run:
            self.out("DRY {}: {}{}".format(label, "CUDA_VISIBLE_DEVICES= " if cpu_only else "", " ".join(argv[1:])))
            return 0
        self.say("RUN {}: python {}".format(label, " ".join(argv[1:])))
        with open(os.path.join(self.logs, log), "a", encoding="utf-8") as fh:
            rc = subprocess.run(argv, stdout=fh, stderr=subprocess.STDOUT, env=self.env(cpu_only), cwd=self.logs).returncode
        self.say("{} exit {}".format(label, rc))
        if rc != 0:
            self.say("STOP: {} failed - see {}".format(label, os.path.join(self.logs, log)))
            raise SystemExit(1)
        return rc

    def run(self):
        cfg, run = self.cfg, self.cfg["run"]
        if not self.dry_run:
            os.makedirs(self.family_dir, exist_ok=True)
            os.makedirs(self.logs, exist_ok=True)
            self.say("config {} (paths file {}); python {}".format(cfg.source, cfg.paths.source, sys.version.split()[0]))
            self.say("machine " + " ".join("{}={}".format(k, v) for k, v in sorted(self.machine_env().items())))
        qrot, qtot = self.source_levelling()
        vscale = self.vscale()
        self.say("LEVEL {} from {} {}: quads {}, vscale {}, basis turns {}, total turns {}, B map {}, particles {}, {} mA".format(
            run["prefix"], run["source"], "(optimizer winner {})".format(run["best_knobs"]) if run.get("best_knobs") else
            "(its levelled setting)", run["quads_json"], vscale, " ".join(qrot), " ".join(qtot), cfg["inputs"]["bfield"],
            cfg["inputs"]["particles"], self.current()))
        prev = run["source"]
        for k in range(1, int(run["max_passes"]) + 1):
            lname = "{}{}".format(run["prefix"], k)
            lrun = os.path.join(self.family_dir, lname)
            self.step("rotate " + lname, "rotate", self.rotate_cmd(prev, lname, k == 1))
            self.step("level " + lname, "level", self.level_cmd(lrun, qrot, vscale))
            self.step("export " + lname, "export", self.export_cmd(lrun, qtot), cpu_only=True)
            d_r = 0.0
            if not self.dry_run:
                ex = _jload(os.path.join(lrun, "export_baseline_levelled", "levelled_export.json"))
                d_r = float(ex["rotation_deg"]) - float(ex["rotation_run_deg"])
                self.say("pass {} ({}): run R {:.4f}, export R {:.4f}, dR {:+.4f}".format(
                    k, lname, float(ex["rotation_run_deg"]), float(ex["rotation_deg"]), d_r))
            if abs(d_r) < float(run["dr_tol_deg"]):
                self.say("pass {} converged (|dR| < {:g} deg): final runs".format(k, float(run["dr_tol_deg"])))
                self.step("level final " + lname, "level", self.level_cmd(lrun, qrot, vscale, final=True))
                self.step("export final " + lname, "export", self.export_cmd(lrun, qtot, final=True), cpu_only=True)
                for label, module, args, log in self.final_cmds(lname, lrun, qtot):
                    self.step(label, module, args, log=log, cpu_only=True)
                self.say("levelled: " + lrun)
                return lrun
            prev = lrun
        self.say("STOP: {} did not converge within pass {}".format(run["prefix"], run["max_passes"]))
        raise SystemExit(1)


def main(argv=None):
    p = argparse.ArgumentParser(description="the levelling chain of one design config")
    p.add_argument("--config", required=True, help="design config (TOML)")
    p.add_argument("--paths", default=None, help="per-machine paths file (default: the search of workflows.paths)")
    p.add_argument("--dry-run", action="store_true", help="print the command lines of pass 1 and the final steps")
    p.add_argument("--results", default=None,
                   help="write the runs (and the chain logs) under this results root instead of the config's [run] results; "
                        "the source and reference runs are still read where the config says")
    a = p.parse_args(argv)
    return Chain(load_config(a.config, a.paths, required=REQUIRED), dry_run=a.dry_run, results=a.results).run()


if __name__ == "__main__":
    main()
