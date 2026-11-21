"""The ``spyral-inflector`` command: the workflow steps, the chain of a design config, and the helpers.

    spyral-inflector chain --config design.toml [--paths paths.toml] [--dry-run]
    spyral-inflector [--config design.toml] [--paths paths.toml] <step> [step arguments ...]
    spyral-inflector paths [--paths paths.toml]          which paths file is read, and its roots
    spyral-inflector config --config design.toml          the config with every ${root} resolved
    spyral-inflector compare <run A> <run B>              bit-for-bit comparison of two run folders

With --config, a single step runs with that config's machine settings ([machine]) in its environment, as it would
inside the chain. ``python -m spyral_inflector`` is the same command.
"""
import argparse
import importlib
import os
import sys

STEPS = {
    "rotate": ("rotate", "rotation fixed point of the corner rule and the basis solves of a new run"),
    "level": ("level", "exit truncation + dz levelling on the space-charge centroid; --final: every particle"),
    "export": ("export", "Baseline-frame export of a levelled run (STEP files, E-field, hand-off specs)"),
    "ab": ("ab", "A/B of one geometry change against a levelled run: one every-particle run"),
    "export-ab": ("export_ab", "Baseline-frame export of an A/B run"),
    "gap": ("electrode_gap", "closest approach of the two spiral electrodes"),
    "report": ("beam_report", "beam report of bunch runs (every particle, tail included)"),
    "phase": ("phase_coupling", "phase coupling of bunch runs"),
    "render": ("render", "3-D and side-view renders of a levelled run"),
    "probe": ("probe", "particle files at the start plane and in the entrance aperture of a levelled run"),
    "probe-files": ("probe_files", "particle files from a bunch npz with probe planes"),
    "handoff-direction": ("handoff_direction", "mean direction of the beam at the hand-off"),
    "compare": ("compare", "bit-for-bit comparison of two run folders"),
}


def _usage():
    width = max(len(k) for k in STEPS)
    lines = ["commands:", "  {:<{w}}  {}".format("chain", "the levelling chain of a design config", w=width),
             "  {:<{w}}  {}".format("paths", "the paths file that is read, and its roots", w=width),
             "  {:<{w}}  {}".format("config", "a design config with every root resolved", w=width)]
    lines += ["  {:<{w}}  {}".format(k, v[1], w=width) for k, v in STEPS.items()]
    return "\n".join(lines)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    p = argparse.ArgumentParser(prog="spyral-inflector", description=__doc__.split("\n\n")[0], epilog=_usage(),
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=None, help="design config (TOML); a step then runs with its [machine] settings")
    p.add_argument("--paths", default=None, help="per-machine paths file (default: MACHINE_PATHS, SPYRAL_INFLECTOR_PATHS, "
                                                 "~/.config/spyral_inflector/paths.toml)")
    p.add_argument("command", help="chain, paths, config or a step (below)")
    p.add_argument("args", nargs=argparse.REMAINDER)
    a = p.parse_args(argv)
    from .paths import load_paths
    from .config import load_config

    if a.command == "chain":
        from . import chain
        extra = (["--config", a.config] if a.config else []) + (["--paths", a.paths] if a.paths else [])
        chain.main(extra + a.args)
        return 0
    if a.command == "paths":
        q = argparse.ArgumentParser(prog="spyral-inflector paths")
        q.add_argument("--paths", default=a.paths)
        paths = load_paths(q.parse_args(a.args).paths)
        print("paths file: {}".format(paths.source or "none found"))
        for label, path in paths.searched:
            print("  searched {}: {}".format(label, path or "(not set)"))
        for name, root in sorted(paths.roots.items()):
            print("  ${{{}}} = {}{}".format(name, root, "" if os.path.isdir(root) else "   (not reachable now)"))
        return 0
    if a.command == "config":
        q = argparse.ArgumentParser(prog="spyral-inflector config")
        q.add_argument("--config", default=a.config, required=a.config is None)
        q.add_argument("--paths", default=a.paths)
        b = q.parse_args(a.args)
        print(load_config(b.config, b.paths).as_json())
        return 0
    if a.command not in STEPS:
        p.error("unknown command {!r}\n{}".format(a.command, _usage()))
    if a.config:
        from .. import machine
        os.environ.update(machine.to_env(load_config(a.config, a.paths)["machine"]))
    module = importlib.import_module("spyral_inflector.workflows." + STEPS[a.command][0])
    rc = module.main(a.args)
    return rc if isinstance(rc, int) else 0


if __name__ == "__main__":
    sys.exit(main())
