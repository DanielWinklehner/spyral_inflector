"""Bit-for-bit comparison of two run folders (a regression check: same inputs, two code versions).

Compared by content, file by file (relative paths):

    .npz          every array: dtype, shape and every bit (NaN equal to NaN)
    .json         every number exactly; strings after replacing each folder's own root by <ROOT>; keys that hold
                  wall-clock times (--ignore-key) are skipped
    .h5           every dataset bit for bit and every attribute (strings rooted as for json)
    .pickle       the bytes; if they differ, the unpickled objects (arrays bit for bit, strings rooted, the memo
                  uuids of cached interpolators skipped, compiled numba functions compared by name)
    .stp/.step    the DATA section (the header carries the file name and a time stamp)
    .csv, .dat    the text, rooted
    log*.txt, .png  listed, not compared (time stamps, renderer metadata)

    python -m spyral_inflector.workflows.compare <folder A> <folder B> [--ignore-key REGEX ...] [--only GLOB ...]

Exit status 0 when every compared file is identical.
"""
import argparse
import fnmatch
import json
import os
import pickle
import re
import sys

import numpy as np

# wall-clock keys (setup_time_s, solve_s, potential_s, wall_s, ...); t_ref_s, t_rms_ns and the like are physics
TIME_KEYS = (r"(^|_)(time|times|elapsed|wall|timestamp|date|started|finished|duration)(_s|_min|_h)?$"
             r"|^(potential|solve|setup|assembly)_s$|^t_(setup|solve|run)")
# identities that are new in every process: the memo uuid of PyPATools' cached interpolators
UNSTABLE_KEYS = r"__uuid$"
SKIP = ("log*.txt", "*.png", "*.log", "*.pyc")


def _rooted(s, root):
    if not isinstance(s, str):
        return s
    for r in {root, root.replace("\\", "/"), root.replace("/", "\\")}:
        s = s.replace(r, "<ROOT>")
    return s.replace("\\", "/")


def _is_compiled(x):
    return type(x).__module__.split(".")[0] == "numba"


def _compiled_name(x):
    f = getattr(x, "py_func", None)
    return "{}.{}".format(f.__module__, f.__qualname__) if f is not None else type(x).__name__


class Diff:
    def __init__(self, root_a, root_b, ignore):
        self.root_a, self.root_b = root_a, root_b
        self.ignore = [re.compile(r) for r in ignore]
        self.out = []
        self.new_file()

    def new_file(self):
        self._seen = set()          # (id A, id B) of the objects already walked: object graphs have cycles
        self._keep = []             # ... kept alive, so that an id is never reused within one file

    def note(self, where, msg):
        self.out.append("{}: {}".format(where, msg))

    def value(self, a, b, where):
        if isinstance(a, (dict, list)) or (hasattr(a, "__dict__") and not isinstance(a, type)):
            pair = (id(a), id(b))
            if pair in self._seen:
                return
            self._seen.add(pair)
            self._keep.append((a, b))
        if isinstance(a, dict) and isinstance(b, dict):
            for k in sorted(set(a) | set(b), key=str):
                if any(r.search(str(k)) for r in self.ignore):
                    continue
                if k not in a or k not in b:
                    self.note(where, "key {!r} only in {}".format(k, "A" if k in a else "B"))
                else:
                    self.value(a[k], b[k], "{}.{}".format(where, k))
        elif isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
            if len(a) != len(b):
                self.note(where, "length {} vs {}".format(len(a), len(b)))
            else:
                for i, (x, y) in enumerate(zip(a, b)):
                    self.value(x, y, "{}[{}]".format(where, i))
        elif isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            self.array(np.asarray(a), np.asarray(b), where)
        elif isinstance(a, str) and isinstance(b, str):
            if _rooted(a, self.root_a) != _rooted(b, self.root_b):
                self.note(where, "{!r} vs {!r}".format(a[:120], b[:120]))
        elif isinstance(a, float) and isinstance(b, float):
            if not (a == b or (a != a and b != b)) or np.signbit(a) != np.signbit(b):
                self.note(where, "{!r} vs {!r}".format(a, b))
        elif _is_compiled(a) or _is_compiled(b):
            # a compiled function cached in an object (numba dispatcher): code, compared by the function it compiles
            if _compiled_name(a) != _compiled_name(b):
                self.note(where, "compiled {} vs {}".format(_compiled_name(a), _compiled_name(b)))
        elif type(a) is not type(b) and not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
            self.note(where, "type {} vs {}".format(type(a).__name__, type(b).__name__))
        elif hasattr(a, "__dict__") and hasattr(b, "__dict__") and not isinstance(a, type):
            self.value(vars(a), vars(b), where)
        elif a != b:
            self.note(where, "{!r} vs {!r}".format(a, b))

    def array(self, a, b, where):
        if a.dtype != b.dtype or a.shape != b.shape:
            self.note(where, "dtype/shape {} {} vs {} {}".format(a.dtype, a.shape, b.dtype, b.shape))
        elif a.dtype.kind in "fc":
            if a.tobytes() != b.tobytes():
                same_nan = np.isnan(a) & np.isnan(b) if a.dtype.kind == "f" else np.zeros(a.shape, bool)
                bad = (a != b) & ~same_nan
                if bad.any():
                    self.note(where, "{} of {} values differ, max |d| {:.3e}".format(int(bad.sum()), a.size,
                                                                                    float(np.max(np.abs(a[bad] - b[bad])))))
        elif a.dtype.kind in "OUS":
            for i, (x, y) in enumerate(zip(a.ravel().tolist(), b.ravel().tolist())):
                self.value(x, y, "{}[{}]".format(where, i))
        elif a.tobytes() != b.tobytes():
            self.note(where, "{} of {} values differ".format(int((a != b).sum()), a.size))


def _step_data(path):
    with open(path, "rb") as fh:
        data = fh.read()
    i = data.find(b"\nDATA;")
    return data[i:] if i >= 0 else data


def _h5_tree(path):
    import h5py
    out = {}

    def visit(name, obj):
        attrs = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in obj.attrs.items()}
        out[name + "@attrs"] = attrs
        if isinstance(obj, h5py.Dataset):
            out[name] = obj[()]
    with h5py.File(path, "r") as fh:
        out["/@attrs"] = dict(fh.attrs.items())
        fh.visititems(visit)
    return out


def compare_file(d, rel, fa, fb):
    ext = os.path.splitext(rel)[1].lower()
    if ext == ".npz":
        with np.load(fa, allow_pickle=True) as za, np.load(fb, allow_pickle=True) as zb:
            d.value({k: za[k] for k in za.files}, {k: zb[k] for k in zb.files}, rel)
    elif ext == ".json":
        with open(fa, encoding="utf-8") as ha, open(fb, encoding="utf-8") as hb:
            d.value(json.load(ha), json.load(hb), rel)
    elif ext in (".h5", ".hdf5"):
        d.value(_h5_tree(fa), _h5_tree(fb), rel)
    elif ext in (".stp", ".step"):
        if _step_data(fa) != _step_data(fb):
            d.note(rel, "STEP DATA sections differ")
    elif ext in (".pickle", ".pkl"):
        with open(fa, "rb") as ha, open(fb, "rb") as hb:
            if ha.read() == hb.read():
                return
        with open(fa, "rb") as ha, open(fb, "rb") as hb:
            d.value(pickle.load(ha), pickle.load(hb), rel)
    else:
        with open(fa, "rb") as ha, open(fb, "rb") as hb:
            ba, bb = ha.read(), hb.read()
        if ba == bb:
            return
        try:
            ta, tb = ba.decode("utf-8"), bb.decode("utf-8")
        except UnicodeDecodeError:
            d.note(rel, "bytes differ")
            return
        if _rooted(ta, d.root_a) != _rooted(tb, d.root_b):
            d.note(rel, "text differs")


def compare_dirs(a, b, ignore=(TIME_KEYS, UNSTABLE_KEYS), only=None, skip=SKIP):
    """Returns (differences, compared files, skipped files, files only in A, files only in B)."""
    a, b = os.path.abspath(a), os.path.abspath(b)
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))     # deep object graphs of the field pickles

    def files(root):
        out = set()
        for dp, _, fns in os.walk(root):
            for fn in fns:
                out.add(os.path.relpath(os.path.join(dp, fn), root).replace("\\", "/"))
        return out
    fa, fb = files(a), files(b)
    d = Diff(a, b, ignore)
    compared, skipped = [], []
    for rel in sorted(fa & fb):
        if only and not any(fnmatch.fnmatch(rel, g) for g in only):
            continue
        if any(fnmatch.fnmatch(os.path.basename(rel), g) for g in skip):
            skipped.append(rel)
            continue
        d.new_file()
        try:
            compare_file(d, rel, os.path.join(a, rel), os.path.join(b, rel))
        except Exception as err:            # an unreadable file is a difference, not a crash of the comparison
            d.note(rel, "could not compare: {}: {}".format(type(err).__name__, err))
        compared.append(rel)
    keep = (lambda r: not only or any(fnmatch.fnmatch(r, g) for g in only))
    return d.out, compared, skipped, sorted(r for r in fa - fb if keep(r)), sorted(r for r in fb - fa if keep(r))


def main(argv=None):
    p = argparse.ArgumentParser(description="bit-for-bit comparison of two run folders")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--ignore-key", nargs="*", default=[TIME_KEYS, UNSTABLE_KEYS],
                   help="regexes of json/h5/pickle keys not compared (default: wall-clock times, memo uuids)")
    p.add_argument("--only", nargs="*", default=None, help="compare only these relative-path globs")
    p.add_argument("--max-lines", type=int, default=200)
    a = p.parse_args(argv)
    diffs, compared, skipped, only_a, only_b = compare_dirs(a.a, a.b, ignore=a.ignore_key, only=a.only)
    print("compared {} files ({} skipped: logs and images); {} only in A, {} only in B; {} differences".format(
        len(compared), len(skipped), len(only_a), len(only_b), len(diffs)))
    for rel in only_a:
        print("  only in A: " + rel)
    for rel in only_b:
        print("  only in B: " + rel)
    for line in diffs[:a.max_lines]:
        print("  DIFF " + line)
    if len(diffs) > a.max_lines:
        print("  ... {} more".format(len(diffs) - a.max_lines))
    ok = not diffs and not only_a and not only_b
    print("IDENTICAL" if ok else "DIFFERENT")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
