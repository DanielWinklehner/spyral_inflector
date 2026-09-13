"""Where along the line does the transverse emittance grow? (follow-up 2 of the vfocus campaign)

Post-processes bunch runs made with `bunch --snapshot-every N` (full 6-D snapshots of every
particle every N steps): at each snapshot the transverse rms emittances of the CORE particles are
computed in the frame of their mean velocity -- eps_u (horizontal transverse, u = z x t), eps_w
(w = t x u, the vertical after the bend / y-like in the axial line) and the 4-D sqrt(det) -- once
over the live particles and once over the subset that is finally transmitted (fixed set: its
emittance changes only by real dynamics, not by losses). Landmarks: the quad centres (from
reload_spiral.json's quad field check), the entrance aperture (mean z of the tail hits there
when known) and the exit (where the transmitted set starts leaving through the exit plane); the
55 mm asymptotic state is appended from the run's asym_state.

    python EmittanceAlongLine.py --npz <run>\\emittance\\bunch_vf_vac.npz <run>\\emittance\\bunch_vf_sc.npz --labels "vfocus1 vac" "vfocus1 8 mA" --out <run>\\emittance\\emittance
Writes <out>.md and <out>.png.
"""
import argparse
import json
import os

import numpy as np

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--npz", nargs="+", required=True, help="bunch npz files with snap_* arrays")
p.add_argument("--labels", nargs="+", default=None)
p.add_argument("--reload-json", nargs="+", default=None, help="reload_spiral.json per npz (quad centres); default: <npz dir>/reload_spiral.json")
p.add_argument("--out", required=True, help="output prefix (.md, .png)")
p.add_argument("--title", default="transverse rms emittance along the line (core particles)")
p.add_argument("--min-n", type=int, default=50)
a = p.parse_args()
labels = a.labels or [os.path.basename(f).replace("bunch_", "").replace(".npz", "") for f in a.npz]


def frame(v):
    t = v.mean(0)
    t /= np.linalg.norm(t)
    u = np.cross([0.0, 0.0, 1.0], t)
    if np.linalg.norm(u) < 0.3:                       # beam along the axis: u = lab x projected
        u = np.array([1.0, 0.0, 0.0]) - t[0] * t
    u /= np.linalg.norm(u)
    w = np.cross(t, u)
    return t, u, w


def emittances(r, v):
    """eps_u, eps_w, sqrt(det 4D) in mm mrad (4D: (mm mrad)^2 -> reported as its square root)."""
    t, u, w = frame(v)
    vl = v @ t
    rc = r - r.mean(0)
    X = np.column_stack([1e3 * (rc @ u), 1e3 * (v @ u) / vl, 1e3 * (rc @ w), 1e3 * (v @ w) / vl])
    S = np.cov(X.T)
    eu = np.sqrt(max(np.linalg.det(S[:2, :2]), 0.0))
    ew = np.sqrt(max(np.linalg.det(S[2:, 2:]), 0.0))
    e4 = np.sqrt(np.sqrt(max(np.linalg.det(S), 0.0)))       # (mm mrad): 4th root of det, comparable to eu, ew
    return eu, ew, e4


runs = []
for i, fn in enumerate(a.npz):
    d = np.load(fn, allow_pickle=True)
    if "snap_r" not in d.files:
        raise SystemExit("{} has no snapshots (run bunch with --snapshot-every)".format(fn))
    is_tail = d["is_tail"] if "is_tail" in d.files else np.zeros(len(d["crossed"]), dtype=bool)
    core = ~is_tail
    trans = d["crossed"] & core
    rows = []
    for k, (r, v, act) in enumerate(zip(d["snap_r"], d["snap_v"], d["snap_active"])):
        live = act & core
        sub = act & trans
        if live.sum() < a.min_n:
            continue
        eu, ew, e4 = emittances(r[live].astype(float), v[live].astype(float))
        if sub.sum() >= a.min_n and sub.sum() >= 0.9 * trans.sum():
            su, sw, s4 = emittances(r[sub].astype(float), v[sub].astype(float))
        else:
            su = sw = s4 = np.nan
        rows.append([int(d["snap_steps"][k]), 1e3 * r[live, 2].mean(), int(live.sum()), eu, ew, e4, int(sub.sum()), su, sw, s4])
    rows = np.array(rows)
    # asymptotic state 55 mm past the exit (transmitted core)
    st = d["asym_state"]
    ok = d["crossed"] & np.all(np.isfinite(st), axis=1) & core
    asym = emittances(st[ok, :3], st[ok, 3:6]) if ok.sum() >= a.min_n else (np.nan,) * 3
    rj = (a.reload_json[i] if a.reload_json else os.path.join(os.path.dirname(fn), "reload_spiral.json"))
    quads = {}
    if os.path.exists(rj):
        with open(rj) as fh:
            qc = json.load(fh).get("quad_field_check") or {}
        quads = {k: 1e3 * v["z_m"] for k, v in qc.items()}
    jf = fn.replace(".npz", ".json")
    ent = None
    if os.path.exists(jf):
        with open(jf) as fh:
            ent = (json.load(fh).get("loss_z_mean_by_electrode") or {}).get("Entrance_Aperture")
    runs.append(dict(label=labels[i], rows=rows, asym=asym, quads=quads, entrance_mm=(1e3 * ent if ent is not None else None),
                     n_trans=int(trans.sum()), n_core=int(core.sum())))


def at_z(rows, z):
    """the row whose mean z is closest to z (rows sorted by step; z decreases along the axial line)"""
    j = int(np.argmin(np.abs(rows[:, 1] - z)))
    return rows[j]


L = ["# " + a.title, "",
     "Per snapshot: transverse rms emittances of the core particles in the frame of their mean velocity, in mm mrad "
     "(eps_u horizontal transverse, eps_w the other transverse, 4D = 4th root of the 4x4 determinant, i.e. the geometric "
     "mean of the two eigen-emittances). 'live' = all particles still tracked; 'transmitted' = the fixed subset that "
     "finally leaves through the exit plane (its emittance changes by dynamics only, not by losses). Values at 55 mm "
     "past the exit from the asymptotic state.", "",
     "CAVEAT: inside the electrostatic fields (the quads sit back to back with the inflector, so everything between the "
     "RFQ exit plane and the asymptotic state is inside a potential) the kinetic velocities carry position-correlated "
     "potential energy, and the rms emittance of kinetic coordinates is NOT an invariant there: the values along the line "
     "are apparent emittances. The real growth is the field-free pair start -> 55 mm past the exit; the values in between "
     "are comparable BETWEEN runs at the same location (same fields), which is what locates the differences.", ""]
for run in runs:
    rows = run["rows"]
    L += ["## {}  ({} core particles, {} transmitted)".format(run["label"], run["n_core"], run["n_trans"]), "",
          "| where | mean z [mm] | n live | live: eps_u / eps_w / 4D | n transmitted | transmitted: eps_u / eps_w / 4D |",
          "|---|---|---|---|---|---|"]
    marks = [("start (RFQ exit plane)", rows[0, 1])]
    for name, zq in sorted(run["quads"].items(), key=lambda kv: -kv[1]):
        marks.append(("centre of {}".format(name), zq))
        marks.append(("after {} (+35 mm)".format(name), zq + 35.0))
    if run["entrance_mm"] is not None:
        marks.append(("inflector entrance aperture", run["entrance_mm"]))
    # exit: the last row with the transmitted subset complete
    ok_rows = rows[np.isfinite(rows[:, 9])]
    if len(ok_rows):
        marks.append(("inflector mid (half way entrance -> exit)", 0.5 * ((run["entrance_mm"] if run["entrance_mm"] is not None else ok_rows[0, 1]) + ok_rows[-1, 1])))
        marks.append(("inflector exit (transmitted set starts leaving)", ok_rows[-1, 1]))
    for name, z in marks:
        r = at_z(rows, z)
        L.append("| {} | {:+.0f} | {} | {:.1f} / {:.1f} / {:.1f} | {} | {} |".format(
            name, r[1], int(r[2]), r[3], r[4], r[5], int(r[6]),
            "-" if not np.isfinite(r[9]) else "{:.1f} / {:.1f} / {:.1f}".format(r[7], r[8], r[9])))
    L.append("| 55 mm past the exit (asymptotic, transmitted) | - | - | - | {} | {:.1f} / {:.1f} / {:.1f} |".format(
        run["n_trans"], *run["asym"]))
    L.append("")
with open(a.out + ".md", "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
for run in runs:
    rows = run["rows"]
    for ax, col, colt, name in ((axes[0], 5, 9, "4D (4th root of det)"), (axes[1], 3, 7, "eps_u (horizontal transverse)"),
                                (axes[2], 4, 8, "eps_w (vertical after the bend)")):
        ln, = ax.plot(rows[:, 1], rows[:, col], "-", lw=1.2, label="{} live".format(run["label"]))
        ax.plot(rows[:, 1], rows[:, colt], "--", lw=1.2, color=ln.get_color(), label="{} transmitted".format(run["label"]))
        ax.set_title(name)
        ax.set_xlabel("mean z of the live core particles [mm, deck frame]")
        ax.set_ylabel("rms emittance [mm mrad]")
        ax.grid(alpha=0.3)
for ax, colt in zip(axes, (9, 7, 8)):
    top = np.nanmax([np.nanmax(run["rows"][:, colt]) for run in runs if np.isfinite(run["rows"][:, colt]).any()] or [1.0])
    ax.set_ylim(0, 1.15 * top)          # the depleted live set at the exit would dominate the axis otherwise
    for run in runs[:1]:
        for name, zq in run["quads"].items():
            ax.axvline(zq, color="0.6", ls=":", lw=1)
            ax.text(zq, ax.get_ylim()[1] * 0.97, name, rotation=90, va="top", ha="right", fontsize=8, color="0.4")
        if run["entrance_mm"] is not None:
            ax.axvline(run["entrance_mm"], color="0.3", ls="--", lw=1)
            ax.text(run["entrance_mm"], ax.get_ylim()[1] * 0.97, "entrance", rotation=90, va="top", ha="right", fontsize=8, color="0.3")
    ax.invert_xaxis()
axes[0].legend(fontsize=8)
fig.suptitle(a.title)
fig.tight_layout()
fig.savefig(a.out + ".png", dpi=130)
print("wrote", a.out + ".md", a.out + ".png")
for run in runs:
    r0, r1 = run["rows"][0], run["rows"][-1]
    print("{}: start 4D {:.1f} -> last live {:.1f} (z {:+.0f} mm); 55 mm asymptotic 4D {:.1f} mm mrad".format(
        run["label"], r0[5], r1[5], r1[1], run["asym"][2]))
