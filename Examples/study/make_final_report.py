"""Morning report of the night chain v2 (RunNight2.ps1): Results/final/report.md + report.png.

Reads whatever exists (candidate summaries, stage scans, exit metrics, the 44k runs) and
skips the rest, so it can be run at any time during the chain.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DECK = os.environ.get("SI_DECK") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
FIN = os.path.join(DECK, "Results", "final")
WIG = os.path.join(DECK, "Results", "wiggle")
CANDIDATES = ["base", "bore16", "bore16_slot19", "bore18_slot19"]
_ap = argparse.ArgumentParser()
_ap.add_argument("--n-traj", type=int, default=300, help="trajectories drawn in the 3D geometry figures")
_ap.add_argument("--fin", default="final", help="results folder under Results/")
_ap.add_argument("--run1", default="final1", help="Phase B geometry run name under Results/wiggle/")
_ap.add_argument("--run2", default="final2", help="Phase D geometry run name under Results/wiggle/")
_ap.add_argument("--title", default="night chain v2 (2026-09-06/07)")
ARGS = _ap.parse_args()
FIN = os.path.join(DECK, "Results", ARGS.fin)
RUN1, RUN2 = ARGS.run1, ARGS.run2
SPIRAL = {"SI_Anode", "SI_Cathode"}
QUADS = {"D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7"}


def jload(path):
    try:
        with open(path, encoding="utf-8-sig") as fh:
            return json.load(fh)
    except Exception:  # noqa: BLE001
        return None


def pct(x):
    return "{:.1f}".format(100 * x) if x is not None else "-"


def loss_split(losses, n):
    sp = sum(c for k, c in losses.items() if k in SPIRAL)
    qd = sum(c for k, c in losses.items() if k in QUADS)
    ap = sum(losses.values()) - sp - qd
    return sp / n, qd / n, ap / n


lines = ["# HCHC-60 spiral inflector: {}".format(ARGS.title), ""]
lines += ["Beam: Bevatech RFQ exit core (`ext3_exit_core_as_txt.txt`, 43,969 particles, 68.45 keV mean, "
          "1.96 % rms energy spread; space charge at 8 mA on this core). Field: `HCHC-60_CentralBField_z-40to5cm_1mm.pickle`. "
          "Goal: transmission > 85 %, minimal losses on powered electrodes, bunch centroid on the median plane.", ""]

# ------------------------------------------------------------------ candidates (5 mm pipeline)
rows = []
for name in CANDIDATES:
    s = jload(os.path.join(WIG, name, "summary.json"))
    if not s or "bunch" not in s:
        continue
    b, k = s["bunch"], s["knobs"]
    n = b["n"]
    sp, qd, ap = loss_split(b["losses"], n)
    rows.append((name, k.get("quad_bore", 0) * 1e3, k.get("aper_hole", 0) * 1e3, k.get("slot_width", 0.015) * 1e3,
                 s["inner_scan"]["best"]["q1"], s["inner_scan"]["best"]["q2"], b["transmission"], sp, qd, ap))
if rows:
    lines += ["## 1. Geometry candidates (5 mm fields, 4x4 quad grid, 5k bunch; ranking only, 5 mm fields read "
              "~18 points low against 2.5 mm on the same geometry)", "",
              "| candidate | bore [mm] | hole r [mm] | slot [mm] | q1/q2 [V] | transmission [%] | spiral | quad el. | apertures |",
              "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append("| {} | {:.0f} | {:.1f} | {:.0f} | {:+.0f}/{:+.0f} | **{}** | {} | {} | {} |".format(
            r[0], r[1], r[2], r[3], r[4], r[5], pct(r[6]), pct(r[7]), pct(r[8]), pct(r[9])))
    lines.append("")

# ------------------------------------------------------------------ chosen geometry and stages
sc_set = jload(os.path.join(FIN, "stageC_settings.json"))
fin_set = jload(os.path.join(FIN, "final_settings.json"))
f1 = jload(os.path.join(WIG, RUN1, "summary.json"))
f2 = jload(os.path.join(WIG, RUN2, "summary.json"))
if f1:
    o = f1.get("optimizer", {})
    lines += ["## 2. Chosen geometry: `{}`".format(f1["name"] if not sc_set else sc_set["cand"]), "",
              "Knobs: {}".format({k: v for k, v in f1["knobs"].items()}),
              "", "Design energy {:.2f} keV (beam mean). Optimizer at 2.5 mm: spiral {:.0f} V, dz {:+.2f} mm, "
              "residuals {}.".format(1e3 * f1.get("energy_mev", 0), o.get("voltage_V", 0), o.get("dz_mm", 0),
                                     {k: round(v, 3) for k, v in o.get("residual_final", {}).items()}), ""]

stage_rows = []


def add_stage(label, T, sp, qd, ap, note=""):
    stage_rows.append((label, T, sp, qd, ap, note))


C1_OVERRIDE = os.path.join(WIG, RUN1, "c1_best_override.txt")
C1_PER_ANGLE = {80.0: 76.1, 85.0: 76.9, 90.0: 76.7}   # best of the 4x4 quad grid per angle, from the log at 22:12 (95/100 deg incomplete then)
for tag, label in (("inner", "B: inner 3x3 grid (final1)"), ("c1_angle", "C1: beam angle"), ("c2_alpha", "C2: quad rotations"),
                   ("c2_alpha_all", "C2x: quad rotations, extended range (merged)"),
                   ("c3_fine", "C3: 375 V grid"), ("c4_vscale", "C4: spiral voltage")):
    if tag == "c1_angle" and os.path.exists(C1_OVERRIDE):
        f = open(C1_OVERRIDE).readline().split()
        if f[3] == "B":
            add_stage("C1: beam angle scan skipped (angle {:.0f} deg taken over; it was flat within a point over 80-90 deg on 2026-09-06)".format(float(f[0])),
                      None, None, None, None, "quad voltages from the Phase B grid")
        else:
            add_stage("C1: beam angle (result from the log; the per-point file was overwritten by a relaunch)",
                      float(f[6]), None, None, None, "phi {:.0f}, q {:+.0f}/{:+.0f} V; angle flat 76-77 % over 80-90 deg".format(float(f[0]), float(f[3]), float(f[4])))
        continue
    d = jload(os.path.join(WIG, RUN1, "scan2_{}.json".format(tag)))
    if d and d.get("best"):
        b = d["best"]
        add_stage(label, b["transmission"], b["lost_spiral"], b["lost_quads"], b["lost_apertures"],
                  "phi {:.0f}, alpha {:.0f}/{:.0f}, q {:+.0f}/{:+.0f} V, V x{:.3f} ({} particles)".format(
                      b["phi"], b["alpha1"], b["alpha2"], b["q1"], b["q2"], b["vscale"], d["n"]))
metrics = {}
for key, label in (("c5", "C5: 10k bunch, superposed field, at spiral scale 1.03 (before the voltage decision)"),
                   ("c6_vs1p03", "C6b: 10k bunch, direct solve at spiral scale 1.03 (rejected: centroid tilted -1.4 deg asymptotically)"),
                   ("c6", "C6: 10k bunch, direct solve 2.5 mm, spiral scale 1.00"),
                   ("c7", "C7: 10k bunch, direct solve 1.25 mm"),
                   ("d", "D: 10k bunch after the dz correction"),
                   ("e_nosc", "E: all 43,969 particles, no space charge"), ("e_sc", "E: all 43,969 particles, PyAMG space charge 8 mA")):
    m = jload(os.path.join(FIN, "exit_metrics_{}.json".format(key)))
    if m:
        metrics[key] = m
        add_stage(label, m["transmission"], m["lost_spiral"], m["lost_quads"], m["lost_apertures"],
                  "z {:+.2f} mm, vert. {:+.2f} deg".format(m["z_exit_mean_mm"], m["vert_angle_mean_deg"]) if "z_exit_mean_mm" in m else "")
if stage_rows:
    lines += ["## 3. Transmission by stage", "",
              "Transmission counts particles through the housing exit opening (a 19 x 50 mm slot aligned with the tilted "
              "exit electrodes; particles that clear the electrodes but hit it are in the aperture column as `Housing_exit`). "
              "The candidate table and the Phase B inner grid were taken at the exit plane, before that check existed, "
              "and read higher by the housing clip (12 % on final1 before tuning).", "",
              "| stage | transmission [%] | spiral | quad el. | apertures | settings / centroid |", "|---|---|---|---|---|---|"]
    for r in stage_rows:
        lines.append("| {} | **{}** | {} | {} | {} | {} |".format(r[0], pct(r[1]), pct(r[2]), pct(r[3]), pct(r[4]), r[5]))
    lines.append("")

if fin_set:
    lines += ["## 4. Final settings", "",
              "Beam orientation: {}phi {:.0f} deg. Quad rotations {:.0f}/{:.0f} deg. Quad voltages {:+.0f}/{:+.0f} V "
              "(D0,D1 = +q1, D2,D3 = -q1, D4,D5 = +q2, D6,D7 = -q2). Spiral {:.0f} V (scale {:.3f}). Inflector shift dz {:+.2f} mm. "
              "Fields at {:.2f} mm. STEP files: `{}`.".format(
                  "x/y swapped, " if fin_set.get("swap_xy") else "", fin_set["phi"], fin_set["alpha1"], fin_set["alpha2"],
                  fin_set["q1"], fin_set["q2"], fin_set["spiral_V"], fin_set["vscale"], fin_set["dz2_mm"],
                  1e3 * fin_set.get("res_m", 0.0025), fin_set["steps"]), ""]

# ------------------------------------------------------------------ exit centroid
cent = [(k, metrics[k]) for k in ("c6_vs1p03", "c6", "d", "e_nosc", "e_sc") if k in metrics and "z_exit_mean_mm" in metrics[k]]
if cent:
    lines += ["## 5. Exit centroid (targets z = 0, vertical angle = 0)", "",
              "The asymptotic columns are taken about 31 mm past the crossing, outside the electric fringe, and are the values "
              "the central region sees; the crossing columns sit inside the fringe (the design particle itself shows about "
              "2 deg there and 0 deg outside).", "",
              "| run | asymptotic z [mm] | asymptotic vert. angle [deg] | z at crossing [mm] | vert. angle at crossing [deg] | r [mm] (design) | azimuth [deg] (design) | horiz. angle [deg] | pr/p (design) | energy at plane rms [%] |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for k, m in cent:
        za = "{:+.2f} +- {:.2f}".format(m["z_asym_mean_mm"], m["z_asym_rms_mm"]) if "z_asym_mean_mm" in m else "-"
        va = "{:+.2f} +- {:.2f}".format(m["vert_angle_asym_mean_deg"], m["vert_angle_asym_rms_deg"]) if "z_asym_mean_mm" in m else "-"
        lines.append("| {} | **{}** | **{}** | {:+.2f} +- {:.2f} | {:+.2f} +- {:.2f} | {:.2f} +- {:.2f} ({:.2f}) | {:.2f} +- {:.2f} ({:.2f}) | {:+.2f} +- {:.2f} | {:+.3f} ({:+.3f}) | {:.2f} |".format(
            k, za, va, m["z_exit_mean_mm"], m["z_exit_rms_mm"], m["vert_angle_mean_deg"], m["vert_angle_rms_deg"],
            m["r_exit_mean_mm"], m["r_exit_rms_mm"], m["design_r_mm"], m["azimuth_mean_deg"], m["azimuth_rms_deg"], m["design_azimuth_deg"],
            m["horiz_angle_mean_deg"], m["horiz_angle_rms_deg"], m["pr_over_p_mean"], m["design_pr_over_p"], m["energy_at_plane_rms_pct"]))
    lines += ["", "Jarrett's reference (X 59.70 mm, Y 0.70 mm, direction 57.15 deg) is not the same point as this exit plane "
              "(design orbit leaves the electrodes at r = 71.0 mm, velocity 75 deg from the radius vector); the mapping needs "
              "the definition of his reference point.", ""]

# ------------------------------------------------------------------ losses per electrode, final runs
fin_runs = [(k, metrics[k]) for k in ("e_nosc", "e_sc", "d") if k in metrics]
if fin_runs:
    names = sorted(set(n for _, m in fin_runs for n in m["losses"]), key=lambda n: -max(m["losses"].get(n, 0) for _, m in fin_runs))
    lines += ["## 6. Losses per electrode [%]", "", "| electrode | " + " | ".join(k for k, _ in fin_runs) + " |", "|---|" + "---|" * len(fin_runs)]
    for n in names:
        lines.append("| {} | ".format(n) + " | ".join("{:.2f}".format(100 * m["losses"].get(n, 0) / m["n"]) for _, m in fin_runs) + " |")
    lines.append("")

sc_json = jload(os.path.join(FIN, "bunch_e_sc_all.json"))
if sc_json and sc_json.get("sc"):
    s = sc_json["sc"]
    lines += ["Space charge: {} Poisson solves, phi {:+.0f}..{:+.0f} V on the mesh, |E_sc| max on the beam {:.1e} V/m, "
              "wall {:.0f} min.".format(s.get("n_solves"), s.get("phi_min_V", 0), s.get("phi_max_V", 0),
                                       s.get("e_sc_max_on_beam_V_per_m", 0), sc_json.get("wall_time_s", 0) / 60), ""]

# ------------------------------------------------------------------ geometry with trajectories (isometric, -z up)
geo_figs = []
if fin_set and os.path.isdir(fin_set.get("steps", "")):
    from plot_geometry_trajectories import electrode_meshes, plot_geometry_trajectories, plot_side_view
    meshes = None
    for key, npz, label in (("e_nosc", os.path.join(WIG, RUN2, "bunch_e_nosc_all.npz"), "no space charge"),
                            ("e_sc", os.path.join(FIN, "bunch_e_sc_all.npz"), "PyAMG space charge, 8 mA"),
                            ("d", os.path.join(WIG, RUN2, "bunch_d_bunch.npz"), "10k bunch, no space charge (Phase D)")):
        if not os.path.exists(npz) or (key == "d" and "e_nosc" in metrics):
            continue
        try:
            if meshes is None:
                meshes = electrode_meshes(fin_set["steps"])
            m = metrics.get(key, {})
            ttl = "{} : {:,d} particles, transmission {:.1f} %".format(label, m.get("n", 0), 100 * m.get("transmission", 0)) if m else label
            out = os.path.join(FIN, "geometry_trajectories_{}.png".format(key))
            n_ok, n_lost = plot_geometry_trajectories(npz, fin_set["steps"], out, title=ttl, n_traj=ARGS.n_traj, meshes=meshes)
            side = os.path.join(FIN, "geometry_side_{}.png".format(key))
            plot_side_view(npz, fin_set["steps"], side, title=ttl, n_traj=ARGS.n_traj, meshes=meshes)
            geo_figs.append((key, label, os.path.basename(out), n_ok, n_lost, os.path.basename(side)))
        except Exception as exc:  # noqa: BLE001
            print("geometry figure for {} failed: {}".format(key, exc))
if geo_figs:
    lines += ["## 7. Geometry with trajectories (isometric, -z upwards: the beam enters from the top)", ""]
    for key, label, fn, n_ok, n_lost, side in geo_figs:
        lines += ["**{}** ({} transmitted and {} lost trajectories drawn, an evenly spaced subsample of the run, "
                  "followed past the exit plane):".format(label, n_ok, n_lost),
                  "", "![{}]({})".format(label, fn), "",
                  "Side view, line of sight perpendicular to the exit velocity, median plane drawn:", "",
                  "![{} side view]({})".format(label, side), ""]

if fin_set and fin_set.get("spiral_V", 0) > 12000:
    lines += ["**WARNING: the spiral electrode voltage {:.0f} V exceeds the 12 kV sparking limit.**".format(fin_set["spiral_V"]), ""]
lines += ["## 8. Caveats", "",
          "- Candidate ranking used 5 mm fields; everything from Phase B on uses 2.5 mm (or 1.25 mm if the C7 check demanded it).",
          "- Scan points use 1,500 particles (+-1.1 points), bunches 10k (+-0.4) or all 43,969 core particles.",
          "- Quad rotations enter the scans by superposition of normal and skew basis fields; the direct solves confirm them.",
          "- The 4,876 unaccelerated RFQ particles are not in the bunch (no usable longitudinal position); they are 99 % lost before the inflector.",
          "- No central-region constraint was applied: r, azimuth and direction are reported against the design orbit of each geometry.",
          "- The spiral voltage scale 1.03 (C4 optimum, 94.1 % on 1,500 particles) was rejected because it tilts the asymptotic centroid "
          "by -1.4 deg; at scale 1.00 the 10k transmission is the same within noise (92.3 vs 92.4 %) and the centroid is level, at the "
          "price of 1.6 points more on the spiral electrodes.",
          "- The C1 per-point results were lost to a file overwrite; the row above carries the result from the log.",
          "- Quad voltages reach +8.0/-10.8 kV at 18 mm bore (about 12 kV/cm at the hyperbola tips); the housing exit opening is the "
          "entrance-aperture rectangle (40 x 19 mm) rotated with the tilted exit electrodes, unchanged by design decision.",
          "- The pipeline runs at 2.5 mm field resolution; the 1.25 mm check agrees to 0.02 points.", ""]

os.makedirs(FIN, exist_ok=True)
with open(os.path.join(FIN, "report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines))

# ------------------------------------------------------------------ figure
fig, axs = plt.subplots(2, 2, figsize=(14, 9))
ax = axs[0, 0]
if stage_rows:
    labels = [r[0].split(":")[0] for r in stage_rows]
    ax.bar(range(len(stage_rows)), [100 * r[1] if r[1] is not None else 0 for r in stage_rows], color="tab:blue")
    ax.set_xticks(range(len(stage_rows)))
    ax.set_xticklabels(labels, fontsize=8, rotation=0)
    ax.axhline(85, color="tab:red", ls="--", lw=1, label="goal 85 %")
    ax.set_ylabel("transmission [%]")
    ax.set_title("transmission by stage")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
ax = axs[0, 1]
d = jload(os.path.join(WIG, RUN1, "scan2_c1_angle.json"))
if d and len(d["results"]) >= 80:
    by = {}
    for r in d["results"]:
        if r["phi"] not in by or r["transmission"] > by[r["phi"]]["transmission"]:
            by[r["phi"]] = r
    ph = sorted(by)
    ax.plot(ph, [100 * by[p]["transmission"] for p in ph], "o-")
    ax.set_title("C1: beam angle")
else:
    ph = sorted(C1_PER_ANGLE)
    ax.plot(ph, [C1_PER_ANGLE[p] for p in ph], "o-")
    ax.set_title("C1: beam angle (reconstructed from the log; per-point file lost)")
ax.set_xlabel("beam angle phi [deg]")
ax.set_ylabel("best transmission over the quad grid [%]")
ax.grid(alpha=0.3)
ax = axs[1, 0]
d = jload(os.path.join(WIG, RUN1, "scan2_c3_fine.json"))
if d:
    q1s = sorted(set(r["q1"] for r in d["results"]))
    q2s = sorted(set(r["q2"] for r in d["results"]))
    grid = np.full((len(q1s), len(q2s)), np.nan)
    for r in d["results"]:
        grid[q1s.index(r["q1"]), q2s.index(r["q2"])] = 100 * r["transmission"]
    im = ax.imshow(grid, origin="lower", aspect="auto", extent=[q2s[0], q2s[-1], q1s[0], q1s[-1]], cmap="viridis")
    fig.colorbar(im, ax=ax, label="transmission [%]")
    ax.set_xlabel("quad 2 [V]")
    ax.set_ylabel("quad 1 [V]")
    ax.set_title("C3: quad grid, alpha {:.0f}/{:.0f}".format(d["best"]["alpha1"], d["best"]["alpha2"]))
ax = axs[1, 1]
for key, label, col in (("e_nosc", "no space charge", "tab:blue"), ("e_sc", "space charge 8 mA", "tab:orange"),
                        ("d", "D: 10k, no SC", "tab:gray")):
    npz = {"e_nosc": os.path.join(WIG, RUN2, "bunch_e_nosc_all.npz"), "e_sc": os.path.join(FIN, "bunch_e_sc_all.npz"),
           "d": os.path.join(WIG, RUN2, "bunch_d_bunch.npz")}[key]
    if os.path.exists(npz) and (key != "d" or "e_nosc" not in metrics):
        z = 1e3 * np.load(npz, allow_pickle=True)["z_exit"]
        ax.hist(z, bins=np.linspace(-25, 25, 101), histtype="step", color=col, density=True,
                label="{}: mean {:+.2f} mm, rms {:.2f}".format(label, z.mean(), z.std()))
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("exit height z [mm]")
ax.set_ylabel("density")
ax.set_title("exit height distribution")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.suptitle("HCHC-60 inflector, night chain v2: {}".format(fin_set["cand"] if fin_set else (sc_set["cand"] if sc_set else "in progress")))
fig.tight_layout()
fig.savefig(os.path.join(FIN, "report.png"), dpi=140)
print("wrote {} and report.png ({} stage rows, {} candidates)".format(os.path.join(FIN, "report.md"), len(stage_rows), len(rows)))
