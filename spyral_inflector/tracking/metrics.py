"""Exit-beam metrics of a bunch run (bunch_<tag>.npz + the si_state pickle): transmission,
losses by category, energy at the exit plane (inside the fringe potential, so not the
asymptotic spread), vertical and horizontal exit angles, exit height, radius, azimuth and
radial momentum fraction of the crossing points against the design orbit, the asymptotic
(fringe-free) height and vertical angle, and the housing-exit clip."""
import json
import sys

import numpy as np

from .deck import CLIGHT, SPIRAL_NAMES, QUAD_NAMES, load_state

MASS_KEV = 1877.0e3


def exit_metrics(npz_path, state_path):
    d = np.load(npz_path, allow_pickle=True) if isinstance(npz_path, str) else npz_path
    st = load_state(state_path)
    trj, vd = st["trj_design"], st["v_design"]
    n_exit = vd[-1] / np.linalg.norm(vd[-1])
    crossed = d["crossed"]
    n_all = len(crossed)
    ex = d["exit_state"]
    r, v = ex[:, :3], ex[:, 3:]
    names = list(d["electrode_names"])
    hit = d["hit_electrode"]
    losses = {n: int((hit == i).sum()) for i, n in enumerate(names) if (hit == i).any()}
    lost_spiral = sum(c for n, c in losses.items() if n in SPIRAL_NAMES)
    lost_quads = sum(c for n, c in losses.items() if n in QUAD_NAMES)
    lost_ap = int((hit >= 0).sum()) - lost_spiral - lost_quads
    out = {"n": int(n_all), "transmission": float(crossed.mean()),
           "lost_spiral": lost_spiral / n_all, "lost_quads": lost_quads / n_all, "lost_apertures": lost_ap / n_all,
           "lost_anode": losses.get("SI_Anode", 0) / n_all, "lost_cathode": losses.get("SI_Cathode", 0) / n_all,
           "losses": losses}
    if len(r) == 0:
        return out

    def ekin(vv):
        b = np.linalg.norm(vv, axis=1) / CLIGHT
        return MASS_KEV * (1 / np.sqrt(1 - b ** 2) - 1)
    e_in, e_out = ekin(d["v0"][crossed]), ekin(v)
    dz = r[:, 2] - trj[-1][2]
    vert = np.degrees(np.arcsin(v[:, 2] / np.linalg.norm(v, axis=1)))
    rad = np.hypot(r[:, 0], r[:, 1])
    az = np.degrees(np.arctan2(r[:, 1], r[:, 0]))
    vr = (v[:, 0] * r[:, 0] + v[:, 1] * r[:, 1]) / rad / np.linalg.norm(v, axis=1)
    vh = v[:, :2] / np.linalg.norm(v[:, :2], axis=1)[:, None]
    hang = np.degrees(np.arctan2(vh[:, 0] * n_exit[1] - vh[:, 1] * n_exit[0], vh @ n_exit[:2]))
    out.update({
        "energy_in_mean_keV": float(e_in.mean()), "energy_in_rms_pct": float(100 * e_in.std() / e_in.mean()),
        "energy_at_plane_rms_pct": float(100 * e_out.std() / e_out.mean()),
        "vert_angle_mean_deg": float(vert.mean()), "vert_angle_rms_deg": float(vert.std()),
        "z_exit_mean_mm": float(1e3 * dz.mean()), "z_exit_rms_mm": float(1e3 * dz.std()),
        "horiz_angle_mean_deg": float(hang.mean()), "horiz_angle_rms_deg": float(hang.std()),
        "r_exit_mean_mm": float(1e3 * rad.mean()), "r_exit_rms_mm": float(1e3 * rad.std()),
        "azimuth_mean_deg": float(az.mean()), "azimuth_rms_deg": float(az.std()),
        "pr_over_p_mean": float(vr.mean()), "pr_over_p_rms": float(vr.std()),
        "design_r_mm": float(1e3 * np.hypot(*trj[-1][:2])),
        "design_azimuth_deg": float(np.degrees(np.arctan2(trj[-1][1], trj[-1][0]))),
        "design_pr_over_p": float((n_exit[0] * trj[-1][0] + n_exit[1] * trj[-1][1]) / np.hypot(*trj[-1][:2])),
    })
    # particles that cleared the exit plane but hit the housing exit opening on the way out
    if "post_exit_hit" in d.files:
        post = d["post_exit_hit"]
        out["lost_housing_after_exit"] = float(post.sum()) / n_all
        out["transmission_through_housing"] = float((crossed & ~post).sum()) / n_all
        out["transmission_exit_plane"] = float((crossed | post).sum()) / n_all
    # asymptotic state: a fixed number of steps after the crossing, outside the electric fringe
    if "asym_state" in d.files:
        a = d["asym_state"]
        ok = crossed & np.all(np.isfinite(a), axis=1)
        if ok.any():
            ra, va = a[ok, :3], a[ok, 3:]
            vert_a = np.degrees(np.arcsin(va[:, 2] / np.linalg.norm(va, axis=1)))
            path_mm = 1e3 * float(d["asym_steps"]) * float(d["dt"]) * float(np.linalg.norm(va, axis=1).mean())
            out.update({
                "asym_path_mm": path_mm, "asym_n": int(ok.sum()),
                "z_asym_mean_mm": float(1e3 * ra[:, 2].mean()), "z_asym_rms_mm": float(1e3 * ra[:, 2].std()),
                "vert_angle_asym_mean_deg": float(vert_a.mean()), "vert_angle_asym_rms_deg": float(vert_a.std()),
                "r_asym_mean_mm": float(1e3 * np.hypot(ra[:, 0], ra[:, 1]).mean()),
                "azimuth_asym_mean_deg": float(np.degrees(np.arctan2(ra[:, 1], ra[:, 0])).mean()),
            })
    return out


def fmt(m):
    if "vert_angle_rms_deg" not in m:
        return "transmission {:.1f} %, nothing transmitted".format(100 * m["transmission"])
    return ("transmission {:.1f} % (spiral {:.1f}, quads {:.1f}, apertures {:.1f} %); exit: vert. angle {:+.2f} +- {:.2f} deg, "
            "z {:+.2f} +- {:.2f} mm, horiz. angle {:+.2f} +- {:.2f} deg, r {:.2f} +- {:.2f} mm (design {:.2f}), "
            "azimuth {:.2f} +- {:.2f} deg (design {:.2f}), pr/p {:+.3f} +- {:.3f} (design {:+.3f})").format(
        100 * m["transmission"], 100 * m["lost_spiral"], 100 * m["lost_quads"], 100 * m["lost_apertures"],
        m["vert_angle_mean_deg"], m["vert_angle_rms_deg"], m["z_exit_mean_mm"], m["z_exit_rms_mm"],
        m["horiz_angle_mean_deg"], m["horiz_angle_rms_deg"], m["r_exit_mean_mm"], m["r_exit_rms_mm"], m["design_r_mm"],
        m["azimuth_mean_deg"], m["azimuth_rms_deg"], m["design_azimuth_deg"], m["pr_over_p_mean"], m["pr_over_p_rms"],
        m["design_pr_over_p"]) + (
        "; asymptotic ({:.0f} mm past the crossing, {} particles): z {:+.2f} +- {:.2f} mm, vert. angle {:+.2f} +- {:.2f} deg, r {:.1f} mm, azimuth {:.1f} deg".format(
            m["asym_path_mm"], m["asym_n"], m["z_asym_mean_mm"], m["z_asym_rms_mm"], m["vert_angle_asym_mean_deg"],
            m["vert_angle_asym_rms_deg"], m["r_asym_mean_mm"], m["azimuth_asym_mean_deg"]) if "z_asym_mean_mm" in m else "") + (
        "; exit plane {:.1f} %, housing exit opening {:.2f} % -> through the housing {:.1f} %".format(
            100 * m["transmission_exit_plane"], 100 * m["lost_housing_after_exit"], 100 * m["transmission_through_housing"]) if "lost_housing_after_exit" in m else "")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2:
        sys.exit("usage: python -m spyral_inflector.tracking.metrics bunch_<tag>.npz si_state_<tag>.pickle [out.json]")
    m = exit_metrics(argv[0], argv[1])
    print(fmt(m))
    if len(argv) > 2:
        with open(argv[2], "w") as fh:
            json.dump(m, fh, indent=2)
    return m


if __name__ == "__main__":
    main()
