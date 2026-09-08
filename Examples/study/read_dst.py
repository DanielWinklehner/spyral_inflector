"""Read TraceWin .dst particle distributions, print beam statistics, compare with
Jarrett's TenThousandRFQParticles.txt, and optionally export in that text format.

TraceWin .dst layout (little-endian): 2 dummy bytes, int32 number of particles,
float64 beam current [mA], float64 RF frequency [MHz], 1 dummy byte, then per particle
six float64: x [cm], x' [rad], y [cm], y' [rad], phase [rad], kinetic energy [MeV];
the file ends with the rest mass [MeV/c^2]. TraceWin's phase is omega*t relative to
the reference particle, so a positive phase is a late particle: z = -phase * beta*lambda / 2pi.

    python read_dst.py "../Particles/MIT RFQ Beamdynamics/ext3_exit.dst" --compare --export --core
"""
import argparse
import json
import os

import numpy as np

CLIGHT = 299792458.0
DECK = os.environ.get("SI_DECK") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
JARRETT = os.path.join(DECK, "Particles", "TenThousandRFQParticles.txt")


def load_dst(path):
    with open(path, "rb") as fh:
        np.fromfile(fh, dtype=np.uint8, count=2)
        n = int(np.fromfile(fh, dtype=np.int32, count=1)[0])
        current = float(np.fromfile(fh, dtype=np.float64, count=1)[0])
        freq = float(np.fromfile(fh, dtype=np.float64, count=1)[0])
        np.fromfile(fh, dtype=np.uint8, count=1)
        data = np.fromfile(fh, dtype=np.float64, count=6 * n).reshape(n, 6)
        mass = float(np.fromfile(fh, dtype=np.float64, count=1)[0])
        rest = fh.read()
    if rest:
        print("  note: {} trailing bytes after the rest mass".format(len(rest)))
    return {"n": n, "current_mA": current, "freq_MHz": freq, "mass_MeV": mass,
            "x": 1e-2 * data[:, 0], "xp": data[:, 1], "y": 1e-2 * data[:, 2], "yp": data[:, 3],
            "phase_rad": data[:, 4], "energy_MeV": data[:, 5]}


def load_jarrett(path=JARRETT):
    d = np.loadtxt(path, skiprows=1)
    d = d[d[:, 9] == 0]
    # columns: x(mm) x'(mrad) y(mm) y'(mrad) z(mm) z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss
    return {"n": len(d), "x": 1e-3 * d[:, 0], "xp": 1e-3 * d[:, 1], "y": 1e-3 * d[:, 2], "yp": 1e-3 * d[:, 3],
            "z": 1e-3 * d[:, 4], "phase_rad": np.deg2rad(d[:, 6]), "time_s": d[:, 7], "energy_MeV": d[:, 8],
            "mass_MeV": 2 * 938.272 + 0.511, "freq_MHz": 32.8, "current_mA": None}


def twiss(u, up):
    u = u - u.mean()
    up = up - up.mean()
    eps = np.sqrt(np.mean(u * u) * np.mean(up * up) - np.mean(u * up) ** 2)
    return eps, np.mean(u * u) / eps, -np.mean(u * up) / eps, -np.mean(u * up) / np.mean(up * up)


def stats(b, label):
    e = b["energy_MeV"]
    em = e.mean()
    gamma = 1 + em / b["mass_MeV"]
    beta = np.sqrt(1 - 1 / gamma ** 2)
    bl = beta * CLIGHT / (b["freq_MHz"] * 1e6)          # beta*lambda [m]
    ph = np.rad2deg(b["phase_rad"])
    z = -b["phase_rad"] / (2 * np.pi) * bl               # late particle -> negative z
    de = (e - em) / em
    out = {"label": label, "n": int(b["n"]), "current_mA": b["current_mA"], "freq_MHz": b["freq_MHz"],
           "mass_MeV": b["mass_MeV"], "energy_mean_MeV": float(em), "energy_median_MeV": float(np.median(e)),
           "energy_rms_pct": float(100 * e.std() / em),
           "energy_percentiles_MeV": {str(p): float(np.percentile(e, p)) for p in (1, 5, 25, 75, 95, 99)},
           "within_pct": {str(w): float(100 * np.mean(np.abs(de) < w / 100)) for w in (2, 5, 8, 10, 11, 15, 20)},
           "beta_lambda_mm": 1e3 * bl,
           "phase_mean_deg": float(ph.mean()), "phase_rms_deg": float(ph.std()),
           "phase_min_deg": float(ph.min()), "phase_max_deg": float(ph.max()),
           "z_rms_mm": float(1e3 * z.std()), "corr_energy_phase": float(np.corrcoef(ph, de)[0, 1])}
    print("\n=== {}: {} particles, I = {} mA, f = {} MHz, m = {:.3f} MeV".format(
        label, b["n"], b["current_mA"], b["freq_MHz"], b["mass_MeV"]))
    print("  energy: mean {:.5f} MeV, median {:.5f}, rms {:.2f} %, 1/99 % = {:.4f}/{:.4f} MeV".format(
        em, np.median(e), out["energy_rms_pct"], out["energy_percentiles_MeV"]["1"], out["energy_percentiles_MeV"]["99"]))
    print("  within +-2/5/8/10/11/15/20 %: " + "  ".join("{}%".format(round(out["within_pct"][k], 1)) for k in ("2", "5", "8", "10", "11", "15", "20")))
    print("  phase: mean {:+.1f} deg, rms {:.1f} deg, range {:+.0f}..{:+.0f} deg; z rms {:.2f} mm (beta*lambda {:.1f} mm); corr(E, phase) {:+.3f}".format(
        ph.mean(), ph.std(), ph.min(), ph.max(), 1e3 * z.std(), 1e3 * bl, out["corr_energy_phase"]))
    for lab, u, up in (("x", b["x"], b["xp"]), ("y", b["y"], b["yp"])):
        eps, bet, alp, sw = twiss(u, up)
        out["twiss_" + lab] = {"rms_mm": float(1e3 * u.std()), "rms_prime_mrad": float(1e3 * up.std()),
                               "eps_rms_mm_mrad": float(1e6 * eps), "eps_norm_rms_mm_mrad": float(1e6 * eps * beta * gamma),
                               "beta_m": float(bet), "alpha": float(alp), "waist_mm": float(1e3 * sw),
                               "mean_mm": float(1e3 * u.mean()), "mean_prime_mrad": float(1e3 * up.mean())}
        print("  {}: rms {:.2f} mm, rms {}' {:.1f} mrad, eps_rms {:.1f} mm mrad (norm. {:.3f}), beta {:.3f} m, alpha {:+.2f}, waist {:+.1f} mm; centroid {:+.2f} mm / {:+.1f} mrad".format(
            lab, 1e3 * u.std(), lab, 1e3 * up.std(), 1e6 * eps, 1e6 * eps * beta * gamma, bet, alp, 1e3 * sw,
            1e3 * u.mean(), 1e3 * up.mean()))
    # longitudinal rms emittance in deg keV and in mm mrad-like units (z, dp/p)
    dw = 1e3 * (e - em)
    eps_l = np.sqrt(np.mean((ph - ph.mean()) ** 2) * np.mean(dw ** 2) - np.mean((ph - ph.mean()) * dw) ** 2)
    dpp = de / (1 + 1 / gamma) if False else de * (gamma / (gamma + 1))   # dp/p = dT/T * gamma/(gamma+1)
    eps_z = np.sqrt(np.mean((z - z.mean()) ** 2) * np.mean(dpp ** 2) - np.mean((z - z.mean()) * dpp) ** 2)
    out["eps_long_rms_deg_keV"] = float(eps_l)
    out["eps_long_rms_mm_mrad"] = float(1e6 * eps_z)
    print("  longitudinal: eps_rms {:.1f} deg keV = {:.1f} mm mrad (z, dp/p); dp/p rms {:.2f} %".format(
        eps_l, 1e6 * eps_z, 100 * dpp.std()))
    return out, z


def export_like_jarrett(b, z, path):
    """Write the distribution in the column format of TenThousandRFQParticles.txt."""
    e = b["energy_MeV"]
    gamma = 1 + e / b["mass_MeV"]
    beta = np.sqrt(1 - 1 / gamma ** 2)
    p = gamma * beta
    dpp = p / p.mean() - 1
    t = b["phase_rad"] / (2 * np.pi * b["freq_MHz"] * 1e6)
    cols = np.column_stack([1e3 * b["x"], 1e3 * b["xp"], 1e3 * b["y"], 1e3 * b["yp"], 1e3 * z, 1e3 * dpp,
                            np.rad2deg(b["phase_rad"]), t, e, np.zeros(len(e))])
    with open(path, "w") as fh:
        fh.write("x(mm) x'(mrad) y(mm) y'(mrad) z(mm) z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss\n")
        np.savetxt(fh, cols, fmt="%.6e")
    print("  wrote {}".format(path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    parser.add_argument("--compare", action="store_true", help="also print Jarrett's file and make a comparison figure")
    parser.add_argument("--export", action="store_true", help="write <file>_as_txt.txt in Jarrett's format")
    parser.add_argument("--core", action="store_true",
                        help="also report the bunch core: |phase| <= 180 deg and energy >= half the median "
                             "(drops unaccelerated stragglers); exported as <file>_core_as_txt.txt")
    args = parser.parse_args()

    beams = []
    for f in args.files:
        b = load_dst(f)
        s, z = stats(b, os.path.basename(f))
        beams.append(("{}\nall {} particles, {:.2f} mA".format(os.path.basename(f), b["n"], b["current_mA"]), b, z))
        with open(os.path.splitext(f)[0] + "_stats.json", "w") as fh:
            json.dump(s, fh, indent=2)
        if args.export:
            export_like_jarrett(b, z, os.path.splitext(f)[0] + "_as_txt.txt")
        if args.core:
            m = (np.abs(b["phase_rad"]) <= np.pi) & (b["energy_MeV"] >= 0.5 * np.median(b["energy_MeV"]))
            core = {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in b.items()}
            core["n"] = int(m.sum())
            core["current_mA"] = b["current_mA"] * m.mean()
            print("  core selection keeps {} of {} particles ({:.1f} %)".format(m.sum(), b["n"], 100 * m.mean()))
            s, zc = stats(core, os.path.basename(f) + " core")
            if m.mean() < 0.999:
                beams.append(("{}\ncore ({:.0f} %, {} particles)".format(os.path.basename(f), 100 * m.mean(), int(m.sum())), core, zc))
            else:
                print("  core is the whole file; no separate row in the figure")
            with open(os.path.splitext(f)[0] + "_core_stats.json", "w") as fh:
                json.dump(s, fh, indent=2)
            if args.export:
                export_like_jarrett(core, zc, os.path.splitext(f)[0] + "_core_as_txt.txt")
                # the rest: unaccelerated stragglers spread over many RF periods. Written
                # with z = 0 so they can be tracked as a beam of their own (static fields
                # do not care about their arrival time; only space charge would)
                strag = {k: (v[~m] if isinstance(v, np.ndarray) else v) for k, v in b.items()}
                strag["n"] = int((~m).sum())
                strag["current_mA"] = b["current_mA"] * (1 - m.mean())
                export_like_jarrett(strag, np.zeros(strag["n"]), os.path.splitext(f)[0] + "_stragglers_as_txt.txt")
                print("  stragglers: {} particles ({:.1f} %), energies {:.1f}..{:.1f} keV, written with z = 0".format(
                    strag["n"], 100 * (1 - m.mean()), 1e3 * strag["energy_MeV"].min(), 1e3 * strag["energy_MeV"].max()))

    if args.compare:
        j = load_jarrett()
        s, zj = stats(j, "TenThousandRFQParticles.txt (Jarrett)")
        print("  (Jarrett's own z column: rms {:.2f} mm; corr with -phase*beta*lambda/2pi: {:+.4f})".format(
            1e3 * j["z"].std(), np.corrcoef(j["z"], zj)[0, 1]))
        beams.append(("TenThousandRFQParticles.txt\n(Jarrett, {} particles)".format(j["n"]), j, j["z"]))

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        nb = len(beams)
        fig, axs = plt.subplots(nb, 4, figsize=(21, 4.6 * nb), squeeze=False)
        for i, (lab, b, z) in enumerate(beams):
            e = b["energy_MeV"]
            axs[i, 0].annotate(lab, xy=(0, 0.5), xytext=(-62, 0), xycoords="axes fraction", textcoords="offset points",
                               rotation=90, va="center", ha="center", fontsize=11, fontweight="bold")
            for k, (u, v, xl, yl) in enumerate(((1e3 * b["x"], 1e3 * b["xp"], "x [mm]", "x' [mrad]"),
                                                (1e3 * b["y"], 1e3 * b["yp"], "y [mm]", "y' [mrad]"),
                                                (np.rad2deg(b["phase_rad"]), 1e3 * e, "phase [deg]", "energy [keV]"))):
                ax = axs[i, k]
                ax.hist2d(u, v, bins=80, cmap="viridis", cmin=1)
                ax.set_xlabel(xl)
                ax.set_ylabel(yl)
                ax.set_title("{} vs {}".format(yl.split(" ")[0], xl.split(" ")[0]), fontsize=10)
            ax = axs[i, 3]
            ax.hist(1e3 * e, bins=np.arange(5, 100, 1.0), color="C0")
            ax.set_yscale("log")
            ax.axvline(1e3 * np.median(e), color="k", ls="--", lw=0.8)
            ax.set_xlabel("energy [keV]")
            ax.set_ylabel("particles (log)")
            core_m = (np.abs(b["phase_rad"]) <= np.pi) & (e >= 0.5 * np.median(e))
            title = "energy: median {:.2f} keV, rms {:.2f} %".format(1e3 * np.median(e), 100 * e.std() / e.mean())
            if core_m.mean() < 0.999:
                title += "\n(core {:.0f} %: mean {:.2f} keV, rms {:.2f} %; stragglers below)".format(
                    100 * core_m.mean(), 1e3 * e[core_m].mean(), 100 * e[core_m].std() / e[core_m].mean())
            ax.set_title(title, fontsize=9)
        fig.tight_layout(rect=(0.03, 0, 1, 1))
        out = os.path.join(os.path.dirname(os.path.abspath(args.files[0])), "dst_comparison.png")
        fig.savefig(out, dpi=130)
        print("wrote", out)
