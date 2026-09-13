"""Vertical budget of the final1 exit beam (8 mA and vacuum): centroid height and tilt vs rms
size and divergence at the asymptotic state (55 mm past the exit), and a ballistic projection
of how much of the core would be outside a vertical half-aperture after another 0..300 mm
of path -- with and without the centroid removed -- to tell offset/tilt from divergence."""
import os
import sys
import numpy as np

R = sys.argv[1]
print("RUN:", R)
for tag in ("vac_all", "sc_all"):
    d = np.load(os.path.join(R, "bunch_{}.npz".format(tag)), allow_pickle=True)
    st = d["asym_state"]
    ok = d["crossed"] & np.all(np.isfinite(st), axis=1) & ~d["is_tail"]
    z = -1e3 * st[ok, 2]                                   # BCS height [mm] (deck z mirrored)
    v = st[ok, 3:6]
    zp = -1e3 * v[:, 2] / np.hypot(v[:, 0], v[:, 1])       # BCS vertical angle [mrad]
    print("=" * 100)
    print("{}: {} core particles at 55 mm past the exit (Baseline frame, +z up)".format(tag, ok.sum()))
    print("  centroid height {:+.2f} mm, mean vertical angle {:+.1f} mrad ({:+.2f} deg)".format(z.mean(), zp.mean(), np.degrees(1e-3 * zp.mean())))
    print("  rms height {:.2f} mm, rms angle {:.1f} mrad, corr(z, z') {:+.2f}".format(z.std(), zp.std(), np.corrcoef(z, zp)[0, 1]))
    print("  ballistic projection, fraction of the core outside a vertical half-aperture (no CR focusing):")
    print("  {:>8s}  {:>26s}  {:>26s}  {:>26s}".format("s [mm]", "+-5 mm: all / centred", "+-8 mm: all / centred", "+-10 mm: all / centred"))
    for s in (0, 50, 100, 150, 200, 300):
        zz = z + 1e-3 * zp * s                             # mm + mrad * mm / 1000
        zc = (z - z.mean()) + 1e-3 * (zp - zp.mean()) * s
        cells = []
        for a in (5, 8, 10):
            cells.append("{:5.1f} % / {:5.1f} %".format(100 * np.mean(np.abs(zz) > a), 100 * np.mean(np.abs(zc) > a)))
        print("  {:8.0f}  {:>26s}  {:>26s}  {:>26s}".format(s, *cells))

print("=" * 100)
print("transverse emittances at 55 mm past the exit (rms, unnormalized, mm mrad); input core: x 30.1, y 40.2 -> 4D 1210")
for tag in ("vac_all", "sc_all"):
    d = np.load(os.path.join(R, "bunch_{}.npz".format(tag)), allow_pickle=True)
    st = d["asym_state"]
    ok = d["crossed"] & np.all(np.isfinite(st), axis=1) & ~d["is_tail"]
    r, v = st[ok, :3], st[ok, 3:6]
    vm = v.mean(0); t_hat = vm / np.linalg.norm(vm)                       # mean direction
    u_hat = np.cross([0.0, 0.0, 1.0], t_hat); u_hat /= np.linalg.norm(u_hat)   # horizontal transverse
    vlong = v @ t_hat
    u = 1e3 * (r - r.mean(0)) @ u_hat; up = 1e3 * (v @ u_hat) / vlong
    z = 1e3 * r[:, 2]; zp = 1e3 * v[:, 2] / vlong
    X = np.column_stack([u, up, z, zp]); S = np.cov(X.T)
    eu = np.sqrt(np.linalg.det(S[:2, :2])); ez = np.sqrt(np.linalg.det(S[2:, 2:])); e4 = np.sqrt(np.linalg.det(S))
    print("  {:8s} eps_u {:6.1f}  eps_z {:6.1f}  product {:7.0f}  sqrt(det 4D) {:7.0f}  -> coupling ratio product/4D {:.2f}".format(tag, eu, ez, eu * ez, e4, eu * ez / e4))
