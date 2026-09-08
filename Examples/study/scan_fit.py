"""Pick the optimum of a BunchScan2 quad grid from a quadratic surrogate instead of the noisy argmax.

For the (phi, alpha1, alpha2, vscale) combination that holds the best point, fit
T(q1, q2) = c0 + c1 q1 + c2 q2 + c3 q1^2 + c4 q1 q2 + c5 q2^2 to all grid points of that
combination (least squares over >= 6 points). If the fit is concave and its maximum lies
inside the grid (a quarter step of margin), print that maximum, rounded to 25 V, with the
fitted transmission; otherwise fall back to the argmax. Same output line as scan_best.py:
'phi alpha1 alpha2 q1 q2 vscale transmission', plus a second line 'source: fit|argmax ...'.
"""
import json
import sys

import numpy as np

d = json.load(open(sys.argv[1]))
res, best = d["results"], d["best"]
key = (best["phi"], best["alpha1"], best["alpha2"], best["vscale"])
pts = [r for r in res if (r["phi"], r["alpha1"], r["alpha2"], r["vscale"]) == key]
q1 = np.array([r["q1"] for r in pts])
q2 = np.array([r["q2"] for r in pts])
t = np.array([r["transmission"] for r in pts])
out_q1, out_q2, out_t, source = best["q1"], best["q2"], best["transmission"], "argmax"

if len(pts) >= 6 and len(set(q1)) >= 3 and len(set(q2)) >= 3:
    # scale to unit grid for conditioning
    s1, s2 = q1.mean(), q2.mean()
    w1, w2 = max(np.ptp(q1), 1.0), max(np.ptp(q2), 1.0)
    u, v = (q1 - s1) / w1, (q2 - s2) / w2
    A = np.column_stack([np.ones_like(u), u, v, u * u, u * v, v * v])
    c, *_ = np.linalg.lstsq(A, t, rcond=None)
    H = np.array([[2 * c[3], c[4]], [c[4], 2 * c[5]]])
    eig = np.linalg.eigvalsh(H)
    if np.all(eig < 0):
        uv = np.linalg.solve(H, -np.array([c[1], c[2]]))
        fq1, fq2 = s1 + uv[0] * w1, s2 + uv[1] * w2
        step1 = np.min(np.diff(np.unique(q1))) if len(set(q1)) > 1 else 0.0
        step2 = np.min(np.diff(np.unique(q2))) if len(set(q2)) > 1 else 0.0
        inside = (q1.min() - 0.25 * step1 <= fq1 <= q1.max() + 0.25 * step1 and
                  q2.min() - 0.25 * step2 <= fq2 <= q2.max() + 0.25 * step2)
        resid = float(np.sqrt(np.mean((A @ c - t) ** 2)))
        MAX_RESID = 0.015   # binomial noise of 1,500 particles is 1.1 points; a worse fit means the grid is not quadratic
        if inside and resid <= MAX_RESID:
            ft = float(c[0] + c[1] * uv[0] + c[2] * uv[1] + c[3] * uv[0] ** 2 + c[4] * uv[0] * uv[1] + c[5] * uv[1] ** 2)
            out_q1, out_q2, out_t, source = 25.0 * round(fq1 / 25.0), 25.0 * round(fq2 / 25.0), ft, \
                "fit (rms residual {:.1f} points, argmax {:+.0f}/{:+.0f} V at {:.1f} %)".format(
                    100 * resid, best["q1"], best["q2"], 100 * best["transmission"])
        elif not inside:
            source = "argmax (fit maximum outside the grid: {:+.0f}/{:+.0f} V)".format(fq1, fq2)
        else:
            source = "argmax (fit residual {:.1f} points exceeds the noise)".format(100 * resid)
    else:
        source = "argmax (fit not concave)"

print("{:.1f} {:.1f} {:.1f} {:.0f} {:.0f} {:.4f} {:.4f}".format(key[0], key[1], key[2], out_q1, out_q2, key[3], out_t))
print("source: {}".format(source))
