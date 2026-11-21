"""Analytic MAGNETIC quadrupoles in the injection line, added to the cyclotron's B map.

Daniel, 2026-10-20 (HCHC-1.5): "do not include the physical models of the quads, instead, use the new analytic
magnetic quad capability to simulate two quads." The lenses are PyPATools' MagneticQuadrupole (Enge fringes,
Maxwell-consistent to the order kept; documents/analytic_quadrupoles.md in PyPATools). Nothing of them is in the
STEP assembly: they only add to the B field, so the electrostatic chain (bases, collisions, space charge) never
sees them, and a grounded beam pipe (geometry knob ``pipe``) stops the beam instead of their poles.

A quad SET is described independently of the system rotation R, exactly as the electrostatic quads were (they
turned with the assembly):

    {"type": "magquad_set",
     "G_T_per_m": [G1, G2],      body gradients, DECK-frame sign: G > 0 focuses x for the beam travelling +z (the
                                 sign the electrostatic voltages had: an ESQ at +V and a magnet at G = 2 V / (a^2 v)
                                 focus the same plane)
     "alpha_deg": [a1, a2],      turn of each quad relative to the inflector assembly [deg]; the tracking frame
                                 sees R + a_i (R = the run's rotation_deg, read from its design state)
     "b_mm": b,                  Baseline z of the centre of the gap between the two irons [mm]
     "q1_top_mm": 253.7, "q2_bottom_mm": 119.7, "gap_mm": 10.0,   iron extents in the Baseline frame [mm]
     "r0_m": 0.021, "leff_per_r0": 1.1, "enge": "cosy_quad", "order": 2, "r_max_m": 0.022}

The iron runs Q1: q1_top .. b + gap/2, Q2: b - gap/2 .. q2_bottom (Baseline mm, the beam travels down). The
effective length is L_eff = L_iron + leff_per_r0 * r0: the PyPATools document gives an iron-to-L_eff relation only
for plate-terminated electrostatic quads, so this is the generic rule for a magnet's fringe extension, and it is an
assumption. The quads are centred on the iron. The iron positions are PRE-levelling Baseline positions, like the ESQ
knobs (quad_z1): the quads move with the assembly's axial shifts exactly as the ESQs did -- the geometry build's lab shift
(the design state's shift_lab, = the build dz) plus a tracking-time --shift-z (Daniel, 2026-10-21: "the analytic magnetic
field has to shift with everything else during the levelling, just like the esqs did"). tracking_spec(..., dz=) takes
that total; the first seed scan and A/B of 2026-10-20 ran with dz = 0, i.e. both magnets 2.90 mm upstream.

The sampled field is added to the B map node by node on the map's own grid, so the result is again one plain
gridded Field on the numba backend and the fused trilinear (GPU) path is unchanged. r_max: the series is meaningless
inside the iron (it grows like r^(2 order + 1)), so the field is zero beyond r_max = r0 + one 1 mm map cell.
"""
import json
import os

import numpy as np

MAGQUAD_DEFAULTS = dict(q1_top_mm=253.7, q2_bottom_mm=119.7, gap_mm=10.0, r0_m=0.021, leff_per_r0=1.1,
                        enge="cosy_quad", order=2, r_max_m=0.022)


def magquad_set(G, alpha_deg, b_mm, **over):
    """A quad-set description (see the module docstring) with the defaults filled in."""
    d = dict(type="magquad_set", G_T_per_m=[float(g) for g in G], alpha_deg=[float(x) for x in alpha_deg],
             b_mm=float(b_mm))
    for k, v in MAGQUAD_DEFAULTS.items():
        d[k] = over.pop(k, v)
    if over:
        raise TypeError("unknown magquad_set keys: {}".format(sorted(over)))
    return d


def load_set(src):
    """A quad set or a ready PyPATools spec from a dict, a JSON file or a JSON string."""
    if isinstance(src, dict):
        return dict(src)
    if os.path.exists(str(src)):
        with open(src) as fh:
            return json.load(fh)
    return json.loads(src)


def iron_layout(desc):
    """[(top, bottom), (top, bottom)]: the Baseline z extents of the two irons [mm], top > bottom."""
    d = dict(MAGQUAD_DEFAULTS, **desc)
    b, g = float(d["b_mm"]), float(d["gap_mm"])
    out = [(float(d["q1_top_mm"]), b + 0.5 * g), (b - 0.5 * g, float(d["q2_bottom_mm"]))]
    for i, (top, bot) in enumerate(out):
        if not top > bot:
            raise ValueError("quad {} has no iron: Baseline {:.2f} .. {:.2f} mm (b = {:.2f} mm)".format(i + 1, top, bot, b))
    return out


def tracking_spec(desc, rotation_deg=0.0, dz=0.0):
    """The PyPATools quadrupoles_from_spec spec of a quad set in the TRACKING frame (deck frame, the inflector turned
    by rotation_deg and shifted along the deck z by dz [m], the build's lab shift plus any tracking shift): z_c =
    -Baseline centre + dz, alpha = R + alpha_i, G as given (deck sign). A dict that is already a spec (has "quads") is
    returned unchanged. The 'layout' block records what it was made from."""
    if "quads" in desc:
        return desc
    d = dict(MAGQUAD_DEFAULTS, **desc)
    r0 = float(d["r0_m"])
    quads, layout = [], []
    for i, (top, bot) in enumerate(iron_layout(d)):
        L_iron = 1e-3 * (top - bot)
        L_eff = L_iron + float(d["leff_per_r0"]) * r0
        zc_bcs = 0.5e-3 * (top + bot)
        alpha = float(rotation_deg) + float(d["alpha_deg"][i])
        quads.append({"name": "Q{}".format(i + 1), "G0": float(d["G_T_per_m"][i]), "L_eff": L_eff, "z_c": -zc_bcs + float(dz),
                      "alpha_deg": alpha})
        layout.append({"quad": "Q{}".format(i + 1), "iron_bcs_mm": [top, bot], "dz_m": float(dz),
                       "iron_bcs_mm_placed": [top - 1e3 * float(dz), bot - 1e3 * float(dz)], "L_iron_mm": 1e3 * L_iron,
                       "L_eff_mm": 1e3 * L_eff, "z_c_bcs_mm": 1e3 * zc_bcs, "G_T_per_m": float(d["G_T_per_m"][i]),
                       "alpha_rel_deg": float(d["alpha_deg"][i]), "alpha_tracking_deg": alpha,
                       "integrated_T": float(d["G_T_per_m"][i]) * L_eff})
    return {"defaults": {"kind": "magnetic", "r0": r0, "order": int(d["order"]), "enge": d["enge"],
                         "r_max": (None if d.get("r_max_m") is None else float(d["r_max_m"]))},
            "quads": quads, "layout": layout, "rotation_deg": float(rotation_deg), "dz_m": float(dz), "source": d}


def baseline_spec(desc, rotation_deg=0.0, dz=0.0):
    """The same quads in the BASELINE frame (for exports), placed with the same deck-frame shift dz: the mirror z -> -z
    keeps x, y and the azimuth, flips the sign of G0 (B is a pseudovector) and mirrors z_c (PyPATools
    test_mirror_frame_rules)."""
    spec = json.loads(json.dumps(tracking_spec(desc, rotation_deg, dz)))
    for q in spec["quads"]:
        q["G0"] = -q["G0"]
        q["z_c"] = -q["z_c"]
    spec["frame"] = ("Baseline (machine) frame: +z up, the beam enters from +z and travels -z; G0 > 0 focuses x for a "
                     "positive ion moving along +z, so these G0 are the deck-frame gradients with the sign flipped")
    return spec


def quads_of(spec):
    from PyPATools.field import quadrupoles_from_spec
    return quadrupoles_from_spec({"defaults": spec.get("defaults", {}), "quads": spec["quads"]})


def add_to_bfield(bfield, spec, log=None):
    """bfield (a gridded Field, tracking frame) plus the quads of spec, sampled on bfield's own nodes. Only the box
    the quads reach (|x|, |y| <= r_max and their profiles' z support) is evaluated; outside it they are exactly zero.
    Returns a new Field on the numba backend; bfield is not modified."""
    from PyPATools.field import Field
    log = log or (lambda m: None)
    quads = quads_of(spec)
    g = bfield.grid
    xs, ys, zs = (np.asarray(g[c], dtype=float) for c in "xyz")
    scale = float(getattr(bfield, "scaling", 1.0) or 1.0)
    vals = {c: scale * np.array(bfield.grid_values[c], dtype=float) for c in "xyz"}
    reach = max((q.r_max if q.r_max is not None else 2.0 * q.r0) for q in quads)
    z_lo = min(q.profile.support()[0] for q in quads)
    z_hi = max(q.profile.support()[1] for q in quads)
    dx = xs[1] - xs[0]
    ix = np.nonzero(np.abs(xs) <= reach + dx)[0]
    iy = np.nonzero(np.abs(ys) <= reach + dx)[0]
    iz = np.nonzero((zs >= z_lo) & (zs <= z_hi))[0]
    if not (ix.size and iy.size and iz.size):
        raise ValueError("the quads do not overlap the B map's grid")
    X, Y, Z = np.meshgrid(xs[ix], ys[iy], zs[iz], indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    add = np.zeros((len(pts), 3))
    for q in quads:
        for i0 in range(0, len(pts), 1 << 20):
            add[i0:i0 + (1 << 20)] += q(pts[i0:i0 + (1 << 20)])
    sl = np.ix_(ix, iy, iz)
    for k, c in enumerate("xyz"):
        vals[c][sl] += add[:, k].reshape(X.shape)
    for q, lay in zip(quads, spec.get("layout") or [{}] * len(quads)):
        log("  magnetic quad  : {} G {:+.3f} T/m, L_eff {:.1f} mm (iron {}), z_c {:+.2f} mm deck, alpha {:+.2f} deg "
            "tracking, r0 {:.1f} mm, r_max {}".format(
                q.label, q.G0, 1e3 * q.L_eff, "Baseline {:.2f}..{:.2f} mm".format(*lay.get("iron_bcs_mm_placed", lay["iron_bcs_mm"]))
                if lay else "-",
                1e3 * q.z_c, q.alpha_deg, 1e3 * q.r0, "-" if q.r_max is None else "{:.1f} mm".format(1e3 * q.r_max)))
    log("  magnetic quads added to the B map on its own grid: {} x {} x {} nodes (z {:+.3f}..{:+.3f} m)".format(
        ix.size, iy.size, iz.size, zs[iz[0]], zs[iz[-1]]))
    return Field.from_arrays(grid={c: np.asarray(g[c], dtype=float) for c in "xyz"}, values=vals,
                             label="{} + {} magnetic quads".format(getattr(bfield, "label", "B-field"), len(quads)),
                             dim=3, units="m", interpolator_backend="numba")
