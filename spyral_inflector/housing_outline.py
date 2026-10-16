"""The housing wall as a shape: a flat exit plate, one spline, and nothing else.

The housing used to be the convex hull of the electrode points offset outwards, written into
the .geo one `Point` and one `Line` per vertex -- 94 of them, with 2-5 mm edges. Those edges
pin the BEM mesh to about 4.7 mm however large `h` is, which on HCHC-60 made the housing 27 %
of the whole model (7108 of 26566 triangles) for a wall that has no detail in it at all.

This module turns that polygon into two curves, the same way `cr_geom_build.py` does it for
the central region:

  * one straight **plate** at the exit, whose outward normal the caller supplies (the design
    orbit at the exit), long enough for the opening plus padding;
  * one **spline** through points spaced evenly along the rest of the wall;
  * at each end of the plate that stops short of the end of the wall, one straight **step**
    square to the plate, from the plate's end back out to the boundary the spline follows.

The plate is a chord across the hull: over its own span it clears every electrode point by
`gap`, but outside that span it cuts across the hull's bridge, which is empty. That is what
makes it free. A plate long enough to span the whole flat side would instead have to stand
tangent to the whole hull, which on the measured geometries pushes the wall 1.6-10 mm further
out. The price of the chord is that the outline is no longer strictly convex: it picks up a
re-entrant kink of a few degrees where the plate ends.

Everything here is plain geometry in the deck frame, in metres; nothing imports gmsh.
"""

import numpy as np
from scipy.interpolate import splev, splprep
from scipy.spatial import ConvexHull, cKDTree

__all__ = ["Wall", "build_wall", "polygon_wall", "offset_convex", "plate_support",
           "clearance_report"]


def _unit(v):
    v = np.asarray(v, dtype=float)[:2]
    n = float(np.linalg.norm(v))
    if n == 0.0:
        raise ValueError("cannot normalise a zero vector")
    return v / n


def _signed_area(poly):
    return 0.5 * float(np.sum(poly[:, 0] * np.roll(poly[:, 1], -1)
                              - np.roll(poly[:, 0], -1) * poly[:, 1]))


def _ccw(poly):
    """The polygon's vertices in counter-clockwise order."""
    poly = np.asarray(poly, dtype=float)[:, :2]
    return poly if _signed_area(poly) > 0.0 else poly[::-1]


def offset_convex(poly, r, n_circle=32):
    """`poly` grown outwards by `r`: the Minkowski sum with a disc, as an `n_circle`-gon.

    This is what the housing has always done -- a circle of points around every hull vertex,
    re-hulled -- kept unchanged so the wall keeps its old clearance. 32 directions guarantee
    r * cos(pi / 32) = 99.5 % of the radius.
    """
    poly = np.asarray(poly, dtype=float)[:, :2]
    ang = 2.0 * np.pi * np.arange(n_circle) / n_circle
    circle = r * np.column_stack([np.cos(ang), np.sin(ang)])
    cloud = (poly[:, None, :] + circle[None, :, :]).reshape(-1, 2)
    return _ccw(cloud[ConvexHull(cloud).vertices])


def _spline(points, n_sample):
    """A cubic spline interpolating `points`, sampled `n_sample` times end to end."""
    points = np.asarray(points, dtype=float)
    if len(points) < 3:
        u = np.linspace(0.0, 1.0, n_sample)[:, None]
        return points[0] * (1.0 - u) + points[-1] * u
    tck, _ = splprep([points[:, 0], points[:, 1]], s=0.0, k=min(3, len(points) - 1), per=False)
    return np.column_stack(splev(np.linspace(0.0, 1.0, n_sample), tck))


def _miter_normals(nodes, closed=True):
    """Outward unit normals at each node, scaled so an offset keeps the faces parallel.

    At a corner, offsetting along the plain edge normal would round the corner off; the miter
    normal is the bisector divided by the cosine of half the turn, which puts the offset corner
    exactly where the two offset faces meet -- so the offset plate stays parallel to the plate.
    """
    nodes = np.asarray(nodes, dtype=float)
    nxt = np.roll(nodes, -1, axis=0) if closed else np.vstack([nodes[1:], nodes[-1]])
    prv = np.roll(nodes, 1, axis=0) if closed else np.vstack([nodes[:1], nodes[:-1]])
    out = []
    for e_in, e_out in zip(nodes - prv, nxt - nodes):
        a = e_in if np.linalg.norm(e_in) > 0 else e_out
        b = e_out if np.linalg.norm(e_out) > 0 else e_in
        a, b = _unit(a), _unit(b)
        n_a, n_b = np.array([a[1], -a[0]]), np.array([b[1], -b[0]])
        bis = n_a + n_b
        bis = n_b if np.linalg.norm(bis) < 1e-12 else bis / np.linalg.norm(bis)
        out.append(bis * min(1.0 / max(float(bis @ n_b), 0.2), 5.0))
    return np.array(out)


def _resample(chain, spacing, end_refine=2):
    """Points along the open polyline `chain`, evenly spaced by arc length, ends included.

    The points land ON the polyline, interpolated along its edges, so they inherit whatever the
    polyline guarantees -- here, the offset polygon's clearance.

    `end_refine` packs the first and last interval that much closer. It used to have to be 6:
    with the plate stopping short of the end of the wall, the spline had to hug the arc that
    rounds the tip of the electrodes, and a cubic cannot hold a turn like that with its nearest
    neighbour 19 mm away. Now the plate runs out to the end of that arc and swallows it, so the
    corner is a plain right angle and the refinement buys nothing: 1, 2, 3 and 6 all give it to
    within a degree. 2 keeps one interpolation point beside the corner, which is what the wall
    needs there and no more.
    """
    chain = np.asarray(chain, dtype=float)
    step = np.linalg.norm(np.diff(chain, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(step)])
    n = max(2, int(round(s[-1] / max(spacing, 1e-9))) + 1)
    want = np.linspace(0.0, s[-1], n)
    if end_refine > 1 and n > 2:
        d = want[1] - want[0]
        extra = [want[0] + j * d / end_refine for j in range(1, end_refine)] + \
                [want[-1] - j * d / end_refine for j in range(1, end_refine)]
        want = np.unique(np.concatenate([want, extra]))
    return np.column_stack([np.interp(want, s, chain[:, 0]), np.interp(want, s, chain[:, 1])])


def _dedupe(points, eps):
    """Drop points closer than `eps` to the one before them; both ends are always kept.

    Once the plate is run out to where the boundary crosses it, the boundary point at the
    plate's own t IS the plate's end, and a spline through a repeated point is not a spline --
    scipy rejects it outright.
    """
    points = np.asarray(points, dtype=float)
    kept = [0]
    for i in range(1, len(points) - 1):
        if np.linalg.norm(points[i] - points[kept[-1]]) > eps:
            kept.append(i)
    while len(kept) > 1 and np.linalg.norm(points[-1] - points[kept[-1]]) <= eps:
        kept.pop()
    return points[kept + [len(points) - 1]]


def _point_to_polyline(points, curve, chunk=2000):
    """Distance from each point to the closed polyline `curve`."""
    curve = np.asarray(curve, dtype=float)
    seg = np.roll(curve, -1, axis=0) - curve
    ll = np.einsum("ij,ij->i", seg, seg)
    ll[ll == 0.0] = 1.0
    pts = np.asarray(points, dtype=float)[:, :2]
    out = np.empty(len(pts))
    for i in range(0, len(pts), chunk):
        w = pts[i:i + chunk, None, :] - curve[None, :, :]
        u = np.clip(np.einsum("ijk,jk->ij", w, seg) / ll[None, :], 0.0, 1.0)
        out[i:i + chunk] = np.linalg.norm(w - u[:, :, None] * seg[None, :, :], axis=2).min(axis=1)
    return out


class Wall:
    """The housing wall as the .geo wants it: an inner and an outer closed loop, each a list
    of `("line" | "spline", points)` edges.

    `points` is the pair of ends for a line and the whole chain of interpolation points, ends
    included, for a spline. Consecutive edges share an end point, and the last edge closes back
    onto the first.
    """

    def __init__(self, inner_edges, outer_edges, meta=None, n_sample=400):
        self.inner_edges = list(inner_edges)
        self.outer_edges = list(outer_edges)
        self.meta = dict(meta or {})
        self._n = int(n_sample)

    @staticmethod
    def _dense(edges, n_sample):
        """The loop as a closed polyline. Straight edges are sampled too, in proportion to
        their length: leaving the exit plate as a single segment is exact for a distance test
        but hides it from anything that asks where the wall runs."""
        lengths = [float(np.linalg.norm(np.diff(np.asarray(p, dtype=float), axis=0), axis=1).sum())
                   for _, p in edges]
        total = sum(lengths) or 1.0
        out = []
        for (kind, pts), ln in zip(edges, lengths):
            pts = np.asarray(pts, dtype=float)
            if kind == "spline":
                piece = _spline(pts, n_sample)
            else:
                n = max(2, int(round(n_sample * ln / total)))
                u = np.linspace(0.0, 1.0, n)[:, None]
                piece = pts[0] * (1.0 - u) + pts[-1] * u
            out.append(piece[:-1])                              # the next edge repeats the end
        return np.vstack(out)

    @property
    def inner(self):
        return self._dense(self.inner_edges, self._n)

    @property
    def outer(self):
        return self._dense(self.outer_edges, self._n)

    def nodes(self, which="inner"):
        """The corners where the edges meet, in order."""
        edges = self.inner_edges if which == "inner" else self.outer_edges
        return np.array([np.asarray(p, dtype=float)[0] for _, p in edges])

    def loop(self, which="inner"):
        """The loop as `(points, curves)`: every defining point once, in order, and the curves
        over them as `(kind, i0, i1)` with indices into `points`, the last one wrapping to 0.

        This is what both the .geo and the hand-off to the central region are written from, so
        the two cannot describe different shapes.
        """
        edges = self.inner_edges if which == "inner" else self.outer_edges
        pts, curves = [], []
        for kind, p in edges:
            p = np.asarray(p, dtype=float)
            i0 = len(pts)
            pts.extend(p[:-1])
            curves.append((kind, i0, i0 + len(p) - 1))
        n = len(pts)
        return np.array(pts), [(k, a, b % n) for k, a, b in curves]

    def n_points(self):
        """How many `Point`s the .geo needs for both loops."""
        return sum(len(p) - 1 for _, p in self.inner_edges + self.outer_edges)

    def wall_thickness(self):
        """(min, max) distance from the inner face to the outer one -- a check on the offset.

        Away from the exit plate's two ends this holds to a couple of per cent. Within about
        2 mm of them it does not: that is the miter at a sharp corner, plus two splines through
        offset points not being exactly parallel, and it leaves the corner a little rounded or a
        little pointed. It never approaches zero, so the wall cannot pinch itself shut.
        """
        d = _point_to_polyline(self.inner, self.outer)
        return float(d.min()), float(d.max())


def polygon_wall(pts_in, pts_out):
    """The wall as it always was: every vertex a corner, every edge a straight line."""
    def edges(p):
        p = _ccw(p)
        return [("line", np.vstack([p[i], p[(i + 1) % len(p)]])) for i in range(len(p))]
    return Wall(edges(pts_in), edges(pts_out), meta={"simplified": False})


def plate_support(points, n_hat, t_hat, t_lo, t_hi, clearance):
    """Where the flat plate has to sit: `s` such that the line p . n_hat = s clears every one
    of `points` whose t falls in [t_lo, t_hi] by `clearance`.

    The window is widened by `clearance` on both sides: a point just outside it is still within
    `clearance` of the plate's end. Points further out than that are nearer the plate's end than
    its face, so the perpendicular distance used here is the conservative one.
    """
    points = np.asarray(points, dtype=float)[:, :2]
    t = points @ t_hat
    sel = (t >= t_lo - clearance) & (t <= t_hi + clearance)
    if not np.any(sel):
        raise ValueError("no electrode point lies behind the exit plate -- check n_hat "
                         "and t_centre")
    return float((points[sel] @ n_hat).max()) + clearance


def _plate_crossings(poly, n_hat, t_hat, s_plate):
    """Where the near arc of `poly` crosses the line p . n_hat = s_plate, as t values.

    The near arc is the one the plate faces; on a counter-clockwise polygon t increases along
    it, from the vertex of least t to the vertex of greatest t.
    """
    t_in, s_in, m = poly @ t_hat, poly @ n_hat, len(poly)
    i_first, i_last = int(np.argmin(t_in)), int(np.argmax(t_in))
    near = [(i_first + k) % m for k in range((i_last - i_first) % m + 1)]
    ta, sa = t_in[near], s_in[near]
    out = []
    for j in range(len(near) - 1):
        if (sa[j] - s_plate) * (sa[j + 1] - s_plate) < 0.0:
            f = (s_plate - sa[j]) / (sa[j + 1] - sa[j])
            out.append(float(ta[j] + f * (ta[j + 1] - ta[j])))
    return out


def build_wall(points_xy, gap, thickness, n_hat, t_centre, plate_half,
               knot_spacing=0.015, n_sample=400, margin=0.0, n_circle=32, hull_xy=None,
               max_bulge=0.002, extend=0.012, min_thickness=0.8):
    """Electrode points in, wall out.

    points_xy   every electrode point the housing has to clear, (N, 2), metres
    gap         clearance from the electrodes to the inner face of the wall
    thickness   wall thickness
    n_hat       outward normal of the flat exit plate, in the x-y plane
    t_centre    where the opening sits along t_hat = perp(n_hat)
    plate_half  half the plate's length: the opening's own half-length plus the padding
    knot_spacing  how far apart the spline's interpolation points sit along the wall [m]
    margin      clearance the spline must keep beyond `gap`, on top of its measured sag [m]
    hull_xy     the convex hull of `points_xy`, if it has already been computed
    max_bulge   how far the spline may stray from the polygon through its own points before
                they are packed closer together and it is built again [m]
    extend      how far the plate may be run PAST the padding to reach the end of the wall,
                where the two meet at a right angle with nothing left to round [m]
    min_thickness  the least the wall may be, as a fraction of `thickness`; the outer face is
                pushed out to hold it where the inner face turns back too tightly for a wall
                of constant thickness to exist at all

    The plate clears the electrodes by `gap` over its own span. The spline's points sit on the
    offset polygon, which carries the same `gap`, and are then pushed outwards until the curve
    BETWEEN them clears too -- a spline through points on a convex boundary cuts the corner
    between them, and a few tenths of a millimetre is a tenth of the clearance here.
    """
    points_xy = np.asarray(points_xy, dtype=float)[:, :2]
    hull = _ccw(points_xy[ConvexHull(points_xy).vertices] if hull_xy is None else hull_xy)
    n_hat = _unit(n_hat)
    t_hat = np.array([-n_hat[1], n_hat[0]])
    t_lo, t_hi = t_centre - plate_half, t_centre + plate_half

    inner = offset_convex(hull, gap, n_circle=n_circle)
    s_plate = plate_support(points_xy, n_hat, t_hat, t_lo, t_hi, gap)

    # --- where the plate should actually stop ---------------------------------------------
    # At the end of the wall the boundary turns through the arc that rounds the tip of the
    # electrodes, and a spline either hugs that arc (points crowded into the corner) or rounds it
    # off into a bulge. Running the plate out to the far end of that arc swallows it: the plate
    # covers the whole of that side, the wall leaves it along n_hat and nothing else -- one shared
    # vertex, a right angle, a sharp edge -- and there is no arc left to spend points on.
    #
    # It must not go PAST that point. Beyond it the wall has to come back in at an acute angle,
    # 45-65 deg, which a 2 mm wall mitres into a 5-8 mm spike. `extend` caps how far the plate may
    # be run: the far end of the arc is a few millimetres away on one side and tens on the other,
    # so this lengthens the plate at the end that needs it and leaves the other alone.
    # The end of the wall is the extreme of the outline in t; `extend` decides whether that is
    # near enough to run the plate to. It is a few millimetres away at the end the plate nearly
    # reaches and the whole width of the housing at the other, so this lengthens one end and
    # leaves the other exactly where the padding put it. Reaching only PART of the way -- towards
    # the extreme but not to it -- would be the worst of both, so it is all or nothing.
    for _ in range(3):
        t_arc = inner @ t_hat
        lo_end, hi_end = float(t_arc.min()), float(t_arc.max())
        # ... and it must not go PAST the end either: a plate end beyond the extreme has the wall
        # coming back at an acute angle, which is the spike. So the plate runs to the end of the
        # wall when the end is within `extend` of where the padding put it, from either side.
        new_lo = lo_end if t_lo - lo_end <= extend else t_lo
        new_hi = hi_end if hi_end - t_hi <= extend else t_hi
        if abs(new_lo - t_lo) < 1e-9 and abs(new_hi - t_hi) < 1e-9:
            break
        t_lo, t_hi = new_lo, new_hi
        s_plate = plate_support(points_xy, n_hat, t_hat, t_lo, t_hi, gap)

    e_lo = s_plate * n_hat + t_lo * t_hat
    e_hi = s_plate * n_hat + t_hi * t_hat

    # --- which vertices of the offset polygon the plate replaces --------------------------
    # The band t_lo..t_hi is an infinite strip, and it crosses the FAR side of the housing as
    # well; only the near arc counts. On a counter-clockwise polygon the outward normal is the
    # edge direction turned by -90 deg, so on the arc whose normals point along n_hat the
    # vertices run in the direction of t_hat: t increases along the near arc, monotonically,
    # from the vertex of least t to the vertex of greatest t. That makes the two ends of the
    # replaced run a plain search, and it stays right when the plate's span contains no vertex
    # at all -- which is the usual case, because the face it lands on is a single long edge.
    t_in, s_in, m = inner @ t_hat, inner @ n_hat, len(inner)
    i_first = int(np.argmin(t_in))                                   # start of the near arc
    i_last = int(np.argmax(t_in))                                    # end of it
    near = [(i_first + k) % m for k in range((i_last - i_first) % m + 1)]
    if int(np.argmax(s_in)) not in near:
        raise ValueError("the exit plate faces away from the wall -- check the sign of n_hat")

    before = near[0]
    for i in near:
        if t_in[i] <= t_lo:
            before = i
    after = near[-1]
    for i in reversed(near):
        if t_in[i] >= t_hi:
            after = i
    clamped = (t_lo < t_in[near[0]], t_hi > t_in[near[-1]])

    keep = [(after + k) % m for k in range((before - after) % m + 1)]
    if len(keep) < 4:
        raise ValueError("the exit plate would swallow the wall ({} vertices left of {})"
                         .format(len(keep), m))

    # The wall picks the boundary up exactly where the plate stops, not at the next vertex: the
    # plate's ends land in the middle of a long edge, so starting the chain at `after` left a
    # 29 mm step where the target spacing was 15, and the spline overshot 3 mm right at the kink.
    #
    # The boundary point at the plate's own t is also what makes the corner a right angle. The
    # step from the plate's end to it is along n_hat and nothing else, so the wall leaves the
    # plate at normal incidence -- which is the corner a 2 mm wall mitres cleanly (2.83 mm) and
    # the one that does not leave a spike. So those two points are KEPT, and only the boundary
    # between them is resampled; letting the resampling step over them left one end at 45-65 deg
    # with an 8 mm spike on the outer face.
    arc = inner[near]
    t_arc = arc @ t_hat

    def _at_t(value):
        j = int(np.clip(np.searchsorted(t_arc, value), 1, len(arc) - 1))
        f = (value - t_arc[j - 1]) / max(t_arc[j] - t_arc[j - 1], 1e-15)
        return arc[j - 1] + np.clip(f, 0.0, 1.0) * (arc[j] - arc[j - 1])

    B_hi, B_lo = _at_t(t_hi), _at_t(t_lo)

    # --- the step from a plate end back out to the boundary ----------------------------------
    # Where the plate is not run out to the end of the wall its end sits a few millimetres INSIDE
    # the offset boundary (it is a chord), and the wall has to step back out -- along n_hat, since
    # the boundary point is taken at the plate's own t, so the step is square to the plate. That
    # step used to be the first interval of the spline, with the right angle at one end and the
    # boundary's own corner at the other, and a cubic cannot turn like that: on the e071 doublet
    # it overshot by 3.4 mm, halving the spacing four times (149 knots) changed nothing, and the
    # outer face wrapped the hook in a 6 mm knob while pinching to 0.9 mm against the plate's end.
    # The step is a straight edge of its own now, like the plate, and the spline starts on the
    # boundary -- a corner may sit where two edges meet, never inside a spline. A step shorter
    # than half the wall is not worth an edge (it would only pin the mesh) and is dropped, the
    # spline then starting at the plate's end, as it always has at the end run out to the wall.
    min_step = 0.5 * thickness
    step_hi = float(np.linalg.norm(B_hi - e_hi)) >= min_step
    step_lo = float(np.linalg.norm(B_lo - e_lo)) >= min_step
    if not step_hi:
        B_hi = e_hi
    if not step_lo:
        B_lo = e_lo
    # CCW order runs ... before -> [the plate] -> after ..., and t increases along it, so the
    # plate is traversed from its t_lo end to its t_hi end and the chain from B_hi round to B_lo.
    chain = np.vstack([B_hi, inner[keep], B_lo])

    # --- the spline's interpolation points, evenly spaced ALONG the boundary ---------------
    # Not every k-th vertex: the offset polygon's vertices are anything but evenly spaced --
    # a long hull edge contributes two of them and a hull corner a 32-gon arc of them -- so
    # every 5th vertex came out 0.6 mm apart in one place and 46 mm apart in the next, and a
    # cubic through points spaced like that overshoots badly (6.2 mm of housing that should
    # not be there). Even arc length removes the overshoot at its source.
    #
    # They are then pushed out until the curve clears the electrodes everywhere. The test is the
    # real one -- distance to the electrode points, not "outside the offset polygon": the plate is
    # deliberately inside that polygon, and so is the spline beside it, because both cut across
    # the hull's empty bridge. A spline through points ON the offset boundary otherwise cuts the
    # corner between them by a few tenths of a millimetre, which is a tenth of the clearance, so
    # the sag is measured and the interpolation points moved out by it. The two plate ends are
    # held, so the plate stays exactly where it was put.
    tree = cKDTree(points_xy)
    spacing = float(knot_spacing)
    for _attempt in range(4):
        # The boundary only: `chain` begins and ends at the points the steps lead to (or at the
        # plate's own ends where there is no step), and those two are held.
        knots = _dedupe(_resample(chain, spacing), 0.05 * spacing)
        free = np.ones(len(knots), dtype=bool)
        free[0] = free[-1] = False
        nrm = _miter_normals(knots, closed=False)
        push, worst = np.zeros(len(knots)), -np.inf
        for _ in range(12):
            moved = knots + push[:, None] * nrm
            curve = _spline(moved, 4 * n_sample)
            deficit = (gap + margin) - tree.query(curve)[0]
            worst = float(deficit.max())
            if worst <= 0.0:
                break
            # blame each offending sample on the nearest point that is allowed to move; blaming
            # a plate end instead would leave the deficit with nothing to fix it
            movable = np.nonzero(free)[0]
            nearest = movable[np.argmin(
                np.linalg.norm(curve[:, None, :] - moved[None, movable, :], axis=2), axis=1)]
            need = np.zeros(len(knots))
            for i in movable:
                hit = deficit[nearest == i]
                if hit.size:
                    need[i] = max(0.0, float(hit.max()))
            if not np.any(need > 0.0):
                break
            push[free] += need[free] * 1.05
        knots_out = knots + push[:, None] * nrm
        # How far does the spline stray from the polygon through its own points? That is the
        # interpolation error alone -- measuring against the offset polygon instead would mix it
        # up with the outward push, which is deliberate. Clearance is safe either way (the push
        # sees to that), but a bulge is housing that should not be there, and it means the points
        # are too far apart for the shape they have to follow. The polygon is the open chain,
        # forth and back, so no closing segment short-cuts the measure; the plate and the steps
        # are straight edges and have no bulge to measure.
        bulge = float(_point_to_polyline(_spline(knots_out, 4 * n_sample),
                                         np.vstack([knots_out, knots_out[-2:0:-1]])).max())
        if bulge <= max_bulge:
            break
        spacing *= 0.5

    inner_edges = ([("line", np.vstack([e_hi, knots_out[0]]))] if step_hi else []) + \
        [("spline", knots_out)] + \
        ([("line", np.vstack([knots_out[-1], e_lo]))] if step_lo else []) + \
        [("line", np.vstack([e_lo, e_hi]))]

    # --- the outer loop: the same shape, `thickness` further out --------------------------
    # The interior points are taken off the TRUE offset of the inner curve, not off the miter
    # of the polygon through its points: the miter is right where two straight faces meet, but
    # between them the inner curve bends, and a spline through mitered points cut the corner and
    # left the wall 0.87 mm thick where 2.00 was asked. The two ENDS keep the miter, because
    # there the wall really is two straight faces meeting -- the plate and the wall leaving it --
    # and the miter is what keeps the outer plate parallel to the inner one and the loop closed.
    _dense_in = _spline(knots_out, 8 * n_sample)
    _tan = np.gradient(_dense_in, axis=0)
    _tan /= np.linalg.norm(_tan, axis=1)[:, None]
    _out_nrm = np.column_stack([_tan[:, 1], -_tan[:, 0]])          # outward, for a CCW loop
    _at = np.argmin(np.linalg.norm(knots_out[:, None, :] - _dense_in[None, :, :], axis=2), axis=1)
    out_pts = knots_out + thickness * _out_nrm[_at]

    # Every corner is mitred, so the outer faces stay parallel to the inner ones and meet
    # exactly: the plate and the step at their right angle, the step (or the plate, where there
    # is no step) and the spline at whatever angle the boundary leaves at. The spline's own
    # tangent sets the latter, not the chord to its next knot.
    def _mitre(prev, node, nxt):
        return node + thickness * _miter_normals(np.vstack([prev, node, nxt]), closed=False)[1]

    B_hi, B_lo = knots_out[0], knots_out[-1]
    o_hi = _mitre(e_hi if step_hi else e_lo, B_hi, B_hi + _tan[0])
    o_lo = _mitre(B_lo - _tan[-1], B_lo, e_lo if step_lo else e_hi)
    out_pts[0], out_pts[-1] = o_hi, o_lo
    o_e_hi = _mitre(e_lo, e_hi, B_hi) if step_hi else o_hi
    o_e_lo = _mitre(B_lo, e_lo, e_hi) if step_lo else o_lo

    # Where the inner face turns back on itself more tightly than the wall is thick -- it does,
    # just past the plate end that is not run out to the end of the wall -- a wall of constant
    # thickness does not exist: the outer face would cross itself. Left alone it pinched to 8 um.
    # So the outer face is pushed out until it keeps `min_thickness` of the wall everywhere,
    # which fills that one corner out a little. The inner face is untouched, so the clearance and
    # the plate are exactly as computed; only the back of the shell moves.
    def _seg(p0, p1):
        return p0 + np.linspace(0.0, 1.0, 64)[:, None] * (p1 - p0)
    _loop_in = np.vstack([_dense_in] + ([_seg(B_lo, e_lo)] if step_lo else []) + [_seg(e_lo, e_hi)]
                         + ([_seg(e_hi, B_hi)] if step_hi else []))
    _free = np.ones(len(out_pts), dtype=bool)
    _free[0] = _free[-1] = False
    _nrm_o = _miter_normals(out_pts, closed=True)
    # The push is capped at one wall thickness and the best iterate is kept: right beside a plate
    # end the only point that could fix the pinch is the plate's own outer corner, which is held
    # so the outer plate stays parallel, and an uncapped loop just pushes its neighbours until the
    # wall is 7 mm thick somewhere else without ever closing the gap.
    _push_o, thin = np.zeros(len(out_pts)), np.inf
    _best, _movable = np.zeros(len(out_pts)), np.nonzero(_free)[0]
    for _ in range(12):
        _moved_o = out_pts + _push_o[:, None] * _nrm_o
        _curve_o = _spline(_moved_o, 4 * n_sample)
        _short = float((min_thickness * thickness
                        - _point_to_polyline(_curve_o, _loop_in)).max())
        if _short < thin:
            thin, _best = _short, _push_o.copy()
        if _short <= 0.0 or not len(_movable):
            break
        _near_o = _movable[np.argmin(
            np.linalg.norm(_curve_o[:, None, :] - _moved_o[None, _movable, :], axis=2), axis=1)]
        _deficit = min_thickness * thickness - _point_to_polyline(_curve_o, _loop_in)
        _need_o = np.zeros(len(out_pts))
        for i in _movable:
            hit = _deficit[_near_o == i]
            if hit.size:
                _need_o[i] = max(0.0, float(hit.max()))
        if not np.any(_need_o > 0.0):
            break
        _push_o[_free] = np.minimum(_push_o[_free] + _need_o[_free] * 1.05, thickness)
    out_pts = out_pts + _best[:, None] * _nrm_o
    outer_edges = ([("line", np.vstack([o_e_hi, out_pts[0]]))] if step_hi else []) + \
        [("spline", out_pts)] + \
        ([("line", np.vstack([out_pts[-1], o_e_lo]))] if step_lo else []) + \
        [("line", np.vstack([o_e_lo, o_e_hi]))]

    wall = Wall(inner_edges, outer_edges, n_sample=n_sample, meta={
        "simplified": True, "n_hull": len(hull), "n_offset": len(inner),
        "n_knots": len(knots_out), "plate_s": float(s_plate),
        "knot_spacing": spacing, "bulge": bulge, "thin": thin,
        "plate_span": (float(t_lo), float(t_hi)), "plate_ends": (e_hi, e_lo),
        "steps": (float(np.linalg.norm(B_hi - e_hi)) if step_hi else 0.0,
                  float(np.linalg.norm(B_lo - e_lo)) if step_lo else 0.0),
        "n_hat": n_hat, "t_hat": t_hat, "worst": worst, "push_max": float(push.max()),
        "clamped": clamped,
        "offset_polygon": inner, "gap": float(gap), "thickness": float(thickness)})
    return wall


def clearance_report(wall, points, gap):
    """Smallest distance from `points` to the wall's inner face, against the `gap` asked for.

    `worst` below zero means the wall comes closer to an electrode than `gap` somewhere.
    """
    d = _point_to_polyline(points, wall.inner)
    i = int(np.argmin(d))
    return {"min": float(d[i]), "worst": float(d[i] - gap),
            "at": np.asarray(points, dtype=float)[i][:2],
            "n_closer": int(np.sum(d < gap - 1e-9))}
