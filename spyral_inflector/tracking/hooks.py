"""Terminators, recorders and interactions for PyPATools' Tracker: electrode collisions,
the inflector exit plane, the hand-off plane crossing, trajectory/envelope/snapshot
recorders and the space-charge re-solve."""
import time

import numpy as np

from PyPATools.pusher import Pusher
from PyPATools.trackers import Tracker, Terminator, Recorder


class ElectrodeCollision(Terminator):
    """Kill particles whose step crosses an electrode surface, and record where."""

    def __init__(self, assembly, skip=None):
        self.assembly = assembly
        self.skip = skip            # an ExitPlane: particles past it are not collision-checked ...
        # ... except against these electrodes (index in the assembly, electrode) for this many
        # steps after the crossing: the housing exit opening sits ~10-14 mm past the electrode
        # end, so a beam that clears the exit plane can still clip it on the way out.
        self.post_exit_electrodes = []
        self.post_exit_steps = 0
        self.post_exit_hit = None
        self.post_exit_electrode = None
        self.post_exit_point = None
        self.names = [e.name for e in assembly.electrodes.values()]
        self.hit_step = None
        self.hit_electrode = None
        self.hit_point = None

    def update(self, step, r_prev, v_prev, r, v, active, t):
        if self.hit_step is None:
            n = len(r)
            self.hit_step = np.full(n, -1, dtype=int)
            self.hit_electrode = np.full(n, -1, dtype=int)
            self.hit_point = np.full((n, 3), np.nan)

        idx = np.where(active)[0]
        if self.skip is not None and self.skip.crossed is not None and idx.size:
            active = self._post_exit_check(step, r_prev, r, active)
            idx = idx[~self.skip.crossed[idx]]

        if idx.size == 0:
            return active

        data = self.assembly.segment_intersects_surface(r_prev[idx], r[idx])
        hit = data["hit_mask"]
        if hit.any():
            gone = idx[hit]
            self.hit_step[gone] = step
            self.hit_electrode[gone] = data["electrode_ids"][hit]
            self.hit_point[gone] = data["hit_points"][hit]
            active = active.copy()
            active[gone] = False
        return active

    def _post_exit_check(self, step, r_prev, r, active):
        """Collision check of coasting particles against the post-exit electrodes only."""
        if self.post_exit_steps <= 0 or not self.post_exit_electrodes:
            return active
        if self.post_exit_hit is None:
            self.post_exit_hit = np.zeros(len(r), dtype=bool)
            self.post_exit_electrode = np.full(len(r), -1, dtype=int)
            self.post_exit_point = np.full((len(r), 3), np.nan)
        sel = active & self.skip.crossed & ~self.post_exit_hit & (step - self.skip.step <= self.post_exit_steps)
        idx = np.where(sel)[0]
        if idx.size == 0:
            return active
        for i, e in self.post_exit_electrodes:
            hit, pts, _fr = e.segment_intersects_surface(r_prev[idx], r[idx])
            if hit.any():
                gone = idx[hit]
                self.post_exit_hit[gone] = True
                self.post_exit_electrode[gone] = i
                self.post_exit_point[gone] = pts[hit]
                active = active.copy()
                active[gone] = False
                idx = idx[~hit]
                if idx.size == 0:
                    break
        return active

    def losses_by_electrode(self, index_map):
        out = {}
        for i in np.unique(self.hit_electrode[self.hit_electrode >= 0]):
            uuid = index_map.get(int(i))
            label = self.assembly.electrodes[uuid].name if uuid else "electrode_{}".format(i)
            out[label] = int((self.hit_electrode == i).sum())
        return out


class ExitPlane(Terminator):
    """Retire particles once they cross the design exit plane, recording their state.

    Crossing the plane outward near the design exit point means the particle has left
    the inflector; what happens afterwards is the central region's business. A pure
    "radius >= r_exit" test would also fire for particles that undershoot the exit and
    reach that radius later on their first turn.

    coast_steps > 0 keeps a crossed particle alive for that many more steps (for
    trajectory plots and the asymptotic state); pair it with
    ElectrodeCollision(skip=<this object>) so coasting particles are not counted as
    losses when they meet the housing on the first turn. asym_steps after the crossing
    the state is stored again (asym_state): outside the electric fringe, the height and
    vertical angle the central region actually sees.
    """

    def __init__(self, r_exit, point=None, normal=None, window=0.05, coast_steps=0, asym_steps=120):
        self.r_exit = r_exit
        self.window = window   # max in-plane distance from the design exit point
        self.point = None if point is None else np.asarray(point, dtype=float)
        self.normal = None
        if normal is not None:
            n = np.asarray(normal, dtype=float)
            self.normal = n / np.linalg.norm(n)
        self.coast_steps = int(coast_steps)
        self.asym_steps = min(int(asym_steps), self.coast_steps) if self.coast_steps > 0 else 0
        self.asym_state = None
        self.coast_until = None
        self.crossed = None
        self.state = None
        self.step = None
        self.offset = None

    def update(self, step, r_prev, v_prev, r, v, active, t):
        if self.crossed is None:
            self.crossed = np.zeros(len(r), dtype=bool)
            self.state = np.full((len(r), 6), np.nan)
            self.step = np.full(len(r), -1, dtype=int)
            self.offset = np.full(len(r), np.nan)
            self.coast_until = np.full(len(r), -1, dtype=int)
            self.asym_state = np.full((len(r), 6), np.nan)

        if self.asym_steps > 0:
            due = active & self.crossed & (step - self.step == self.asym_steps)
            if due.any():
                self.asym_state[due, :3] = r[due]
                self.asym_state[due, 3:] = v[due]

        if self.normal is None:
            reached = np.linalg.norm(r[:, :2], axis=1) >= self.r_exit
        else:
            d = r - self.point
            ahead = d @ self.normal >= 0.0
            outward = (v @ self.normal) > 0.0
            in_plane = np.linalg.norm(d - np.outer(d @ self.normal, self.normal), axis=1)
            reached = ahead & outward & (in_plane <= self.window)

        new = active & ~self.crossed & reached
        if new.any():
            self.crossed[new] = True
            self.state[new, :3] = r[new]
            self.state[new, 3:] = v[new]
            self.step[new] = step
            d = r[new] - self.point
            self.offset[new] = np.linalg.norm(d - np.outer(d @ self.normal, self.normal), axis=1)
            active = active.copy()
            if self.coast_steps > 0:
                self.coast_until[new] = step + self.coast_steps
            else:
                active[new] = False

        if self.coast_steps > 0:
            expired = active & self.crossed & (self.coast_until >= 0) & (step >= self.coast_until)
            if expired.any():
                active = active.copy()
                active[expired] = False
        return active


class PlaneCrossing(Terminator):
    """Record each particle's first crossing of a plane, interpolated within the step
    (crossing time good to a small fraction of dt). Particles are not terminated; with
    gate = an ExitPlane only particles it has marked as crossed are watched."""

    def __init__(self, point, normal, gate=None, dt=None):
        self.point = np.asarray(point, dtype=float)
        n = np.asarray(normal, dtype=float)
        self.normal = n / np.linalg.norm(n)
        self.gate = gate
        self._dt = None if dt is None else float(dt)
        self.crossed = None
        self.state = None      # (N, 6) r, v at the crossing
        self.time = None       # (N,) crossing time [s]
        self.step = None

    def update(self, step, r_prev, v_prev, r, v, active, t):
        if self.crossed is None:
            self.crossed = np.zeros(len(r), dtype=bool)
            self.state = np.full((len(r), 6), np.nan)
            self.time = np.full(len(r), np.nan)
            self.step = np.full(len(r), -1, dtype=int)
        s_prev = (r_prev - self.point) @ self.normal
        s_now = (r - self.point) @ self.normal
        new = active & ~self.crossed & (s_prev < 0.0) & (s_now >= 0.0)
        if self.gate is not None and getattr(self.gate, "crossed", None) is not None:
            new &= self.gate.crossed
        if new.any():
            frac = np.clip(s_prev[new] / (s_prev[new] - s_now[new]), 0.0, 1.0)[:, None]
            self.state[new, :3] = r_prev[new] + frac * (r[new] - r_prev[new])
            self.state[new, 3:] = v_prev[new] + frac * (v[new] - v_prev[new])
            self.time[new] = t + frac[:, 0] * self._dt if self._dt is not None else np.nan
            self.step[new] = step
            self.crossed[new] = True
        return active


class TrajectoryRecorder(Recorder):
    """Positions of an evenly spaced subsample of the bunch every few steps, NaN once a
    particle is gone (for 3-D plots of the geometry with trajectories)."""

    def __init__(self, n_particles, n_record=400, every=4):
        self.idx = np.unique(np.linspace(0, n_particles - 1, min(max(n_record, 1), n_particles)).astype(int))
        self.every = every
        self.rows = []

    def record(self, step, r_prev, v_prev, r, v, active, t):
        if step % self.every == 0:
            x = np.array(r[self.idx], dtype=np.float32)
            x[~active[self.idx]] = np.nan
            self.rows.append(x)
        return None

    @property
    def data(self):
        if not self.rows:
            return np.zeros((0, len(self.idx), 3), dtype=np.float32)
        return np.array(self.rows, dtype=np.float32)


class SnapshotRecorder(Recorder):
    """Full lab-frame state of every live particle at given steps (6-D snapshots)."""

    def __init__(self, steps):
        self.steps = set(int(k) for k in steps)
        self.snapshots = {}

    def record(self, step, r_prev, v_prev, r, v, active, t):
        if step in self.steps:
            self.snapshots[step] = (r.copy(), v.copy(), active.copy(), float(t))
        return None


class Envelope(Recorder):
    """Centroid and rms size of the live particles every few steps."""

    def __init__(self, every=4):
        self.every = every
        self.rows = []

    def record(self, step, r_prev, v_prev, r, v, active, t):
        if step % self.every == 0 and active.any():
            ra = r[active]
            self.rows.append([step, t, int(active.sum()), ra[:, 0].mean(), ra[:, 1].mean(), ra[:, 2].mean(),
                              ra[:, 0].std(), ra[:, 1].std(), ra[:, 2].std()])
        return None


class SpaceCharge:
    """Re-solve the space-charge field of the live particles every n steps and swap it
    into slot 1 of the CompositeField [vacuum, space charge]."""

    def __init__(self, solver, composite, charges, every, min_particles=0):
        self.solver = solver
        self.composite = composite
        self.charges = charges
        self.every = every
        # below this many depositing particles the field is not re-solved (the last
        # stragglers alone carry no space charge worth a Poisson solve)
        self.min_particles = int(min_particles)
        self.n_skipped = 0
        self.n_solves = 0
        self.solve_time = 0.0
        self.log = []
        self.exclude = None         # an ExitPlane: particles past it no longer deposit charge

    def solve_now(self, r, active, step, verbose=False):
        t0 = time.time()
        if self.exclude is not None and self.exclude.crossed is not None:
            active = active & ~self.exclude.crossed
        if not np.any(active):
            return
        if active.sum() < self.min_particles:
            self.n_skipped += 1
            return
        phi, e_sc = self.solver.solve(r[active], self.charges[active])
        self.composite.fields[1] = e_sc
        wall = time.time() - t0
        self.solve_time += wall
        self.n_solves += 1
        e_on_beam = np.linalg.norm(e_sc(r[active]), axis=1)
        row = [step, int(active.sum()), float(phi.min()), float(phi.max()),
               float(e_on_beam.max()), float(np.median(e_on_beam)), wall]
        self.log.append(row)
        if verbose or self.n_solves % 20 == 1:
            print("    SC solve {:4d} @ step {:5d}: {:5d} live, phi {:+.1f}..{:+.1f} V, "
                  "|E_sc| on beam max {:.2e} / median {:.2e} V/m, {:.2f} s".format(
                      self.n_solves, step, row[1], row[2], row[3], row[4], row[5], wall), flush=True)

    def apply(self, step, r_prev, v_prev, r, v, active, t, dt):
        # called after the push of `step`; the new field serves the next `every` steps
        if (step + 1) % self.every == 0 and np.any(active):
            self.solve_now(r, active, step + 1)
        return r, v, active


class DelayedInjection:
    """Inject particles at their own arrival times.

    The RFQ's unaccelerated tail arrives at the entrance plane up to many RF periods
    after the core, so it cannot be placed as a spatial bunch: its particles start with
    alive False, parked at their injection state, and are switched on at their injection
    step. In the static inflector fields the delay only matters through space charge,
    which is exactly where it has to be right. Pair with
    Tracker.run(stop_on_all_lost=False), so the run outlives the core."""

    def __init__(self, inject_step, r_inject, v_inject):
        self.k = np.asarray(inject_step, dtype=int)
        self.r_inject = np.asarray(r_inject, dtype=float)
        self.v_inject = np.asarray(v_inject, dtype=float)
        self.pending = self.k > 0
        self.injected_step = np.where(self.pending, -1, 0)

    @property
    def alive0(self):
        """The initial alive mask: everything that is not waiting to be injected."""
        return ~self.pending

    def apply(self, step, r_prev, v_prev, r, v, active, t, dt):
        # placed during step k-1 (after that step's push), so the first push of step k moves
        # the particle from exactly its injection time k * dt
        due = self.pending & (self.k <= step + 1)
        if due.any():
            r[due] = self.r_inject[due]
            v[due] = self.v_inject[due]
            r_prev[due] = self.r_inject[due]      # no phantom segment for the collision check
            v_prev[due] = self.v_inject[due]
            active = active.copy()
            active[due] = True
            self.pending[due] = False
            self.injected_step[due] = step
        return r, v, active


class _PD(object):
    """The minimal particle-distribution interface Tracker.run expects."""

    def __init__(self, r, v):
        self.x_vec, self.v_vec = r, v
        self.alive = np.ones(len(r), dtype=bool)

    def set_p_from_v_vec(self, v):
        return None


def plane_axes(normal):
    """In-plane unit vectors (u, v) of a plane with unit normal w: u horizontal (z x w),
    v = w x u (close to +z, i.e. pointing along the beam's original direction of
    travel, which is DOWN in the deck's -z-up pictures); (u, v, w) right-handed."""
    w = np.asarray(normal, dtype=float)
    w = w / np.linalg.norm(w)
    u = np.cross([0.0, 0.0, 1.0], w)
    if np.linalg.norm(u) < 1e-9:
        u = np.array([1.0, 0.0, 0.0])
    u = u / np.linalg.norm(u)
    return u, np.cross(w, u), w


def continue_design(ion, efield, bfield, r_exit, v_exit, dt, path_m):
    """Track the design particle on from its exit state through the same fields and return
    (point, unit velocity) where its path length past the exit reaches path_m: the hand-off
    plane goes through that point, perpendicular to the design orbit there (the orbit bends
    by ~35 deg over 30 mm in the cyclotron field, so the exit direction is the wrong normal)."""
    r0 = np.asarray(r_exit, dtype=float).reshape(1, 3)
    v0 = np.asarray(v_exit, dtype=float).reshape(1, 3)
    if path_m <= 0.0:                                   # hand-off at the electrode exit itself
        return r0[0].copy(), v0[0] / np.linalg.norm(v0[0])
    n_steps = int(1.5 * path_m / (float(np.linalg.norm(v0)) * dt)) + 20
    snap = SnapshotRecorder(range(0, n_steps))
    Tracker(Pusher(ion, algorithm="rk4_rel"), efield, bfield, terminators=[], recorders=[snap]).run(
        _PD(r0.copy(), v0.copy()), dt, n_steps, show_progress=False, sync_back=False)
    ks = sorted(snap.snapshots)
    R = np.vstack([r0] + [snap.snapshots[k][0][:1] for k in ks])
    V = np.vstack([v0] + [snap.snapshots[k][1][:1] for k in ks])
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(R, axis=0), axis=1))])
    if s[-1] < path_m:
        raise RuntimeError("design particle only travelled {:.1f} mm in {} steps".format(1e3 * s[-1], n_steps))
    i = int(np.searchsorted(s, path_m))
    f = (path_m - s[i - 1]) / (s[i] - s[i - 1])
    point = R[i - 1] + f * (R[i] - R[i - 1])
    vel = V[i - 1] + f * (V[i] - V[i - 1])
    return point, vel / np.linalg.norm(vel)
