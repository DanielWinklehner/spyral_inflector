"""openPMD hand-off files for the central region: the transmitted particles at the
crossing of a plane (position, momentum, time, RF phase) or a 6-D snapshot at one time.
Format: PyPATools' documents/openpmd_plane_crossing_handoff.md."""
import numpy as np

from .deck import CLIGHT
from .hooks import plane_axes

FRAME = ("spyral_inflector deck frame: origin on the cyclotron axis in the median plane, z along the axis "
         "(beam travels +z, Bz < 0 at the centre), x/y of the deck; SI units, momenta in eV/c. "
         "Machine frame (+z up, beam from the top, orbits counter-clockwise from above, +x = theta 0 on a hill): "
         "the mirror image through the median plane, z -> -z and pz -> -pz, x/y/azimuths unchanged")

# openPMD-beamphysics only knows a few species names; map the PyPATools names onto them so
# ParticleGroup(h5=...) can compute energies (the exact PyPATools mass is stored alongside)
PMD_SPECIES_NAMES = {"H2_1+": "H2+", "H_1+": "proton", "proton": "proton", "H_1-": "H-", "electron": "electron"}


def _pmd_species(ion):
    name = getattr(ion, "name", "H2_1+")
    species = {"name": PMD_SPECIES_NAMES.get(name, name), "mass_mev": float(ion.mass_mev),
               "a": getattr(ion, "a", None), "z": getattr(ion, "z", None),
               "charge_state": getattr(ion, "q", getattr(ion, "charge_state", 1))}
    return {k: val for k, val in species.items() if val is not None}


def _betagamma(v):
    gamma = 1.0 / np.sqrt(1.0 - np.sum(v * v, axis=1) / CLIGHT ** 2)
    return v * gamma[:, None] / CLIGHT


def save_handoff_openpmd(path, plane, sel, ion, rf_hz, current_ma, trj, vdes, handoff_distance_m,
                         phase_reference="mean", extra_meta=None, mode="plane_crossing"):
    """Write the particles that crossed `plane` (a PlaneCrossing) to an openPMD-beamphysics
    HDF5 file: lab-frame position and momentum at the interpolated crossing, the crossing
    time, the RF phase relative to the reference crossing time, in-plane coordinates u, v
    and the plane definition as attributes. Returns a dict with the numbers written."""
    from PyPATools.particles_src.particle_io import save_openpmd
    sel = np.asarray(sel, dtype=bool) & plane.crossed & np.isfinite(plane.time)
    n_tracked = int(len(sel))
    r, v, t = plane.state[sel, :3], plane.state[sel, 3:], plane.time[sel]
    t_ref = float(np.median(t)) if phase_reference == "median" else float(np.mean(t))
    phase = np.angle(np.exp(1j * 2.0 * np.pi * rf_hz * (t - t_ref)))     # (-pi, pi]
    u_ax, v_ax, w_ax = plane_axes(plane.normal)
    d = r - plane.point
    n = int(sel.sum())
    q_injected = current_ma * 1e-3 / rf_hz if rf_hz > 0 else 0.0     # charge per RF bucket of the injected beam
    q_bunch = q_injected * n / max(n_tracked, 1)                       # the fraction that reached the plane
    L, T = (1, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0)
    meta = dict(bunch_charge=q_bunch, bunch_freq=rf_hz, time=t, status=np.ones(n, dtype=int),
                spec_version=1, injected_bunch_charge_c=q_injected, transmitted_current_ma=float(current_ma * n / max(n_tracked, 1)),
                pypatools_species_name=getattr(ion, "name", "H2_1+"),
                extra_records={"phase": phase, "u": d @ u_ax, "v": d @ v_ax, "t_rel": t - t_ref},
                extra_units={"phase": (1.0, (0,) * 7), "u": (1.0, L), "v": (1.0, L), "t_rel": (1.0, T)},
                mode=mode, plane_origin_m=plane.point, plane_normal=w_ax, plane_axis_u=u_ax, plane_axis_v=v_ax,
                handoff_distance_m=float(handoff_distance_m), phase_reference=phase_reference,
                handoff_definition=("plane through the design particle's continued orbit at handoff_distance_m of path "
                                    "past its exit point, perpendicular to the design velocity there (plane_normal)"),
                phase_reference_time_s=t_ref, rf_frequency_hz=float(rf_hz), beam_current_ma=float(current_ma),
                design_exit_point_m=np.asarray(trj[-1], dtype=float), design_exit_velocity_mps=np.asarray(vdes[-1], dtype=float),
                frame=FRAME,
                momentum_note="px,py,pz are lab-frame momenta in eV/c (beamphysics convention); position in m at the crossing",
                n_particles=n)
    if extra_meta:
        meta.update(extra_meta)
    save_openpmd(path, r, _betagamma(v), _pmd_species(ion), **meta)
    return {"n": n, "t_ref_s": t_ref, "phase_rms_deg": float(np.degrees(phase.std())), "t_rms_ns": float(1e9 * t.std()),
            "u_rms_mm": float(1e3 * (d @ u_ax).std()), "v_rms_mm": float(1e3 * (d @ v_ax).std())}


def save_snapshot_openpmd(path, snap, sel, ion, rf_hz, current_ma, extra_meta=None):
    """Write a 6-D lab-frame snapshot (r, v, active, t) of the selected particles to openPMD."""
    from PyPATools.particles_src.particle_io import save_openpmd
    r, v, active, t = snap
    m = np.asarray(sel, dtype=bool) & active
    n_tracked = int(len(m))
    n = int(m.sum())
    q_injected = current_ma * 1e-3 / rf_hz if rf_hz > 0 else 0.0
    q_bunch = q_injected * n / max(n_tracked, 1)
    meta = dict(bunch_charge=q_bunch, bunch_freq=rf_hz, time=np.full(n, t), status=np.ones(n, dtype=int),
                spec_version=1, injected_bunch_charge_c=q_injected, transmitted_current_ma=float(current_ma * n / max(n_tracked, 1)),
                pypatools_species_name=getattr(ion, "name", "H2_1+"),
                mode="lab6d_snapshot", snapshot_time_s=float(t), rf_frequency_hz=float(rf_hz), beam_current_ma=float(current_ma),
                frame=FRAME, n_particles=n)
    if extra_meta:
        meta.update(extra_meta)
    save_openpmd(path, r[m], _betagamma(v[m]), _pmd_species(ion), **meta)
    return {"n": n, "t_s": float(t)}
