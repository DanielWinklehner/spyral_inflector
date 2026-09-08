# HCHC-60 spiral-inflector study: the scripts behind the results

These are the driver scripts of the 2026 HCHC-60 injection-line study, copied from the
working deck. They sit on top of `spyral_inflector.tracking` (the package does the
building, solving, tracking, metrics and plots); the scripts organize scans, retunes and
reports. Every script has a docstring with its command line.

## The deck folder

The scripts expect a *deck* folder with

```
<deck>/Fields/      B-field map pickle (PyPATools Field), e.g. HCHC-60_CentralBField_z-40to5cm_1mm.pickle
<deck>/Particles/   RFQ output: TraceWin .dst or the 10-column text file (x(mm) x'(mrad) y(mm) y'(mrad)
                    z(mm) z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss); .dst is read directly
<deck>/Geometry/    STEP exports of the geometries (one file per electrode), written by the scripts
<deck>/Results/     everything the scripts write (json, npz, png, md reports, openPMD hand-off files)
<deck>/Scripts/     these scripts (or set SI_DECK to the deck folder and run them from anywhere)
```

By default a script takes the deck to be the parent of the folder it lives in; the
environment variable `SI_DECK` overrides that.

## Environment

`environment-windows.yml` (exported from the working machine) or `environment-ubuntu.yml`
(derived, untested) at the repository root, then the three project packages editable:
`pip install -e PyPATools -e py_electrodes -e spyral_inflector`. A CUDA GPU is assumed
(warp collision kernels, cupy in the space-charge solver); without one everything falls
back to the CPU and is several times slower.

## One geometry, one run

* `RunFinalOnce.py`: direct BEM solve of a STEP folder at given quad voltages/rotations and
  spiral voltage, full core beam with and without space charge, exit metrics, figures,
  openPMD hand-off files, `report.md`. Same as `python -m spyral_inflector.tracking.run_final`.
* `../hchc60_pg5L_final.py`: the 20-line version that also builds the geometry from its knobs.

## Building a geometry and tuning its quads

* `GeometryPoint.py`: one geometry point: build from knobs, design-particle optimization of
  spiral voltage and axial shift, STEP export, three basis fields with the quads rotated,
  quad-voltage retune by superposition of the basis fields, a bunch, `summary.json`.
* `BunchScan2.py`: transmission on a grid of quad voltages, rotations, spiral-voltage scales
  and beam angles by superposing basis fields (spiral, q1, q1 skew, q2, q2 skew).
* `scan_best.py`, `scan_fit.py`, `scan_complete.py`, `merge_scans.py`, `fine_phi_setup.py`,
  `combine_moves.py`: helpers for picking, fitting, guarding and merging scan results.

## Scans and optimizations

* `FocusScan.py`: electrode-shape grid (V-depth `sigma`, `anglingAng`, `gammaAng`) with a
  GeometryPoint per shape, refinement of the best at finer fields, space-charge confirmation.
  Keep the gamma wedge at or below about 11 deg: 14 deg gave a mesh the BEM solver crawls on.
* `SCRetune.py`: quad voltages and rotations retuned with space charge (superposed fields,
  8 mA, 2 mm cells), then a full-beam confirmation. Solves its own unrotated basis fields:
  basis solves of a finished chain carry that chain's rotation and must not be reused.
* `FringeOptBunch.py`: exit truncation and spiral voltage optimized on the bunch with the
  centroid kept level; voltage and axial shift are re-levelled together (they are coupled).
  Its truncation choice is still made before the levelling; see the docstring.
* `RunNight2.ps1` (template) and `RunPG5L_run.ps1` (filled): the overnight chain that
  produced the pg5L geometry (beam angle, quad rotations, voltage grids, spiral scale, direct
  solves, resolution check, dz correction, full-beam runs, report). Windows PowerShell.
* `make_final_report.py`: the chain's report.

## Wrappers kept for the scripts above

`TrackFromStep.py`, `BunchTrack.py`, `BunchTrackSC.py`, `exit_metrics.py`,
`plot_geometry_trajectories.py` call the package; `track_inflector.py` is a shim exposing
the package under the old module name.

## Inputs

* `read_dst.py`: statistics of a TraceWin .dst, comparison with an older text file, export
  to the text format (with the bunch-core selection: |phase| <= 180 deg, energy >= half the
  median). The package reads .dst files directly with the same selection.
* `convert_bfield.py`, `flip_bfield_z.py`, `bfield_axis.py`, `validate_bfield.py`: field-map
  conversion to a PyPATools Field pickle, axis flip, on-axis profile and checks.

## Conventions

Deck frame: origin on the cyclotron axis in the median plane, z along the axis, the beam
travels in +z and arrives at z = 0; plots draw -z upwards. Transmission counts particles
through the housing exit opening. Quad voltages: D0,D1 = +q1, D2,D3 = -q1, D4,D5 = +q2,
D6,D7 = -q2. Space-charge results at 2 mm cells read about 1.9 points, at 1.5 mm 1.3 and at
1 mm 0.8 points above the converged value.
