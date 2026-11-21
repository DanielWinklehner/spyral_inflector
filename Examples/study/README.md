# Study scripts: scans, optimizations and retunes on top of the package

Driver scripts of an injection-line study. They sit on top of `spyral_inflector.tracking` (the package does the
building, solving, tracking, metrics and plots); the scripts organize scans, retunes and reports. Every script
has a docstring with its command line. The levelling chain itself (rotation fixed point, space-charge levelling,
Baseline-frame export, reports) is in the package: `spyral_inflector.workflows`, driven by a design config, see the
top-level README.

The scripts carry no machine: the beam file (`--particles`), the B map (`--bfield`), the gap ray and half gap
(`--gap-azimuth`, `--half-gap`) are arguments, and the species, RF and start plane of the injected bunch come from
the machine settings (`spyral_inflector.machine`: `SI_SPECIES`, `SI_RF_MHZ`, `SI_Z_START`, `SI_INJECT_CORE`) in
the environment.

## The deck folder

The scripts expect a *deck* folder with

```
<deck>/Fields/      B-field map pickles (PyPATools Field)
<deck>/Particles/   injected beam: TraceWin .dst or the 10-column text file (x(mm) x'(mrad) y(mm) y'(mrad)
                    z(mm) z'(mrad) Phase(deg) Time(s) Energy(MeV) Loss); .dst is read directly
<deck>/Geometry/    STEP exports of the geometries (one file per electrode), written by the scripts
<deck>/Results/     everything the scripts write (json, npz, png, md reports, openPMD hand-off files)
<deck>/Scripts/     these scripts (or set SI_DECK to the deck folder and run them from anywhere)
```

By default a script takes the deck to be the parent of the folder it lives in; the environment variable `SI_DECK`
overrides that, and `SI_RESULTS` puts the results elsewhere.

## Environment

`environment-windows.yml` or `environment-ubuntu.yml` (derived, untested) at the repository root, then the three
project packages editable: `pip install -e PyPATools -e py_electrodes -e spyral_inflector`. A CUDA GPU is assumed
(warp collision kernels, cupy in the space-charge solver); without one everything falls back to the CPU and is
several times slower.

## One geometry, one run

* `RunFinalOnce.py`: direct BEM solve of a STEP folder at given quad voltages/rotations and spiral voltage, full
  core beam with and without space charge, exit metrics, figures, openPMD hand-off files, `report.md`. Same as
  `python -m spyral_inflector.tracking.run_final`.
* `RunFinalPush.py`: field intake, geometry, exit-plane rotation, basis solves, quad retune, bunch runs and
  Baseline-frame exports of one design in eleven resumable phases (`../final_push_minimal.py` is the same on
  generated inputs).

## Building a geometry and tuning its quads

* `GeometryPoint.py`: one geometry point: build from knobs, design-particle optimization of spiral voltage and
  axial shift, STEP export, basis fields with the quads rotated, quad-voltage retune by superposition, a bunch,
  `summary.json`.
* `QuadBases.py`: the basis fields of an existing run (spiral, each quad and its skew).
* `BunchScan2.py`: transmission on a grid of quad voltages, rotations, spiral-voltage scales and beam angles by
  superposing basis fields.
* `TripletScan.py`, `QuadOptN.py`, `SCGridN.py`: N-quad geometries, their quad optimization (DFO-LS on superposed
  bases) and space-charge grids.
* `scan_best.py`, `scan_fit.py`, `scan_complete.py`, `merge_scans.py`, `fine_phi_setup.py`, `combine_moves.py`:
  helpers for picking, fitting, guarding and merging scan results.

## Scans and optimizations

* `FullOpt.py`, `FullOptRadius.py`: shape optimization of a levelled design (knobs, longitudinal profiles, exit
  height and direction targets) with the quads retuned in the loop.
* `FocusScan.py`, `VFocusScan.py`, `VFocusRefine.py`, `VFocusRank.py`: electrode-shape grids with a geometry point
  per shape, refinement of the best, space-charge confirmation. Keep the gamma wedge at or below about 11 deg: 14
  deg gave a mesh the BEM solver crawls on.
* `SCRetune.py`, `SCRetuneV.py`: quad voltages and rotations retuned with space charge (superposed fields), then a
  full-beam confirmation. Basis solves of a finished chain carry that chain's rotation and must not be reused.
* `FringeOptBunch.py`: exit truncation and spiral voltage optimized on the bunch with the centroid kept level.
* `RotationScan2.py`, `DoubletRotScan.py`: rotation scans; `InjectCompare.py`: spatial bunch vs timed injection.

## Wrappers kept for the scripts above

`TrackFromStep.py`, `BunchTrack.py`, `BunchTrackSC.py`, `exit_metrics.py`, `plot_geometry_trajectories.py` call
the package; `track_inflector.py` is a shim exposing the package under the old module name (set `ti.BFIELD`,
`ti.PARTICLES`, `ti.STEP_DIR` before use).

## Inputs

* `read_dst.py`: statistics of a TraceWin .dst, comparison with a text file, export to the text format (with the
  bunch-core selection: |phase| <= 180 deg, energy >= half the median). The package reads .dst files directly.
* `convert_bfield.py`, `flip_bfield_z.py`, `validate_bfield.py`: field-map conversion to a PyPATools Field pickle,
  the 180-degree turn about x, and checks. `SI_PREVIOUS_BFIELD` names an earlier map that `RunFinalPush.py` and
  `PlotIntermediate.py` compare the axis with.

## Conventions

Deck frame: origin on the cyclotron axis in the median plane, z along the axis, the beam travels in +z and arrives
at z = 0; plots draw -z upwards. Transmission counts particles through the housing exit opening. Quad voltages:
D0,D1 = +q1, D2,D3 = -q1, D4,D5 = +q2, D6,D7 = -q2. Space-charge results at 2 mm cells read about 1.9 points, at
1.5 mm 1.3 and at 1 mm 0.8 points above the converged value (one design's resolution study).
