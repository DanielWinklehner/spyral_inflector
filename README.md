# spyral_inflector

Design and simulation of cyclotron spiral inflectors: analytic electrode geometry,
BEM field solve of the built electrodes (with apertures, housing and electrostatic
quadrupoles), fringe-field optimization of the design orbit, particle tracking, and
STEP export of the result.

## Requirements

- Python 3.10+, numpy, scipy, matplotlib
- [PyPATools](https://github.com/DanielWinklehner/PyPATools) (ions, particle
  distributions, fields, pushers, trackers)
- [py_electrodes](https://github.com/DanielWinklehner/py_electrodes), branch
  `step-import-via-gmsh` or later (gmsh meshing, STEP import/export with
  transformations, collision detection)
- [bempp-cl](https://github.com/bempp/bempp-cl) with an OpenCL driver; optional
  `cupy` for the GPU GMRES, optional `dfols` for the trust-region optimizer

```bash
pip install -e .
```

## Quick start

```python
from spyral_inflector import *
from PyPATools.beam import ParticleDistribution
from PyPATools.species import IonSpecies
from PyPATools.field import Field

ion = ParticleDistribution(species=IonSpecies("H2_1+"))
ion.set_mean_energy_z_mev(0.070)

si = SpiralInflector(ion=ion, method="numerical", solver="bempp",
                     volt=12000.0,        # +V on the anode, -V on the cathode
                     gap=0.019,           # electrode gap [m]
                     tilt=31.0,           # k' tilt of the exit [deg]
                     dx=0.010,            # electrode thickness [m]
                     sigma=0.0022,        # V-shape depth [m], 0 for flat faces
                     vee_shape="parabolic",
                     ns=100,              # points along the design orbit
                     aspect_ratio=2.4,    # electrode width / gap
                     gammaAng=5.0,        # exit wedge cut [deg]
                     anglingAng=11.0)     # inner face angling [deg]

# Beam travels in +z from negative z to the median plane at z = 0; Bz is negative there.
si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -0.901}))
si.initialize()

si.set_parameter(key="h", value=0.005)                  # surface mesh size [m]
si.set_parameter(key="make_aperture", value=True)
si.set_parameter(key="aperture_params", value={"thickness": 4e-3, "radius": 50e-3,
                                               "length": 40e-3, "width": 15e-3,
                                               "top_distance": 5e-3, "bottom_distance": 10e-3,
                                               "hole_type": "rectangle", "voltage": 0.0})
si.set_parameter(key="make_housing", value=True)
si.set_parameter(key="housing_params", value={"zmin": -0.12, "zmax": 0.03, "span": True,
                                              "gap": 6e-3, "thickness": 4e-3, "voltage": 0.0,
                                              "experimental": True})
si.generate_geometry()

# Full-trajectory fringe correction: entrance/exit truncation, axial shift and the
# electrode voltage are set so that one test particle exits flat, centred and on the
# median plane. Rebuilds and re-solves the BEM model as needed.
result = si.optimize_trajectory(maxiter=15, res=0.005)

si.calculate_potential(res=0.0025, limits=((-0.1, 0.1), (-0.1, 0.1), (-0.29, 0.06)),
                       domain_decomp=(4, 4, 4))
si.calculate_efield()                                    # -> si.numerical_variables["ef_itp"]

r, v = si.fast_track(r_start=[0.0, 0.0, -0.15], v_start=[0.0, 0.0, ion.v_mean_m_per_s],
                     nsteps=2000, dt=1e-10)

for i, e in enumerate(si.numerical_variables["objects"].electrodes.values()):
    e.export("{:03d}_{}.step".format(i, e.name))        # transformations included
```

## Parameters

Constructor: `ion`, `method` (`"analytical"` or `"numerical"`), `solver` (`"bempp"`),
`volt`, `gap`, `tilt`, `dx`, `sigma`, `vee_shape`, `ns`, `aspect_ratio`, `rotation`,
`gammaAng`, `anglingAng`, `b_lim`, `debug`.

`set_parameter(key=..., value=...)`: `h`; `make_aperture` / `aperture_params`;
`make_housing` / `housing_params`; `make_cylinder` / `cylinder_params`;
`make_quadrupoles` / `quadrupole_params` (`a`, `b` hyperbola vertex radii, `radius`,
`z_starts`, `lengths`, `voltages`, `aper_rad` = hole diameter of the grounded
apertures). Quadrupole electrodes are named `D0..D3` (quad 1) and `D4..D7`
(quad 2); `D0, D1` sit on the x axis and get `+V`, `D2, D3` on the y axis get `-V`.

`optimize_trajectory(maxiter, solver="auto"|"dfols"|"broyden", res, vary_voltage=True,
exclude_quadrupoles=True, fixed={knob: value}, bounds, tol_angle, tol_offset,
tol_width, ...)`: knobs 0 entrance truncation [deg], 1 exit truncation [deg], 2 axial
shift dz [m], 3 voltage scale. Returns a dict (`db_entrance`, `db_exit`, `dz`,
`voltage`, `residual_final`, `history`, `measurements`) and leaves the object at the
best point with the full assembly solved. `optimize_fringe` is the older entrance/exit
angle iteration and is kept for compatibility.

## Conventions

- Units are SI (m, V, T); energies in MeV where the argument name says so.
- Field maps are `PyPATools.field.Field` objects (`Field.from_file`, `Field.from_arrays`,
  or `Field(dim=0, field={...})` for a uniform field). The map must cover the bore at
  negative z with the median plane at z = 0.
- The BEM solution is linear in the electrode voltages: fields for other voltage
  settings can be superposed from basis solves without re-solving.

## Examples

- `Examples/hchc60_deck.py`: the IsoDAR HCHC-60 deck (apertures, housing, quadrupole
  doublet, field map, optimizer, STEP export) with comments where it differs from
  the earlier generation script.
- `Examples/optimize_trajectory_minimal.py`: the optimizer on a uniform field, with a
  knob held fixed.
- `Examples/track_bunch_from_step.py`: reload exported STEP files with voltages,
  solve, and track a bunch with collision detection.
- `Examples/analytical_test.py`, `bempp_test.py`, `bempp_test_jm.py`: older tests.
