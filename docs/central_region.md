# How to make a central region

### Initialization

```python
from spyral_inflector import *
import numpy as np

h2p = IonSpecies("H2_1+", 0.1)
h2p.calculate_from_energy_mev(0.2 / h2p.a())

cr = CentralRegion(r_cr=[0.1, 0.3],
                   dee_voltage=70e3,
                   dee_opening_angle=42.5,
                   rf_phase=0.0,
                   ion=h2p)

cr.initialize(xi=np.array([0.05, 0.0, 0.0]),
              vi=h2p.v_m_per_s() * np.array([0.0, 1.0, 0.0]))

cr.load_bfield(bfield='/home/philip/work/C-44_AIMA_Poles_wVP_large_r_z=pm1cm.comsol',
               bf_scale=10.36,
               spatial_unit="m",
               extents=np.array([[-0.3, 0.3], [-0.3, 0.3], [-0.01, 0.01]]),
               extents_dims=["X", "Y", "Z"])
```

### Static Equilibrium Orbits

```python
initial_r, v_angle = orbit_finder(cr, energy_mev=0.4, verbose=False)

orbit_h2p = IonSpecies("H2_1+", 0.4)

r_init = np.array([initial_r, 1e-12, 0.0])
v_init = np.array([orbit_h2p.v_m_per_s() * np.sin(v_angle),
                   orbit_h2p.v_m_per_s() * np.cos(v_angle),
                   0.0])

r, v = cr_track(orbit_h2p,
                r_init=r_init,
                v_init=v_init,
                dt=5e-11,
                nterms=1,
                input_bfield=cr.analytic_parameters["bf_itp"])
```


### Making dees

```python
my_dee = AbstractDee(opening_angle=42.5,
                     r_init=0.05,
                     char_len=0.05,
                     gap=gap,
                     thickness=thickness)
my_dee.initialize()
...
...
dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]

for dee in dees:
    dee.make_transforms()
    for _ in range(3):
        dee.next_bottom_segment()
        dee.next_top_segment()

cr._abstract_dees = dees
```

### Tracking

```python
dth = 0.0
# dth = -np.pi / 4.0

simple_tracker(cr,
               r_start=np.array([0.09 * np.cos(dth),
                                 0.09 * np.sin(dth),
                                 0.0]),
               v_start=h2p.v_m_per_s() * np.array([np.cos(dth + np.pi / 2.0),
                                                   np.sin(dth + np.pi / 2.0),
                                                   0.0]))
```

### Interfacing with a spiral inflector

```python

```