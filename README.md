# spyral_inflector
A python module to help with the design of cyclotron spiral inflectors.

## Installation

**Prerequisites:**
This modules has several required prerequisite modules to be installed. Some of the common packages that you may already have are:
- numpy
- scipy
- matplotlib

In addition to those, you will need to have:

- [PyPATools](https://github.com/DanielWinklehner/PyPATools)
- [py_electrodes](https://github.com/DanielWinklehner/py_electrodes)

As of May 2019, two solvers can be used to calculate fields. You only need one installed to be able to utilize the features in the spyral_inflector modules.
- [bempp-cl](https://github.com/bempp/bempp-cl)
- [FEniCS](https://fenicsproject.org/)

_Note: At present, the FEniCS solver does not fully work. Fixing this is high on
our priority list!_

### Installation using [Anaconda3](https://www.anaconda.com/)
[PyPATools](https://github.com/DanielWinklehner/PyPATools) comes with an 
_environment.yml_ file that can be used to directly create a conda environment
in both Windows and Ubuntu 18 (other Linux distributions may well work, but
haven't been tested yet). 

The PyPATools conda environment is meant to include everything needed to run
programs in the 
[Python Particle Accelerator Tools](https://github.com/users/DanielWinklehner/projects/1) 
project and includes installation of py_electrodes and bempp-cl through pip+git.

_Note: Unfortunately, bempp-cl requires OpenCL drivers 
installed and the Windows Subsystem for Linux (WSL) does not support that kind
of hardware access (yet)._

### Installing the module
```bash
git clone https://github.com/DanielWinklehner/spyral_inflector.git
cd spyral_inflector
```

The module can be installed using the setup.py file:
```bash
python setup.py install
```
or through pip (or pip3):
```bash
pip install .
```

or directly from git:

```bash
pip install git+https://github.com/DanielWinklehner/spyral_inflector.git
```

The git installation works without cloning the repository first. Examples can be
directly downloaded from the
[examples](https://github.com/DanielWinklehner/spyral_inflector/tree/master/Examples) 
folder.

_Note: If you are not using a separate environment (conda or env), we recommend using
the ``--user`` option with pip to install modules locally._

## Generating an analytical model

Required Parameters:
```
"ion": Ion species (object from dans_pymodules)
"volt": Voltage difference applied to the electrodes
"gap": Distance between the electrodes (assuming they are parallel)
"tilt": Tilt angle of the exit electrodes
"dx": Thickness of the electrodes
"ns": Number of steps in the analytical solution
```
Optional Parameters:
```
"b_lim": Range of values for the b angle.
"rotation": Axial rotation of the spiral inflector
"debug": Default is false
```
The analytical model of the spiral inflector is the simplest one the code can create. It only depends on the listed parameters above and takes very little time to create. The geometry is created by calculating the central trajectory through the spiral inflector using the equations of motion that describe the system.

**Simplest possible example:**

First start by importing the spyral_inflector module, which will include imports of dans_pymodules, numpy, matplotlib, and bempp.
```python
from spyral_inflector import *
```

The next step is to define the IonSpecies from dans_pymodules, which will contain all of the information about the ion that's being injected into the spiral inflector.
```python
h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.07 / h2p.a())
```

The spiral inflector object itself is only instantiated with a few parameters, and others (for BEMPP, optimization, etc.) are set through another method (see the numerical inflector section).
```python
si = SpiralInflector(ion=h2p,
                     method="analytical",
                     solver="bempp",    # Can be either "bempp" or "fenics"
                     volt=12000,        # Electrode voltages
                     gap=18e-3,         # Electrode gap distance
                     tilt=0,            # Tilt angle
                     dx=10e-3,          # Electrode thickness
                     sigma=0,           # v-shape parameter
                     ns=60,             # Number of trajectory points
                     debug=False)
```

The magnetic field of the cyclotron can be defined in a few different ways, but for this situation we are only concerned with a constant magnetic field pointing in the negative z direction. This is done with the field object from dans_pymodules.
```python
si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))
```

The spiral inflector can now be initialized. The purpose of initializing is to perform some of the basic calculations for variables used in other methods (i.e. rotations) and to check if required parameters are missing. After initializing, the analytical geometry of the spiral inflector can be generated. In the geometry generation, the method for calculating the central trajectory is called internally and used to calculate the edge geometry of the electrodes.
```python
si.initialize()
si.generate_geometry()
```
Lastly, all this information is used to create a macro that can be used in AutoDesk Inventor that creates a full solid model of the spiral inflector generated by the code.
**(Note: soon this will be replaced with exports from OpenCASCADE)**
```python
export_electrode_geometry(si, fname='electrode_macro.ivb') 
```

## Generating a numerical model
Generating the numerically solved spiral inflector is a significantly more computationally intensive process than creating the analytical model. The idea is to create a meshed model of the spiral inflector using the information about the analytical model and perform an electrostatic finite element analysis to calculate the electric fields. These electric fields are then used the calculate the true trajectory of an ion being injected straight into the inflector. However, fringing electric fields at the entrance and exit of the inflector will cause this ion to deviate from a "central trajectory". To mediate this, segments of the entrance and exit can be "cut off". The optimization routine repeats this process until the ion leaves the inflector on the mid/median plane of the cyclotron.

**Example:**

The process for creating the numerical model starts similar to that of the analytical model. The only main difference is that the method is now "numerical".
```python
from spyral_inflector import *

h2p = IonSpecies("H2_1+", 0.035)
h2p.calculate_from_energy_mev(0.07 / h2p.a())

si = SpiralInflector(ion=h2p,
                     method="numerical",
                     solver="bempp",
                     volt=12000,
                     gap=18e-3,
                     tilt=27.0,
                     dx=10e-3,
                     sigma=1.5E-3,
                     ns=60,
                     debug=False)

si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

si.initialize()

si.generate_geometry()
```

Now, the parameters that will be used in the meshing and finite element analysis need to be set. The first parameter is "h" which is a characteristic length used by the meshing program (gmsh) to create the full meshed model. The smaller the h, the more fine of a model you will generate. Warning: it is very easy to set this parameter too small and the finite element analysis will never finish.

```python
si.set_parameter(key="h", value=0.005)
```

To create an even more realistic model of the spiral inflector system as a whole, you can add apertures at the entrance and exit of the device to reduce fringe field effects. Additionally, a cylindrical boundary may be created to simulate a housing for the electrodes. At the end, the apertures can be exported as an Inventor macro in the same way the inflector electrodes are.

```python
si.set_parameter(key="make_aperture", value=True)
si.set_parameter(key="aperture_params", value={"thickness": 4e-3,
                                               "radius": 50e-3,
                                               "length": 45e-3,
                                               "width": 18e-3,
                                               "top_distance": 5e-3,
                                               "bottom_distance": 10e-3,
                                               "voltage": 0.0})
si.set_parameter(key="make_cylinder", value=True)
si.set_parameter(key="cylinder_params", value={"radius": 120e-3,
                                               "zmin": -150e-3,
                                               "zmax": 80e-3,
                                               "voltage": 0.0})
```

Now, the numerical model can be meshed and the optimization routine can begin. maxiter sets the maximum number of iterations allowed if the exit angle does fall between the tolerance bounds. The optimization will automatically regenerate the meshed model after it finishes.

```python
generate_meshed_model(si)
optimize_fringe(si, maxiter=5, tol=0.02, res=0.005)
```

With the optimized geometry ready, the next thing to do is get the electric field and central trajectory.

```python
calculate_efield(si, res=0.002,
                     limits=((-0.08, 0.08), (-0.08, 0.08), (-0.12, 0.05)),
                     domain_decomp=(7, 7, 7))

track(si, r_start=np.array([0.0, 0.0, -0.15]),
          v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
          nsteps=15000,
          dt=1e-11)
```

The electric field may be handy if you need to import it into another program (in that case, don't use the domain decomposition). Now, the process is nearly complete, and all that needs to be done is exporting the electrode and aperture geometries.

```python
draw_geometry(si, show=True, filename='auto')
export_electrode_geometry(si, fname='electrode_macro.ivb')
export_aperture_geometry(si, fname='aperture_macro.ivb')
```