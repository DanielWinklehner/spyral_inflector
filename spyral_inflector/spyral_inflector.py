from .field_solving import *
from .geometry import *
from .optimization import *
from .plotting import *
from .trajectories import *

from dans_pymodules import *
import numpy as np

__author__ = "Daniela Campo, Philip Weigel, Daniel Winklehner"
__doc__ = """
This module calculates the geometry of a spiral inflector for given initial parameters:
ion: mass, charge, energy
spiral inflector: voltage, external b-field
geometric parameters: thickness, gap, v-shape angle

This module is based on a MATLAB code written by Daniela Campo, which in turn is based on the
paper: Toprek NIM A 440 (2000).
"""

# Initialize some global constants
amu = const.value("atomic mass constant energy equivalent in MeV")
echarge = const.value("elementary charge")
clight = const.value("speed of light in vacuum")
colors = MyColors()


def orthogonal_proj(zfront, zback):
    a = (zfront + zback) / (zfront - zback)
    b = -2 * (zfront * zback) / (zfront - zback)
    # -0.0001 added for numerical stability as suggested in:
    # http://stackoverflow.com/questions/23840756
    return np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, a, b],
                     [0, 0, -0.0001, zback]])


class SpiralInflector(object):
    def __init__(self,
                 method="analytical",
                 solver="bempp",
                 debug=False,
                 outp_folder="",
                 **kwargs):

        # --- Program Variables -------------------------------------------------------------------------------------- #
        self._method = method  # Either analytical or numerical
        self._solver = solver  # Either 'bempp' or 'fenics'

        # assert self._method in ["analytical", "numerical"], \
        #     "Spiral inflector method must be 'analytical' or 'numerical'!"

        self._debug = debug
        self._outp_folder = outp_folder
        self._temp_folder = None
        self._initialized = False

        # --- Semi-fixed Spiral Inflector Parameters ----------------------------------------------------------------- #
        # Note: After changing one of these parameters, the spiral inflector has to be recalculated from scratch
        self._params_analytic = {"ion": None,  # ion species to design the spiral inflector for (from dans_pymodules)
                                 "bf_itp": None,  # type: Field
                                 "volt": None,  # The voltage applied to the electrodes (V)
                                 # Geometry:
                                 "gap": None,  # The gap between the two electrodes (m)
                                 "tilt": None,  # The tilt angle of the exit of the spiral inflector (degrees)
                                 "dx": None,  # The thickness of the spiral inflector electrodes (m)
                                 "aspect_ratio": None,  # Aspect ratio between electrode width and thickness
                                 "sigma": None,  # The "v-shape" of the electrodes (m)
                                 "ns": None,  # Resolution of the analytical solution (steps along beam trajectory s)
                                 "b_lim": np.deg2rad(np.array([0.0, 90.0])),  # Limits of the curvature
                                 "rotation": 0.0,  # Clockwise rotation of the spiral inflector, assuming that at
                                 # 0.0 deg the entrance E-field points in the x-dir
                                 }

        for key in self._params_analytic.keys():
            if key in kwargs.keys():
                self._params_analytic[key] = kwargs[key]

        # --- Parameters used in the spiral inflector calculation ---------------------------------------------------- #
        self._variables_analytic = {"bf_design": 0.0,
                                    # An avg/max? field for the analytical formula calculated from bf_itp
                                    "ef_design": 0.0,  # Design E-Field (2*V/dx) (V/m)
                                    "r_cyc": 0.0,  # Cyclotron radius of the design particle
                                    "b": None,  # Array holding the possible values of curvature
                                    "trj_design": None,  # type: np.ndarray
                                    "v_design": None,  # type: np.ndarray
                                    "height": 0.0,  # Height of the inflector (extension in z-direction) (m)
                                    "kp": 0.0,  # Tilt parameter
                                    "k": 0.0,
                                    "d": None,  # Inner gap size vs. deflection angle
                                    "orbit_center": (None, None)  # x,y position of orbit center
                                    }

        # --- Parameters used by the BEM++ potential and field calculation ------------------------------------------- #
        self._params_numerical = {"h": None,  # the desired mesh spacing for BEM++ mesh generation
                                  "make_aperture": False,  # Make apertures at the exit and entrance
                                  "gmres_tol": 1E-5,  # Tolerance used to calculate the BEMPP solution
                                  "aperture_params": {"thickness": None,
                                                      "radius": None,
                                                      "length": None,
                                                      "width": None,
                                                      "top_distance": None,
                                                      "bottom_distance": None,
                                                      "hole_type": None,
                                                      "voltage": 0.0},
                                  "make_cylinder": False,  # Cylindrical boundary
                                  "cylinder_params": {"radius": None,
                                                      "zmin": None,
                                                      "zmax": None,
                                                      "voltage": 0.0},
                                  "make_housing": False,  # Convex hull housing
                                  "housing_params": {"zmin": None,
                                                     "zmax": None,
                                                     "span": None,
                                                     "gap": None,
                                                     "thickness": None,
                                                     "voltage": 0.0,
                                                     "experimental": None}
                                  }

        for key in self._params_numerical.keys():
            if key in kwargs.keys():
                self._params_numerical[key] = kwargs[key]

        self._variables_numerical = {"objects": {},  # A dictionary of objects (apertures, electrodes)
                                     "limits": None,  # BEM++ Potential/E-Field Limits
                                     "full mesh": None,  # Full mesh of the geometry
                                     "i": None,  # A running 3-index for mesh generation
                                     "f_space": None,  # BEM++ Function Space
                                     "operator": None,  # BEM++ Operator
                                     "grid_fun": None,  # BEM++ Grid Function
                                     "solution": None,  # The BEM++ solution
                                     "n_fun_coeff": None,  # BEM++ Neumann coefficients from the solution
                                     "ef_phi": None,  # Object containing electric potential
                                     "ef_itp": None,  # Interpolator object for E-Field
                                     }

        # --- Additional parameters used for particle tracking ------------------------------------------------------- #
        self._params_track = {"dt": 1e-10,
                              "nsteps": 10000}
        self._variables_track = {"trj_tracker": None,
                                 "shift": None,  # type: np.ndarray
                                 }

        # --- Additional parameters used for optimization ------------------------------------------------------------ #
        self._variables_optimization = {"x_rot": 0
                                        }
        # --- Experimental parameters -------------------------------------------------------------------------------- #
        self._params_exp = {"y_opt": False
                            }
        for key in self._params_exp.keys():
            if key in kwargs.keys():
                self._params_exp[key] = kwargs[key]

    def __str__(self):

        printstr = "Spiral Inflector: {}\n".format(self.__repr__())
        printstr += " * Current Parameters:\n"
        printstr += " ** Fixed (need complete recalculation):\n"

        for key, item in self._params_analytic.items():
            printstr += "     * {}: {}\n".format(key, item)

        printstr += " ** Derived Parameters:\n"

        for key, item in self._variables_analytic.items():
            printstr += "     * {}: {}\n".format(key, item)

        return printstr

    def generate_design_trajectory(self):
        r = None
        if self._method == "analytical":
            r = self.generate_analytical_trajectory()
        elif self._method == "numerical":
            r = self.generate_numerical_trajectory()
        return r

    def generate_geometry(self):
        if self._method == "analytical":
            return self.generate_analytical_geometry()
        elif self._method == "numerical":
            return self.generate_numerical_geometry()

        return 1

    def _calculate_bf_max(self):

        # TODO: think about all of this!

        assert isinstance(self._params_analytic["bf_itp"],
                          Field), "It seems no B-field interpolator has been loaded yet."

        _x = np.zeros(10000)
        _y = np.zeros(10000)
        _z = np.linspace(-1.0, 1.0, 10000)
        _points = np.vstack([_x, _y, _z]).T

        _bf = self._params_analytic["bf_itp"]  # type: Field
        _, _, _bz = _bf(_points)

        self._variables_analytic["bf_design"] = np.round(np.sign(np.mean(_bz)) * np.max(abs(_bz)), 2)

    def get_parameter(self, key):

        if key in self._params_analytic.keys():

            return self._params_analytic[key]

        elif key in self._variables_analytic.keys():

            return self._variables_analytic[key]

        else:
            return None

    @property
    def debug(self):
        return self._debug

    @property
    def initialized(self):
        return self._initialized

    @property
    def method(self):
        return self._method

    @property
    def solver(self):
        return self._solver

    @property
    def analytic_parameters(self):
        return self._params_analytic

    @property
    def analytic_variables(self):
        return self._variables_analytic

    @property
    def bempp_parameters(self):
        return self._params_numerical

    @property
    def bempp_variables(self):
        return self._variables_numerical

    @property
    def numerical_parameters(self):
        return self._params_numerical

    @property
    def numerical_variables(self):
        return self._variables_numerical

    @property
    def track_parameters(self):
        return self._params_track

    @property
    def track_variables(self):
        return self._variables_track

    @analytic_parameters.setter
    def analytic_parameters(self, analytic_parameters):
        self._params_analytic = analytic_parameters

    @analytic_variables.setter
    def analytic_variables(self, analytic_variables):
        self._variables_analytic = analytic_variables

    @bempp_parameters.setter
    def bempp_parameters(self, bempp_parameters):
        self._params_numerical = bempp_parameters

    @bempp_variables.setter
    def bempp_variables(self, bempp_variables):
        self._variables_numerical = bempp_variables

    @numerical_parameters.setter
    def numerical_parameters(self, numerical_parameters):
        self._params_numerical = numerical_parameters

    @numerical_variables.setter
    def numerical_variables(self, numerical_variables):
        self._variables_numerical = numerical_variables

    @track_parameters.setter
    def track_parameters(self, track_parameters):
        self._params_track = track_parameters

    @track_variables.setter
    def track_variables(self, track_variables):
        self._variables_track = track_variables

    def initialize(self):

        abort_flag = False
        for key, item in self._params_analytic.items():
            if item is None:
                print("Item {} is not set in analytic parameters!".format(key))
                abort_flag = True
        if abort_flag:
            exit(1)

        print("Initialization...", end="")

        _ion = self._params_analytic["ion"]  # type: IonSpecies
        volt = self._params_analytic["volt"]  # type: float
        gap = self._params_analytic["gap"]  # type: float

        self._calculate_bf_max()

        self._variables_analytic["b"] = np.linspace(self._params_analytic["b_lim"][0],
                                                    self._params_analytic["b_lim"][1],
                                                    self._params_analytic["ns"])
        self._variables_analytic["r_cyc"] = abs(_ion.b_rho() / self._variables_analytic["bf_design"])  # (m)
        self._variables_analytic["ef_design"] = 2.0 * volt / gap  # (V/m)
        self._variables_analytic["height"] = abs(2.0 * _ion.total_kinetic_energy_ev() / (_ion.q() * 2.0 * volt / gap))

        # Reset the design trajectory to None
        self._variables_analytic["trj_design"] = None  # type: np.ndarray

        # Reset the analytical geometry
        self._variables_analytic["geo"] = None  # type: np.ndarray

        # Create global rotation matrix
        _rot = np.deg2rad(self._params_analytic["rotation"])
        self._variables_analytic["rot"] = np.array([[np.cos(_rot), -np.sin(_rot), 0.0],
                                                    [np.sin(_rot), np.cos(_rot), 0.0],
                                                    [0.0, 0.0, 1.0]])

        self._initialized = True

        print("Done!")

        if self._debug:
            print(self)

    def load_bfield(self, bfield=None, bf_scale=1.0, spatial_unit="cm"):

        bf_load_from_file = False

        if isinstance(bfield, int):
            bfield = float(bfield)

        if isinstance(bfield, float):

            self._params_analytic["bf_itp"] = Field(label="Cyclotron B-Field",
                                                    dim=0,
                                                    field={"x": 0.0, "y": 0.0, "z": bfield * bf_scale},
                                                    units=spatial_unit)

        elif isinstance(bfield, str):

            if os.path.exists(bfield):

                print("Loading B-Field from file {}".format(os.path.split(bfield)[1]))
                bf_load_from_file = True

            else:

                print("Couldn't find B-Field file at {} opening file dialog.".format(bfield))
                bfield = None

        elif isinstance(bfield, Field):

            self._params_analytic["bf_itp"] = bfield

        if bfield is None:

            _fd = FileDialog()
            bfield = _fd.get_filename("open")

            if bfield is not None:

                bf_load_from_file = True

            else:

                exit(0)

        if bf_load_from_file:
            self._params_analytic["bf_itp"] = Field(label="Cyclotron B-Field",
                                                    scaling=bf_scale,
                                                    units=spatial_unit)

            self._params_analytic["bf_itp"].load_field_from_file(filename=bfield)

            print("Successfully loaded B-Field from file")

        if self._initialized:
            self.initialize()

    def _set_blim(self, b_min=None, b_max=None):

        reinit = False

        if b_min is not None:
            assert b_min >= 0.0, "b_min has to be >= 0, found {}".format(b_min)
            self._params_analytic["b_lim"][0] = np.deg2rad(b_min)
            reinit = True

        if b_max is not None:
            assert b_max <= 90.0, "b_max has to be <= 90, found {}".format(b_max)
            self._params_analytic["b_lim"][1] = np.deg2rad(b_max)
            reinit = True

        if reinit:
            self.initialize()

        return 0

    def set_parameter(self, group=None, key=None, value=None):

        if key in self._params_analytic.keys():

            self._params_analytic[key] = value

            # If fixed parameter is changed, recalculation is necessary!
            if self._initialized:
                self.initialize()

        elif key in self._params_numerical.keys():

            self._params_numerical[key] = value

        elif group in self._params_numerical.keys():

            if key in self._params_numerical[group].keys():
                self._params_numerical[group][key] = value

        elif key in self._params_exp:

            self._params_exp[key] = value

    def save(self, fname):
        import pickle

        saved_numerical_vars = {}
        saved_numerical_vars["full mesh"] = self._variables_numerical["full mesh"]
        saved_numerical_vars["n_fun_coeff"] = self._variables_numerical["n_fun_coeff"]
        saved_numerical_vars["limits"] = self._variables_numerical["limits"]

        save_obj = {"method": self._method,
                    "solver": self._solver,
                    "params_analytic": self._params_analytic,
                    "variables_analytic": self._variables_analytic,
                    "params_numerical": self._params_numerical,
                    "variables_numerical": saved_numerical_vars,
                    "params_track": self._params_track,
                    "variables_track": self._variables_track,
                    "variables_optimization": self._variables_optimization}

        with open(fname, 'wb') as outfile:
            pickle.dump(save_obj, outfile)

        return 0

    def load(self, fname):
        import pickle

        with open(fname, 'rb') as infile:
            save_obj = pickle.load(infile)

            self._method = save_obj["method"]
            self._solver = save_obj["solver"]

            for key, val in save_obj["params_analytic"].items():
                self._params_analytic[key] = val

            for key, val in save_obj["variables_analytic"].items():
                self._variables_analytic[key] = val

            for key, val in save_obj["params_numerical"].items():
                self._params_numerical[key] = val

            for key, val in save_obj["variables_numerical"].items():
                self._variables_numerical[key] = val

            for key, val in save_obj["params_track"].items():
                self._params_track[key] = val

            for key, val in save_obj["variables_track"].items():
                self._variables_track[key] = val

            for key, val in save_obj["variables_optimization"].items():
                self._variables_optimization[key] = val

        self.initialize()

        return 0

    # Function wrappers below
    def calculate_efield(self):
        if self._solver == 'bempp':
            return calculate_efield_bempp(self)
        elif self._solver == 'fenics':
            return calculate_efield_fenics(self)

    def calculate_potential(self, **kwargs):
        return calculate_potential(self, **kwargs)

    def solve(self):
        if self._solver == 'bempp':
            return self.solve_bempp()
        elif self._solver == 'fenics':
            return self.solve_fenics()

    def solve_bempp(self):
        return solve_bempp(self)

    def solve_fenics(self):
        return solve_fenics(self)

    # def generate_aperture_geometry(self, *args):
    #     return generate_aperture_geometry(self, *args)
    #
    # def generate_cylinder_geometry(self):
    #     return generate_cylinder_geometry(self)
    #
    # def generate_spiral_electrode_geometry(self, *args):
    #     return generate_spiral_electrode_geometry(self, *args)

    def generate_vacuum_space(self):
        return generate_vacuum_space(self)

    def generate_solid_assembly(self, **kwargs):
        return generate_solid_assembly(self, **kwargs)

    def generate_meshed_model(self, **kwargs):
        return generate_meshed_model(self, **kwargs)

    def save_geo_files(self):
        return save_geo_files(self)

    def aperture_geometry_macro(self, **kwargs):
        return aperture_geometry_macro(self, **kwargs)

    def electrode_geometry_macro(self, **kwargs):
        return electrode_geometry_macro(self, **kwargs)

    def optimize_fringe(self, **kwargs):
        return optimize_fringe(self, **kwargs)

    def draw_geometry(self, **kwargs):
        return draw_geometry(self, **kwargs)

    def draw_mesh(self):
        return draw_mesh(self)

    # def plot_potential(self, **kwargs):
    #     return plot_potential(self, kwargs)

    def plot_bfield(self, **kwargs):
        return plot_bfield(self, kwargs)

    def plot_trajectories(self, **kwargs):
        return plot_trajectories(self, **kwargs)

    def generate_analytical_trajectory(self):
        return generate_analytical_trajectory(self)

    def generate_numerical_trajectory(self):
        return generate_numerical_trajectory(self)

    def generate_analytical_geometry(self):
        return generate_analytical_geometry(self)

    def generate_numerical_geometry(self):
        return generate_numerical_geometry(self)

    def track(self, **kwargs):
        return track(self, **kwargs)

    def fast_track(self, **kwargs):
        return fast_track(self, **kwargs)

    def fast_track_with_termination(self, **kwargs):
        return fast_track_with_termination(self, **kwargs)

    def calculate_orbit_center(self, k=None, kp=None, height=None):
        if k is None:
            k = self.analytic_variables["k"]
        if kp is None:
            kp = self.analytic_variables["kp"]
        if height is None:
            height = self.analytic_variables["height"]
        return calculate_orbit_center(k, kp, height)


if __name__ == "__main__":
    h2p = IonSpecies("H2_1+", 0.035)
    h2p.calculate_from_energy_mev(0.07 / h2p.a())

    si = SpiralInflector(ion=h2p,
                         method="numerical",
                         volt=12000,
                         gap=18e-3,
                         tilt=27.0,
                         dx=10e-3,
                         aspect_ratio=2.5,
                         sigma=1.5E-3,
                         ns=60,
                         debug=False)

    si.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

    si.initialize()

    si.generate_geometry()
    draw_geometry(si, freq=10, show=True)

    si.set_parameter(key="h", value=0.01)  # Mesh characteristic length
    si.set_parameter(key="make_aperture", value=True)
    si.set_parameter(key="aperture_params", value={"thickness": 4e-3,
                                                   "radius": 50e-3,
                                                   "length": 45e-3,
                                                   "width": 18e-3,
                                                   "top_distance": 5e-3,
                                                   "bottom_distance": 10e-3,
                                                   "voltage": 0.0})
    si.set_parameter(key="make_cylinder", value=False)
    si.set_parameter(key="cylinder_params", value={"radius": 120e-3,
                                                   "zmin": -150e-3,
                                                   "zmax": 80e-3,
                                                   "voltage": 0.0})

    generate_meshed_model(si)

    ts = time.time()
    optimize_fringe(si, initial_guess=(2.8683, -4.9849), maxiter=5, tol=0.02, res=0.005)
    print("Optimizing took {:.4f} s".format(time.time() - ts))

    ts = time.time()
    print("Calculating electric field...")
    si.calculate_potential(res=0.002,
                           limits=((-0.08, 0.08), (-0.08, 0.08), (-0.12, 0.05)),
                           domain_decomp=(7, 7, 7))
    si.calculate_efield()
    print("Calculating field took {:.4f} s".format(time.time() - ts))

    with open('timing.txt', 'a') as outfile:
        outfile.write("Generating electric field took {:.4f} s\n".format(time.time() - ts))

    ts = time.time()

    si.track(r_start=np.array([0.0, 0.0, -0.15]),
             v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
             nsteps=15000,
             dt=1e-11)

    print("Tracking took {:.4f} s".format(time.time() - ts))

    with open('timing.txt', 'a') as outfile:
        outfile.write("Tracking took {:.4f} s\n".format(time.time() - ts))

    si.draw_geometry(show=True, filename='auto')
    si.save_geo_files()
