from dans_pymodules import *
import numpy as np
from scipy import constants as const
from mpl_toolkits.mplot3d import proj3d

# Try importing BEMPP
import bempp.api
# noinspection PyUnresolvedReferences
from bempp.api.shapes.shapes import __generate_grid_from_geo_string as generate_from_string

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

# Define the directions:
X = 0
Y = 1
Z = 2


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
                 debug=False,
                 outp_folder="",
                 **kwargs):

        # --- Program Variables -------------------------------------------------------------------------------------- #
        self._method = method  # Either analytical or numerical

        assert self._method in ["analytical", "numerical"], \
            "Spiral inflector method must be 'analytical' or 'numerical'!"

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
                                 "sigma": None,  # The "v-shape" of the electrodes (m)
                                 "ns": None,  # Resolution of the analytical solution (steps along beam trajectory s)
                                 "b_lim": np.deg2rad(np.array([0.0, 90.0])),  # Limits of the curvature
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
                                    }

        # --- Parameters used by the BEM++ potential and field calculation ------------------------------------------- #
        self._params_bempp = {"h": None,  # the desired mesh spacing for BEM++ mesh generation
                              "make_aperture": False,  # Make apertures at the exit and anetrnace
                              "aperture_params": {"thickness": None,
                                                  "radius": None,
                                                  "length": None,
                                                  "width": None,
                                                  "distance": None,
                                                  "voltage": 0.0},
                              "make_cylinder": False,  # Cylindrical housing
                              "cylinder_params": {"radius": None,
                                                  "zmin": None,
                                                  "zmax": None,
                                                  "voltage": 0.0}
                              }

        for key in self._params_bempp.keys():
            if key in kwargs.keys():
                self._params_bempp[key] = kwargs[key]

        self._variables_bempp = {"objects": {},  # A dictionary of objects (apertures, electrodes)
                                 "full mesh": None,  # Full mesh of the geometry
                                 "i": None,  # A running 3-index for mesh generation
                                 "f_space": None,  # BEM++ Function Space
                                 "operator": None,  # BEM++ Operator
                                 "grid_fun": None,  # BEM++ Grid Function
                                 "solution": None,  # The BEM++ solution
                                 "ef_itp": None,  # Interpolator object for E-Field
                                 }

        # --- Additional parameters used for particle tracking ------------------------------------------------------- #
        self._params_track = {}
        self._variables_track = {"trj_tracker": None,
                                 "shift": None}

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
        if self._method is "analytical":
            r = self.__generate_analytical_trajectory()
        elif self._method is "numerical":
            r = self.__generate_numerical_trajectory()
        return r

    def generate_geometry(self):
        if self._method is "analytical":
            r = self.__generate_analytical_geometry()
        elif self._method is "numerical":
            r = self.__generate_numerical_geometry()
        return r

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

    def calculate_efield(self,
                         limits=((None, None), (None, None), (None, None)),
                         res=0.002,
                         domain_decomp=(4, 4, 4),
                         overlap=0):
        """
        Calculates the E-Field from the BEM++ solution using the user defined cube or
        the cube corresponding to the cyclindrical outer boundary.
        :param limits: tuple, list or np.ndarray of shape (3, 2)
                       containing xmin, xmax, ymin, ymax, zmin, zmax
                       use None to use the individual limit from the electrode system.
        :param res: resolution of the 3D mesh
        :param domain_decomp: how many subdomains to use for calculation in the three directions x, y, z
                              Note: it can significantly increase computation speed to use more subdomains,
                              up to a point...
        :param overlap: overlap of the subdomains in cell numbers, does not have effect at the moment.
                        Note: There is a minimum overlap of one cell at overlap = 0
        :return:
        """

        limits = np.array(limits)

        if limits.shape != (3, 2):
            print("Wrong shape of limits: {}. "
                  "Must be ((xmin, xmax), (ymin, ymax), (zmin, zmax)) = (3, 2).".format(limits.shape))
            return 1

        if self._variables_bempp["solution"] is None:
            print("Please solve with BEM++ before calculating the E-Field")
            return 1

        _ts = time.time()

        sol = self._variables_bempp["solution"]
        fsp = self._variables_bempp["f_space"]

        # noinspection PyUnresolvedReferences
        all_vert = self._variables_bempp["full mesh"].leaf_view.vertices

        # get limits from electrodes
        limits_elec = np.array([[np.min(all_vert[i, :]), np.max(all_vert[i, :])] for i in range(3)])

        # replace None limits with electrode limits
        limits[np.where(limits is None)] = limits_elec[np.where(limits is None)]

        _n = np.array((limits[:, 1] - limits[:, 0]) / res, int)

        # Recalculate resolution to match integer n's
        _d = (limits[:, 1] - limits[:, 0]) / _n

        # Generate a full mesh to be indexed later
        _r = np.array([np.linspace(limits[i, 0], limits[i, 1], _n[i]) for i in range(3)])
        mesh = np.meshgrid(_r[X], _r[Y], _r[Z], indexing='ij')  # type: np.ndarray

        # Initialize potential array
        pot = np.zeros(mesh[0].shape)

        # Index borders (can be float)
        borders = np.array([np.linspace(0, _n[i], domain_decomp[i] + 1) for i in range(3)])

        start_idxs = np.array(borders[:, :-1], int) - overlap
        start_idxs[:, 0] = 0

        end_idxs = np.array(borders[:, 1:], int) + overlap
        end_idxs[:, -1] = np.array(borders[:, -1], int)

        # Print out domain information
        if self._debug:
            print("E-Field Calculation. "
                  "Grid spacings: ({:.4f}, {:.4f}, {:.4f}), number of meshes: {}".format(_d[0], _d[1], _d[2], _n))
            print("Number of Subdomains: {}, "
                  "Domain decomposition {}:".format(np.product(domain_decomp), domain_decomp))

            for i, dirs in enumerate(["x", "y", "z"]):
                print("{}: Indices {} to {}".format(dirs, start_idxs[i], end_idxs[i]))

        # Iterate over all the dimensions, calculate the subset of e-field
        domain_idx = 1
        for x1, x2 in zip(start_idxs[X], end_idxs[X]):
            for y1, y2 in zip(start_idxs[Y], end_idxs[Y]):
                for z1, z2 in zip(start_idxs[Z], end_idxs[Z]):

                    if self._debug:
                        print("[{}] Domain {}/{}, "
                              "Index Limits: x = ({}, {}), "
                              "y = ({}, {}), "
                              "z = ({}, {})".format(time.strftime('%H:%M:%S', time.gmtime(int(time.time() - _ts))),
                                                    domain_idx, np.product(domain_decomp), x1, x2, y1, y2, z1, z2))

                    grid_pts = np.vstack([_mesh[x1:x2, y1:y2, z1:z2].ravel() for _mesh in mesh])

                    # TODO: We can decide to omit certain regions here to speed up the
                    # TODO: process since the pot array was initialized as zeros... -DW
                    sl_pot = bempp.api.operators.potential.laplace.single_layer(fsp, grid_pts)
                    _pot = sl_pot * sol
                    pot[x1:x2, y1:y2, z1:z2] = _pot.reshape([x2 - x1, y2 - y1, z2 - z1])

                    domain_idx += 1

                    del grid_pts
                    del sl_pot
                    del _pot

        ex, ey, ez = np.gradient(pot, _d[X], _d[Y], _d[Z])

        del pot
        gc.collect()

        self._variables_track["ef_itp"] = Field("Spiral Inflector E-Field",
                                                dim=3,
                                                field={"x": RegularGridInterpolator(points=_r, values=-ex,
                                                                                    bounds_error=False, fill_value=0.0),
                                                       "y": RegularGridInterpolator(points=_r, values=-ey,
                                                                                    bounds_error=False, fill_value=0.0),
                                                       "z": RegularGridInterpolator(points=_r, values=-ez,
                                                                                    bounds_error=False, fill_value=0.0)
                                                       })

        return 0

    def draw_geometry(self, freq = 10, show=False, filename=None):

        if self._variables_analytic["geo"] is None or self._variables_analytic["trj_design"] is None:
            print("No geometry yet ... generating")
            self.generate_geometry()

        trj_design = self._variables_analytic["trj_design"]  # type: np.ndarray
        geo = self._variables_analytic["geo"]  # type: np.ndarray
        shift = self._variables_track["shift"]  # type: np.ndarray
        if shift is None:
            shift = np.zeros(3)

        fig = plt.figure()
        ax = Axes3D(fig)

        proj3d.persp_transformation = orthogonal_proj

        ns = self._params_analytic["ns"]  # type: int

        # Plot the beam trajectory
        ax.plot(trj_design[:, 0] + shift[0],
                trj_design[:, 1] + shift[1],
                trj_design[:, 2] + shift[2],
                color=colors[3], linewidth=3.0)

        # Plot the splines/edge lines
        for i in range(10):
            ax.plot(geo[i, :, 0] + shift[0],
                    geo[i, :, 1] + shift[1],
                    geo[i, :, 2] + shift[2],
                    color=colors[4], linewidth=2.0)

        # Lines between splines
        for i in range(ns):

            if (np.mod(i, freq) == 0) or (i == ns - 1):
                ax.plot([geo[2, i, 0] + shift[0], geo[3, i, 0] + shift[0]],
                        [geo[2, i, 1] + shift[1], geo[3, i, 1] + shift[1]],
                        [geo[2, i, 2] + shift[2], geo[3, i, 2] + shift[2]],
                        color=colors[1], linewidth=2.0)

                ax.plot([geo[7, i, 0] + shift[0], geo[8, i, 0] + shift[0]],
                        [geo[7, i, 1] + shift[1], geo[8, i, 1] + shift[1]],
                        [geo[7, i, 2] + shift[2], geo[8, i, 2] + shift[2]],
                        color=colors[1], linewidth=2.0)

                ax.plot([geo[0, i, 0] + shift[0], geo[4, i, 0] + shift[0]],
                        [geo[0, i, 1] + shift[1], geo[4, i, 1] + shift[1]],
                        [geo[0, i, 2] + shift[2], geo[4, i, 2] + shift[2]],
                        color=colors[2], linewidth=2.0)

                ax.plot([geo[4, i, 0] + shift[0], geo[1, i, 0] + shift[0]],
                        [geo[4, i, 1] + shift[1], geo[1, i, 1] + shift[1]],
                        [geo[4, i, 2] + shift[2], geo[1, i, 2] + shift[2]],
                        color=colors[2], linewidth=2.0)

                ax.plot([geo[5, i, 0] + shift[0], geo[9, i, 0] + shift[0]],
                        [geo[5, i, 1] + shift[1], geo[9, i, 1] + shift[1]],
                        [geo[5, i, 2] + shift[2], geo[9, i, 2] + shift[2]],
                        color=colors[2], linewidth=2.0)

                ax.plot([geo[9, i, 0] + shift[0], geo[6, i, 0] + shift[0]],
                        [geo[9, i, 1] + shift[1], geo[6, i, 1] + shift[1]],
                        [geo[9, i, 2] + shift[2], geo[6, i, 2] + shift[2]],
                        color=colors[2], linewidth=2.0)

        if self._variables_track["trj_tracker"] is not None:

            trj_tracked = self._variables_track["trj_tracker"]

            # Plot the tracked beam trajectory
            ax.plot(trj_tracked[:, 0] + shift[0],
                    trj_tracked[:, 1] + shift[1],
                    trj_tracked[:, 2] + shift[2],
                    color=colors[1], linewidth=3.0)

        plt.axis("equal")
        plt.xlabel("x")
        plt.ylabel("y")
        ax.set_xlim(-0.1, 0.1)
        ax.set_ylim(-0.1, 0.1)
        ax.set_zlabel("z")
        ax.set_aspect("equal")

        ax.view_init(elev=-180.0, azim=-90.0)

        if filename == 'dialog':
            print("Please select figure filename in dialog...")
            fd = FileDialog()
            filename = fd.get_filename('open')
        elif filename == 'auto':
            filename = os.path.join(self._outp_folder, "geometry.png")

        if filename is not None:
            try:
                fig.savefig(filename=filename)
            except Exception as ex:
                print("Something went wrong when trying to save the file: {}".format(ex))
                return 1

        if show:
            plt.show()

        return 0

    def draw_mesh(self):

        if self._variables_bempp["full mesh"] is not None:
            self._variables_bempp["full mesh"].plot()

    def generate_aperture_geo(self, electrode_type):

        if electrode_type == "top_aperture":
            i = 0
        elif electrode_type == "bottom_aperture":
            i = -1
        else:
            print("You must specify either 'top_aperture' or 'bottom_aperture'!")
            return 1

        # Check if all parameters are set
        abort_flag = False
        for key, item in self._params_bempp["aperture_params"].items():
            if item is None:
                print("Item 'aperture_params/{}' is not set in BEM++ parameters!".format(key))
                abort_flag = True
        if abort_flag:
            return 1

        h = self._params_bempp["h"]
        geo = self._variables_analytic["geo"]
        trj = self._variables_analytic["trj_design"]  # type: np.ndarray
        thickness = self._params_bempp["aperture_params"]["thickness"]
        radius = self._params_bempp["aperture_params"]["radius"]
        length = self._params_bempp["aperture_params"]["length"]
        width = self._params_bempp["aperture_params"]["width"]
        aperture_distance = self._params_bempp["aperture_params"]["distance"]
        voltage = self._params_bempp["aperture_params"]["voltage"]

        gmsh_str = """
h = {};
""".format(h * 1.5)

        # Note: These mid-vectors are swapped compared to the ones in the VB macro!
        mid_vec_a = (geo[4, i, :] - geo[9, i, :])
        mid_vec_b = (geo[8, i, :] - geo[7, i, :])

        mid_vec_a /= np.linalg.norm(mid_vec_a)
        mid_vec_b /= np.linalg.norm(mid_vec_b)

        norm_vec = np.cross(mid_vec_b, mid_vec_a)
        norm_vec /= np.linalg.norm(norm_vec)

        offset_a = norm_vec * aperture_distance
        offset_b = norm_vec * (aperture_distance + thickness)

        if i == 0:
            offset_a *= -1.0
            offset_b *= -1.0

        gmsh_str += """
Point(1) = {{ {}, {}, {}, h }};
Point(2) = {{ {}, {}, {}, h }};
Point(3) = {{ {}, {}, {}, h }};
Point(4) = {{ {}, {}, {}, h }};
Point(5) = {{ {}, {}, {}, h }};""".format(
            trj[i, 0] + offset_a[0],
            trj[i, 1] + offset_a[1],
            trj[i, 2] + offset_a[2],
            trj[i, 0] + offset_a[0] + mid_vec_a[0] * radius,
            trj[i, 1] + offset_a[1] + mid_vec_a[1] * radius,
            trj[i, 2] + offset_a[2] + mid_vec_a[2] * radius,
            trj[i, 0] + offset_a[0] + mid_vec_b[0] * radius,
            trj[i, 1] + offset_a[1] + mid_vec_b[1] * radius,
            trj[i, 2] + offset_a[2] + mid_vec_b[2] * radius,
            trj[i, 0] + offset_a[0] - mid_vec_a[0] * radius,
            trj[i, 1] + offset_a[1] - mid_vec_a[1] * radius,
            trj[i, 2] + offset_a[2] - mid_vec_a[2] * radius,
            trj[i, 0] + offset_a[0] - mid_vec_b[0] * radius,
            trj[i, 1] + offset_a[1] - mid_vec_b[1] * radius,
            trj[i, 2] + offset_a[2] - mid_vec_b[2] * radius)

        gmsh_str += """
Point(6) = {{ {}, {}, {}, h }};
Point(7) = {{ {}, {}, {}, h }};
Point(8) = {{ {}, {}, {}, h }};
Point(9) = {{ {}, {}, {}, h }};
Point(10) = {{ {}, {}, {}, h }};""".format(
            trj[i, 0] + offset_b[0],
            trj[i, 1] + offset_b[1],
            trj[i, 2] + offset_b[2],
            trj[i, 0] + offset_b[0] + mid_vec_a[0] * radius,
            trj[i, 1] + offset_b[1] + mid_vec_a[1] * radius,
            trj[i, 2] + offset_b[2] + mid_vec_a[2] * radius,
            trj[i, 0] + offset_b[0] + mid_vec_b[0] * radius,
            trj[i, 1] + offset_b[1] + mid_vec_b[1] * radius,
            trj[i, 2] + offset_b[2] + mid_vec_b[2] * radius,
            trj[i, 0] + offset_b[0] - mid_vec_a[0] * radius,
            trj[i, 1] + offset_b[1] - mid_vec_a[1] * radius,
            trj[i, 2] + offset_b[2] - mid_vec_a[2] * radius,
            trj[i, 0] + offset_b[0] - mid_vec_b[0] * radius,
            trj[i, 1] + offset_b[1] - mid_vec_b[1] * radius,
            trj[i, 2] + offset_b[2] - mid_vec_b[2] * radius)

        gmsh_str += """
    Point(11) = {{ {}, {}, {}, h }};
    Point(12) = {{ {}, {}, {}, h }};
    Point(13) = {{ {}, {}, {}, h }};
    Point(14) = {{ {}, {}, {}, h }};""".format(
            trj[i, 0] + offset_a[0] + mid_vec_a[0] * width / 2.0 + mid_vec_b[0] * length / 2.0,
            trj[i, 1] + offset_a[1] + mid_vec_a[1] * width / 2.0 + mid_vec_b[1] * length / 2.0,
            trj[i, 2] + offset_a[2] + mid_vec_a[2] * width / 2.0 + mid_vec_b[2] * length / 2.0,
            trj[i, 0] + offset_a[0] - mid_vec_a[0] * width / 2.0 + mid_vec_b[0] * length / 2.0,
            trj[i, 1] + offset_a[1] - mid_vec_a[1] * width / 2.0 + mid_vec_b[1] * length / 2.0,
            trj[i, 2] + offset_a[2] - mid_vec_a[2] * width / 2.0 + mid_vec_b[2] * length / 2.0,
            trj[i, 0] + offset_a[0] - mid_vec_a[0] * width / 2.0 - mid_vec_b[0] * length / 2.0,
            trj[i, 1] + offset_a[1] - mid_vec_a[1] * width / 2.0 - mid_vec_b[1] * length / 2.0,
            trj[i, 2] + offset_a[2] - mid_vec_a[2] * width / 2.0 - mid_vec_b[2] * length / 2.0,
            trj[i, 0] + offset_a[0] + mid_vec_a[0] * width / 2.0 - mid_vec_b[0] * length / 2.0,
            trj[i, 1] + offset_a[1] + mid_vec_a[1] * width / 2.0 - mid_vec_b[1] * length / 2.0,
            trj[i, 2] + offset_a[2] + mid_vec_a[2] * width / 2.0 - mid_vec_b[2] * length / 2.0)

        gmsh_str += """
    Point(15) = {{ {}, {}, {}, h }};
    Point(16) = {{ {}, {}, {}, h }};
    Point(17) = {{ {}, {}, {}, h }};
    Point(18) = {{ {}, {}, {}, h }};""".format(
            trj[i, 0] + offset_b[0] + mid_vec_a[0] * width / 2.0 + mid_vec_b[0] * length / 2.0,
            trj[i, 1] + offset_b[1] + mid_vec_a[1] * width / 2.0 + mid_vec_b[1] * length / 2.0,
            trj[i, 2] + offset_b[2] + mid_vec_a[2] * width / 2.0 + mid_vec_b[2] * length / 2.0,
            trj[i, 0] + offset_b[0] - mid_vec_a[0] * width / 2.0 + mid_vec_b[0] * length / 2.0,
            trj[i, 1] + offset_b[1] - mid_vec_a[1] * width / 2.0 + mid_vec_b[1] * length / 2.0,
            trj[i, 2] + offset_b[2] - mid_vec_a[2] * width / 2.0 + mid_vec_b[2] * length / 2.0,
            trj[i, 0] + offset_b[0] - mid_vec_a[0] * width / 2.0 - mid_vec_b[0] * length / 2.0,
            trj[i, 1] + offset_b[1] - mid_vec_a[1] * width / 2.0 - mid_vec_b[1] * length / 2.0,
            trj[i, 2] + offset_b[2] - mid_vec_a[2] * width / 2.0 - mid_vec_b[2] * length / 2.0,
            trj[i, 0] + offset_b[0] + mid_vec_a[0] * width / 2.0 - mid_vec_b[0] * length / 2.0,
            trj[i, 1] + offset_b[1] + mid_vec_a[1] * width / 2.0 - mid_vec_b[1] * length / 2.0,
            trj[i, 2] + offset_b[2] + mid_vec_a[2] * width / 2.0 - mid_vec_b[2] * length / 2.0)

        gmsh_str += """
    Circle(1) = { 2, 1, 3 };
    Circle(2) = { 3, 1, 4 };
    Circle(3) = { 4, 1, 5 };
    Circle(4) = { 5, 1, 2 };

    Circle(5) = { 7, 6, 8 };
    Circle(6) = { 8, 6, 9 };
    Circle(7) = { 9, 6, 10 };
    Circle(8) = { 10, 6, 7 };

    // Skip Lines 9, 10, 11 because I messed up -PW

    Line(11) = { 11, 12 };
    Line(12) = { 12, 13 };
    Line(13) = { 13, 14 };
    Line(14) = { 14, 11 };

    Line(15) = { 15, 16 };
    Line(16) = { 16, 17 };
    Line(17) = { 17, 18 };
    Line(18) = { 18, 15 };

    Line(19) = { 11, 15 };
    Line(20) = { 12, 16 };
    Line(21) = { 13, 17 };
    Line(22) = { 14, 18 };

    Line(23) = { 2, 7 };
    Line(24) = { 3, 8 };
    Line(25) = { 4, 9 };
    Line(26) = { 5, 10 };

    Line Loop(1) = { 1, 2, 3, 4 };
    Line Loop(2) = { 5, 6, 7, 8 };

    Line Loop(3) = { 11, 12, 13, 14 };
    Line Loop(4) = { 15, 16, 17, 18 };

    Line Loop(5) = { 1, 24, -5, -23 };
    Line Loop(6) = { 2, 25, -6, -24 };
    Line Loop(7) = { 3, 26, -7, -25 };
    Line Loop(8) = { 4, 23, -8, -26 };

    Line Loop(9) = { 11, 20, -15, -19 };
    Line Loop(10) = { 12, 21, -16, -20 };
    Line Loop(11) = { 13, 22, -17, -21 };
    Line Loop(12) = { 14, 19, -18, -22 };

    Plane Surface(1) = { 1, 3 };
    Plane Surface(2) = { 2, 4 };

    Ruled Surface(3) = { 5 };
    Ruled Surface(4) = { 6 };
    Ruled Surface(5) = { 7 };
    Ruled Surface(6) = { 8 };

    Plane Surface(7) = { 9 };
    Plane Surface(8) = { 10 };
    Plane Surface(9) = { 11 };
    Plane Surface(10) = { 12 };
    """

        self._variables_bempp["objects"][electrode_type] = {"gmsh_str": gmsh_str, "voltage": voltage}

        if self._debug:
            print(gmsh_str)

        return 0

    def generate_cylinder_geo(self):

        # Check if all parameters are set
        abort_flag = False
        for key, item in self._params_bempp["cylinder_params"].items():
            if item is None:
                print("Item 'aperture_params/{}' is not set in BEM++ parameters!".format(key))
                abort_flag = True
        if abort_flag:
            return 1

        h = self._params_bempp["h"]
        radius = self._params_bempp["cylinder_params"]["radius"]
        zmin = self._params_bempp["cylinder_params"]["zmin"]
        zmax = self._params_bempp["cylinder_params"]["zmax"]
        voltage = self._params_bempp["cylinder_params"]["voltage"]
        electrode_type = "cylinder"

        gmsh_str = """
h = {};
""".format(h * 5)

        gmsh_str += """
Point(1) = {{ {}, {}, {}, h }};
Point(2) = {{ {}, {}, {}, h }};
Point(3) = {{ {}, {}, {}, h }};
Point(4) = {{ {}, {}, {}, h }};
Point(5) = {{ {}, {}, {}, h }};""".format(0.0, 0.0, zmin,
                                          radius, 0.0, zmin,
                                          0.0, radius, zmin,
                                          -radius, 0.0, zmin,
                                          0.0, -radius, zmin)

        gmsh_str += """
Point(6) = {{ {}, {}, {}, h }};
Point(7) = {{ {}, {}, {}, h }};
Point(8) = {{ {}, {}, {}, h }};
Point(9) = {{ {}, {}, {}, h }};
Point(10) = {{ {}, {}, {}, h }};""".format(0.0, 0.0, zmax,
                                           radius, 0.0, zmax,
                                           0.0, radius, zmax,
                                           -radius, 0.0, zmax,
                                           0.0, -radius, zmax)

        gmsh_str += """
Circle(1) = { 2, 1, 3 };
Circle(2) = { 3, 1, 4 };
Circle(3) = { 4, 1, 5 };
Circle(4) = { 5, 1, 2 };

Circle(5) = { 7, 6, 8 };
Circle(6) = { 8, 6, 9 };
Circle(7) = { 9, 6, 10 };
Circle(8) = { 10, 6, 7 };

Line(23) = { 2, 7 };
Line(24) = { 3, 8 };
Line(25) = { 4, 9 };
Line(26) = { 5, 10 };

Line Loop(1) = { 1, 2, 3, 4 };
Line Loop(2) = { 5, 6, 7, 8 };

Line Loop(5) = { 1, 24, -5, -23 };
Line Loop(6) = { 2, 25, -6, -24 };
Line Loop(7) = { 3, 26, -7, -25 };
Line Loop(8) = { 4, 23, -8, -26 };

Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };

Ruled Surface(3) = { 5 };
Ruled Surface(4) = { 6 };
Ruled Surface(5) = { 7 };
Ruled Surface(6) = { 8 };
"""

        self._variables_bempp["objects"][electrode_type] = {"gmsh_str": gmsh_str, "voltage": voltage}

        if self._debug:
            print(gmsh_str)

        return 0

    def export_electrode_geometry(self, fname="electrode_macro.ivb"):

        geo = self._variables_analytic["geo"] * 100.0  # Fix scaling in inventor

        # Generate text for Inventor macro
        header_text = """Sub CreateSpiralElectrode()
    Dim oApp As Application
    Set oApp = ThisApplication

    ' Get a reference to the TransientGeometry object.
    Dim tg As TransientGeometry
    Set tg = oApp.TransientGeometry

    Dim oPart As PartDocument
    Dim oCompDef As PartComponentDefinition
    Dim oSketch As Sketch3D
    Dim oSpline As SketchSplines3D
    Dim vertexCollection1 As ObjectCollection
    Dim vertexCollection2 As ObjectCollection
    Dim vertexCollection3 As ObjectCollection
    Dim vertexCollection4 As ObjectCollection
    Dim vertexCollection5 As ObjectCollection
    Dim oLine As SketchLines3D
    Dim number_of_points As Long
    Dim loft_section_index As Long
    Dim frequency As Integer: frequency = 10
    Dim oLoftDef As LoftDefinition
    Dim oLoftSections As ObjectCollection
    Dim spiral_electrode As LoftFeature
"""
        electrode_texts = ["", ""]
        for j in range(2):
            electrode_texts[j] += """
    Set oPart = oApp.Documents.Add(kPartDocumentObject, , True)
    oPart.UnitsOfMeasure.LengthUnits = kMeterLengthUnits

    Set oCompDef = oPart.ComponentDefinition

"""
            # Loop over the five splines
            spline_texts = ["", "", "", "", ""]
            for k in range(5):
                spline_texts[k] += """
    Set oSketch = oCompDef.Sketches3D.Add
    Set oSpline = oSketch.SketchSplines3D
    Set vertexCollection{} = oApp.TransientObjects.CreateObjectCollection(Null)

""".format(k + 1)

                for i in range(self._params_analytic["ns"]):
                    spline_texts[k] += "    "
                    spline_texts[k] += "Call vertexCollection{}.Add(tg.CreatePoint({:.6f}, {:.6f}, {:.6f}))".format(
                        k + 1,
                        geo[k + (5 * j), i, 0],
                        geo[k + (5 * j), i, 1],
                        geo[k + (5 * j), i, 2])
                    spline_texts[k] += "\n"

                spline_texts[k] += """
    Call oSpline.Add(vertexCollection{})

""".format(k + 1)

                electrode_texts[j] += spline_texts[k]

            electrode_texts[j] += """
    ' Find out total number of points in single spline
    number_of_points = vertexCollection1.Count

    ' Container holding the loft sections (rectangles)
    Set oLoftSections = oApp.TransientObjects.CreateObjectCollection

    For i = 0 To number_of_points - 2

        If i Mod frequency = 0 Then
            Set oSketch = oCompDef.Sketches3D.Add
            Set oLine = oSketch.SketchLines3D

            Call oLine.AddByTwoPoints(vertexCollection1.Item(i + 1), vertexCollection5.Item(i + 1))
            Call oLine.AddByTwoPoints(vertexCollection2.Item(i + 1), vertexCollection5.Item(i + 1))
            Call oLine.AddByTwoPoints(vertexCollection3.Item(i + 1), vertexCollection4.Item(i + 1))
            Call oLine.AddByTwoPoints(vertexCollection3.Item(i + 1), vertexCollection1.Item(i + 1))
            Call oLine.AddByTwoPoints(vertexCollection4.Item(i + 1), vertexCollection2.Item(i + 1))

            loft_section_index = i / frequency + 6

            Call oLoftSections.Add(oCompDef.Sketches3D.Item(loft_section_index).Profiles3D.AddOpen)
        End If
    Next i

    ' Make a new 3D sketch for the rectangle on top
    Set oSketch = oCompDef.Sketches3D.Add
    Set oLine = oSketch.SketchLines3D

    Call oLine.AddByTwoPoints(vertexCollection1.Item(number_of_points), vertexCollection5.Item(number_of_points))
    Call oLine.AddByTwoPoints(vertexCollection2.Item(number_of_points), vertexCollection5.Item(number_of_points))
    Call oLine.AddByTwoPoints(vertexCollection3.Item(number_of_points), vertexCollection4.Item(number_of_points))
    Call oLine.AddByTwoPoints(vertexCollection3.Item(number_of_points), vertexCollection1.Item(number_of_points))
    Call oLine.AddByTwoPoints(vertexCollection4.Item(number_of_points), vertexCollection2.Item(number_of_points))

    Call oLoftSections.Add(oCompDef.Sketches3D.Item(loft_section_index + 1).Profiles3D.AddOpen)

    ' Do more loft stuff
    Set oLoftDef = oCompDef.Features.LoftFeatures.CreateLoftDefinition(oLoftSections, kJoinOperation)

    Call oLoftDef.LoftRails.Add(oCompDef.Sketches3D.Item(1).Profiles3D.AddOpen)
    Call oLoftDef.LoftRails.Add(oCompDef.Sketches3D.Item(2).Profiles3D.AddOpen)
    Call oLoftDef.LoftRails.Add(oCompDef.Sketches3D.Item(3).Profiles3D.AddOpen)
    Call oLoftDef.LoftRails.Add(oCompDef.Sketches3D.Item(4).Profiles3D.AddOpen)
    Call oLoftDef.LoftRails.Add(oCompDef.Sketches3D.Item(5).Profiles3D.AddOpen)

    Set spiral_electrode = oCompDef.Features.LoftFeatures.Add(oLoftDef)

"""

        footer_text = """
    ' oPart.UnitsOfMeasure.LengthUnits = kMillimeterLengthUnits

    ThisApplication.ActiveView.Fit

End Sub
"""

        with open(os.path.join(self._outp_folder, filename), "w") as outfile:

            outfile.write(header_text + electrode_texts[0] + electrode_texts[1] + footer_text)

        print("Done!")

    def export_aperture_geometry(self, fname="aperture_macro.ivb"):

        geo = self._variables_analytic["geo"] * 100.0  # Scaling for inventor
        trj = self._variables_analytic["trj_design"]  # type: np.ndarray
        thickness = self._params_bempp["aperture_params"]["thickness"]
        radius = self._params_bempp["aperture_params"]["radius"]
        length = self._params_bempp["aperture_params"]["length"]
        width = self._params_bempp["aperture_params"]["width"]
        aperture_distance = self._params_bempp["aperture_params"]["distance"]
        voltage = self._params_bempp["aperture_params"]["voltage"]

        aperture_string = """Sub createApertures()
    Dim oApp As Application
    Set oApp = ThisApplication

    ' Get a reference to the TransientGeometry object.
    Dim tg As TransientGeometry
    Set tg = oApp.TransientGeometry

    Dim oPart As PartDocument
    Dim oCompDef As PartComponentDefinition
    Dim oExtrudeDef As ExtrudeDefinition
    Dim xyPlane As Inventor.WorkPlane
    Dim aperture As ExtrudeFeature
    Dim oWorkPlaneSketch As Sketch3D
    Dim oPoints As SketchPoints3D
    Dim oWorkPlane As Inventor.WorkPlane

    Dim sketch As Inventor.PlanarSketch
    Dim sketch_cut As Inventor.PlanarSketch

    ' Some variables from the python code
    Dim aperture_thickness As Double: aperture_thickness = {}
    Dim aperture_radius As Double: aperture_radius = {}
    Dim hole_length As Double: hole_length = {}
    Dim hole_width As Double: hole_width = {}
""".format(thickness, radius, length, width)

        for i in [0, -1]:

            mid_vec_a = (geo[8, i, :] - geo[7, i, :])
            mid_vec_b = (geo[4, i, :] - geo[9, i, :])

            norm_vec = np.cross(mid_vec_b, mid_vec_a)
            norm_vec /= np.linalg.norm(norm_vec)

            mid_vec_a += trj[i, :]
            mid_vec_b += trj[i, :]

            offset = norm_vec * aperture_distance

            if i == -1:
                offset *= -1.0

            aperture_string += """
    Set oPart = oApp.Documents.Add(kPartDocumentObject, , True)
    oPart.UnitsOfMeasure.LengthUnits = kMeterLengthUnits

    Set oCompDef = oPart.ComponentDefinition
    Set xyPlane = oCompDef.WorkPlanes.Item(3)

    Set oWorkPlaneSketch = oCompDef.Sketches3D.Add()
    Set oPoints = oWorkPlaneSketch.SketchPoints3D
    Call oPoints.Add(tg.CreatePoint({}, {}, {})) ' This will be the origin of the part
    Call oPoints.Add(tg.CreatePoint({}, {}, {}))
    Call oPoints.Add(tg.CreatePoint({}, {}, {}))
""".format(trj[i, 0] + offset[0], trj[i, 1] + offset[1],
           trj[i, 2] + offset[2], mid_vec_a[0] + offset[0], mid_vec_a[1] + offset[1],
           mid_vec_a[2] + offset[2], mid_vec_b[0] + offset[0], mid_vec_b[1] + offset[1], mid_vec_b[2] + offset[2])

            if i == 0:
                extrude_dir = "kNegativeExtentDirection"
            elif i == -1:
                extrude_dir = "kPositiveExtentDirection"

            aperture_string += """
    Set oWorkPlane = oCompDef.WorkPlanes.AddByThreePoints(oPoints.Item(1), oPoints.Item(2), oPoints.Item(3))

    ' Create a sketch to create a cylindrical aperture
    Set sketch = oCompDef.Sketches.Add(oWorkPlane, False)

    Call sketch.SketchCircles.AddByCenterRadius(tg.CreatePoint2d(0#, 0#), aperture_radius)
    Call sketch.Profiles.AddForSolid

    Set oExtrudeDef = oCompDef.Features.ExtrudeFeatures.CreateExtrudeDefinition(sketch.Profiles.Item(1), kNewBodyOperation)
    Call oExtrudeDef.SetDistanceExtent(aperture_thickness, {})

    Set aperture = oCompDef.Features.ExtrudeFeatures.Add(oExtrudeDef)

    ' Create a sketch that will cut the hole in the cylinder
    Set sketch_cut = oCompDef.Sketches.Add(oWorkPlane, False)

    Call sketch_cut.SketchLines.AddAsTwoPointCenteredRectangle(tg.CreatePoint2d(0#, 0#), tg.CreatePoint2d(hole_length / 2#, hole_width / 2#))
    Call sketch_cut.Profiles.AddForSolid

    Set oExtrudeDef = oCompDef.Features.ExtrudeFeatures.CreateExtrudeDefinition(sketch_cut.Profiles.Item(1), kCutOperation)
    Call oExtrudeDef.SetDistanceExtent(aperture_thickness, {})

    Call oCompDef.Features.ExtrudeFeatures.Add(oExtrudeDef)
""".format(extrude_dir, extrude_dir)

        aperture_string += """
    ' Fin
    ' oPart.UnitsOfMeasure.LengthUnits = kMillimeterLengthUnits
    ThisApplication.ActiveView.Fit

End Sub
"""

        with open(os.path.join(self._outp_folder, fname), "w") as outfile:
            outfile.writelines(aperture_string)

    def __generate_analytical_trajectory(self):

        if not self._initialized:
            self.initialize()

        print("Calculating Design Trajectory... ", end="")

        h = self._variables_analytic["height"]
        tilt = self._params_analytic["tilt"]  # type: float

        self._variables_analytic["kp"] = np.tan(np.deg2rad(tilt))  # Tilt parameter
        self._variables_analytic["k"] = ((h / self._variables_analytic["r_cyc"]) + self._variables_analytic["kp"]) / 2.0

        cp = self._variables_analytic["c+"] = (2.0 * self._variables_analytic["k"] + 1.0)
        cm = self._variables_analytic["c-"] = -(2.0 * self._variables_analytic["k"] - 1.0)

        # --- Trajectory coordinates --- #
        _x = +0.5 * h * ((2.0 / (1.0 - (4.0 * (self._variables_analytic["k"] ** 2.0)))) -
                         (np.cos(cp * self._variables_analytic["b"]) / cp) - np.cos(
            -cm * self._variables_analytic["b"]) / cm)

        _y = -0.5 * h * (np.sin(cp * self._variables_analytic["b"]) / cp +
                         np.sin(-cm * self._variables_analytic["b"]) / cm)

        _z = - h * (1.0 - np.sin(self._variables_analytic["b"]))

        self._variables_analytic["trj_design"] = np.array([_x, _y, _z]).T

        # Rotation/flip
        if not ((self._variables_analytic["bf_design"] < 0.0) ^ (self._params_analytic["ion"].q() < 0.0)):
            if self._debug:
                print("Flipping direction of cyclotron motion...", end="")
            self._variables_analytic["trj_design"][:, 1] = -self._variables_analytic["trj_design"][:, 1]

        print("Done!")
        if self._debug:
            print("Design Trajectory:")
            print(self._variables_analytic["trj_design"])
            print("")

        return self._variables_analytic["trj_design"]

    def __generate_analytical_geometry(self, **kwargs):
        """
        The process of generating the geometry is as follows:
        Create the inner and outer surface of the spiral electrodes, shift the inside edges according to sigma,
        put everything together in one array.
        :return:
        """

        if self._variables_analytic["trj_design"] is None:
            print("No analytical design trajectory yet, generating...")
            self.generate_design_trajectory()

        print("Generating analytical geometry... ", end="")

        geos = []
        _ion = self._params_analytic["ion"]  # type: IonSpecies
        ns = self._params_analytic["ns"]  # type: int
        gap = self._params_analytic["gap"]  # type: float
        sigma = self._params_analytic["sigma"]  # type: float
        kp = self._variables_analytic["kp"]  # type: float
        b = self._variables_analytic["b"]  # type: np.ndarray
        cp = self._variables_analytic["c+"]  # type: float
        cm = self._variables_analytic["c-"]  # type: float
        trj_design = self._variables_analytic["trj_design"]  # type: np.ndarray

        for thickness in [0.0, self._params_analytic["dx"]]:

            end_distance = 2.0 * thickness + gap  # End to end distance of electrodes in (m)
            aspect_ratio = 2.5 * (gap / end_distance)

            # Distance between electrodes at inflection angle theta
            d = end_distance * np.ones(ns) / (np.sqrt(1.0 + ((kp ** 2.0) * (np.sin(b)) ** 2.0)))

            # Save inner gap size vs deflection angle as class variable
            if thickness == 0.0:
                self._variables_analytic["d"] = d

            # x-component of the velocity vector
            vx = np.array(0.5 * _ion.v_m_per_s() * (np.sin(cp * b) + np.sin(cm * b)))

            # y-component of the velocity vector
            vy = np.array(-0.5 * _ion.v_m_per_s() * (np.cos(cp * b) - np.cos(cm * b)))

            # Rotation/flip
            if not ((self._variables_analytic["bf_design"] > 0.0) ^ (_ion.q() > 0.0)):
                if self._debug:
                    print("Flipping direction of cyclotron motion...", end="")
                vy = -vy

            v2 = np.sqrt((vx ** 2.0) + (vy ** 2.0))  # xy-magnitude of the velocity
            vz = np.sqrt((_ion.v_m_per_s() ** 2.0) - (v2 ** 2.0) + 0.0j)  # z-component of the velocity

            # Checks for imaginary components of the z-component
            for i in range(ns):
                if np.imag(vz[i]) != 0.0:
                    vz[i] = 0.0 + 0.0j  # Redefines that element as 0

            vz = np.real(vz)
            v3 = np.sqrt((vx ** 2.0) + (vy ** 2.0) + (vz ** 2.0))  # 3-d magnitude of the velocity vector (Should = v)

            # Save vx, vy, vz in as class variable in same format as trj_design
            self._variables_analytic["v_design"] = np.array([vx, vy, vz]).T

            # Construction of the normalized vectors of the optical coordinate system
            v_path = np.transpose(np.array([vx, vy, vz]))
            v_optical = np.zeros((3, ns, 3))  # [h, u, v]

            for j in range(ns):

                for k in range(3):
                    v_optical[2, j, k] = v_path[j, k] / v3[j]

                if v2[j] == 0:
                    v_optical[0, j, :] = np.array([0, -1, 0])
                    v_optical[1, j, :] = np.array([-1, 0, 0])

                elif v2[j] > 0:
                    v_optical[0, j, :] = np.array([vy[j] / v2[j], -vx[j] / v2[j], 0])
                    v_optical[1, j, :] = np.array(
                        [-(vx[j] * vz[j]) / (v3[j] * v2[j]), -(vy[j] * vz[j]) / (v3[j] * v2[j]),
                         ((v3[j] ** 2) - (vz[j] ** 2)) / (v3[j] * v2[j])])

            # Rotation to the tilt angle if atilt > 0 || atilt < 0
            v_rh = np.copy(v_optical)  # "Rotated" or "right-handed" optical coordinate system vectors [hr, ur, v]

            # Note: the theta in the commented eq. below is not the same theta as in the code
            nemo = np.arctan(kp * np.sin(b))

            for i in range(ns):
                for j in range(3):
                    v_rh[0, i, j] = ((np.cos(nemo[i])) * v_optical[0, i, j]) - (np.sin(nemo[i])) * v_optical[1, i, j]
                    v_rh[1, i, j] = ((np.cos(nemo[i])) * v_optical[1, i, j]) + (np.sin(nemo[i])) * v_optical[0, i, j]

            # Turn track of unit vectors
            t1 = np.arange(0, 5, 0.01)
            v_er = np.zeros((3, 3, np.size(t1)))

            for i in range(3):
                for j in range(3):
                    for k in range(np.size(t1)):
                        v_er[i, j, k] = trj_design[ns - 1, j] + t1[k] * v_rh[i, ns - 1, j]

            # Construction of the electrodes
            edge_lines = np.zeros((5, ns, 3))

            xi = 0.5 * aspect_ratio * end_distance

            for i in range(ns):
                for j in range(3):
                    edge_lines[0, i, j] = trj_design[i, j] + 0.5 * d[i] * v_rh[1, i, j] + xi * v_rh[0, i, j]
                    edge_lines[1, i, j] = trj_design[i, j] + 0.5 * d[i] * v_rh[1, i, j] - xi * v_rh[0, i, j]
                    edge_lines[2, i, j] = trj_design[i, j] - 0.5 * d[i] * v_rh[1, i, j] + xi * v_rh[0, i, j]
                    edge_lines[3, i, j] = trj_design[i, j] - 0.5 * d[i] * v_rh[1, i, j] - xi * v_rh[0, i, j]

            geos.append(edge_lines)

        # Apply the v-shape modification:
        diff = (geos[0][0, :, :] - geos[1][0, :, :])  # Difference between ext and int
        diff_norm = (np.sqrt(diff[:, 0] ** 2 + diff[:, 1] ** 2 + diff[:, 2] ** 2))  # Magnitude of diff
        diff_hat = np.zeros(np.shape(diff))  # Initialize normalized vector

        for i in range(ns):
            diff_hat[i, :] = diff[i, :] / diff_norm[i]  # Calculate each normalized vector

        geo = np.zeros([10, ns, 3])  # Initialize geo array

        # Upper Electrode
        geo[0, :, :] = geos[0][0, :, :] + diff_hat * sigma
        geo[1, :, :] = geos[0][1, :, :] + diff_hat * sigma
        geo[2, :, :] = geos[1][0, :, :] + diff_hat * sigma  # Should reduce the thickness
        geo[3, :, :] = geos[1][1, :, :] + diff_hat * sigma  # Should reduce the thickness
        geo[4, :, :] = 0.5 * (geos[0][0, :, :] + geos[0][1, :, :])

        # Lower Electrode
        geo[5, :, :] = geos[0][2, :, :] + diff_hat * sigma
        geo[6, :, :] = geos[0][3, :, :] + diff_hat * sigma
        geo[7, :, :] = geos[1][2, :, :]
        geo[8, :, :] = geos[1][3, :, :]
        geo[9, :, :] = 0.5 * (geos[0][2, :, :] + geos[0][3, :, :])

        self._variables_analytic["geo"] = geo

        print("Done!")

        return self._variables_analytic["geo"]

    def generate_meshed_model(self, apertures=None, cylinder=None):

        if apertures is not None:
            self._params_bempp["make_aperture"] = apertures

        if cylinder is not None:
            self._params_bempp["make_cylinder"] = cylinder

        if self._variables_analytic["geo"] is None:
            print("No geometry generated yet... starting now...")
            self.generate_geometry()

        err = 0

        err += self.generate_spiral_electrode_geo("anode")
        err += self.generate_spiral_electrode_geo("cathode")

        if self._params_bempp["make_aperture"]:
            err += self.generate_aperture_geo("top_aperture")
            err += self.generate_aperture_geo("bottom_aperture")

        if self._params_bempp["make_cylinder"]:
            err += self.generate_cylinder_geo()

        if err > 0:
            print("Error occured during geo generation!")
            return 1

        if generate_from_string is not None:

            # Initialize empty arrays of the correct shape (3 x n)
            vertices = np.zeros([3, 0])
            elements = np.zeros([3, 0])
            vertex_counter = 0
            domain_index = 0
            domains = np.zeros([0], int)

            for name, electrode in self._variables_bempp["objects"].items():

                # Create the mesh of each electrode
                electrode["mesh"] = generate_from_string(electrode["gmsh_str"])
                electrode["domain"] = domain_index

                # Add it to the full model
                _vertices = electrode["mesh"].leaf_view.vertices
                _elements = electrode["mesh"].leaf_view.elements

                vertices = np.concatenate((vertices, _vertices), axis=1)
                elements = np.concatenate((elements, _elements + vertex_counter), axis=1)
                domains = np.concatenate((domains, domain_index * np.ones(_elements.shape[1], int)))

                # Increase the running counters
                vertex_counter += _vertices.shape[1]
                domain_index += 1

            # noinspection
            self._variables_bempp["full mesh"] = bempp.api.grid.grid_from_element_data(vertices, elements, domains)

            if self._debug:
                self._variables_bempp["full mesh"].plot()

        else:
            print("BEMPP could not be loaded, the geo strings were generated, but that's as far as it goes...")
            return 1

    def generate_spiral_electrode_geo(self, electrode_type):

        abort_flag = False
        for key, item in self._params_bempp.items():
            if item is None:
                print("Item {} is not set in BEM++ parameters!".format(key))
                abort_flag = True
        if abort_flag:
            return 1

        if self._variables_analytic["geo"] is None:
            print("No analytic geometry generated yet... starting now...")
            self.generate_geometry()

        self._variables_bempp["i"] = [1, 1, 1]

        voltage = 0.0

        if electrode_type is "anode":
            k = 0
            invert = True
            voltage = self._params_analytic["volt"]
        elif electrode_type is "cathode":
            k = 5
            invert = False
            voltage = -self._params_analytic["volt"]

        geo = self._variables_analytic["geo"]  # type: np.ndarray
        h = self._params_bempp["h"]
        gmsh_str = ""

        ly = geo.shape[1]

        for i in range(ly):

            for j in [0, 2, 3, 1, 4]:

                if i < 3 or i > ly - 4:

                    gmsh_str += self._gmsh_point(geo[j + k, i, :], h=h)

                else:

                    gmsh_str += self._gmsh_point(geo[j + k, i, :])

                self._variables_bempp["i"][0] += 1

                if j != 0:
                    gmsh_str += self._gmsh_line(self._variables_bempp["i"][0] - 2,
                                                self._variables_bempp["i"][0] - 1)
                    self._variables_bempp["i"][1] += 1

                if j == 4:

                    gmsh_str += self._gmsh_line(self._variables_bempp["i"][0] - 1,
                                                self._variables_bempp["i"][0] - 5)
                    self._variables_bempp["i"][1] += 1

                    if i == 0:  # Electrode end #1

                        gmsh_str += self._gmsh_surface([-(self._variables_bempp["i"][1] - 5),
                                                        -(self._variables_bempp["i"][1] - 4),
                                                        -(self._variables_bempp["i"][1] - 3),
                                                        -(self._variables_bempp["i"][1] - 2),
                                                        -(self._variables_bempp["i"][1] - 1)],
                                                       'Plane', invert)
                    if i == ly - 1:  # Electrode end #2

                        gmsh_str += self._gmsh_surface([(self._variables_bempp["i"][1] - 8),
                                                        (self._variables_bempp["i"][1] - 6),
                                                        (self._variables_bempp["i"][1] - 4),
                                                        (self._variables_bempp["i"][1] - 2),
                                                        (self._variables_bempp["i"][1] - 1)],
                                                       'Plane', invert)
                if i > 0:

                    gmsh_str += self._gmsh_line(self._variables_bempp["i"][0] - 1,
                                                self._variables_bempp["i"][0] - 6)
                    self._variables_bempp["i"][1] += 1

                    if i == 1:

                        if j == 2:

                            gmsh_str += self._gmsh_surface([self._variables_bempp["i"][1] - 8,
                                                            -(self._variables_bempp["i"][1] - 1),
                                                            -(self._variables_bempp["i"][1] - 2),
                                                            self._variables_bempp["i"][1] - 3],
                                                           'Ruled', invert)
                        elif j == 3:

                            gmsh_str += self._gmsh_surface([self._variables_bempp["i"][1] - 9,
                                                            -(self._variables_bempp["i"][1] - 1),
                                                            -(self._variables_bempp["i"][1] - 2),
                                                            self._variables_bempp["i"][1] - 3],
                                                           'Ruled', invert)
                        elif j == 1:

                            gmsh_str += self._gmsh_surface([self._variables_bempp["i"][1] - 10,
                                                            -(self._variables_bempp["i"][1] - 1),
                                                            -(self._variables_bempp["i"][1] - 2),
                                                            self._variables_bempp["i"][1] - 3],
                                                           'Ruled', invert)
                        elif j == 4:

                            gmsh_str += self._gmsh_surface([self._variables_bempp["i"][1] - 12,
                                                            -(self._variables_bempp["i"][1] - 1),
                                                            -(self._variables_bempp["i"][1] - 3),
                                                            self._variables_bempp["i"][1] - 4],
                                                           'Ruled', invert)
                    elif i > 1:

                        if j == 0:

                            if i == 2:
                                c = 0
                            else:
                                c = 1

                            gmsh_str += self._gmsh_surface([(self._variables_bempp["i"][1] - 12 - c),
                                                            (self._variables_bempp["i"][1] - 2),
                                                            -(self._variables_bempp["i"][1] - 3),
                                                            -(self._variables_bempp["i"][1] - 11)],
                                                           'Ruled', invert)
                        elif j != 0 and j != 4:

                            gmsh_str += self._gmsh_surface([self._variables_bempp["i"][1] - 12,
                                                            -(self._variables_bempp["i"][1] - 1),
                                                            -(self._variables_bempp["i"][1] - 2),
                                                            self._variables_bempp["i"][1] - 3],
                                                           'Ruled', invert)
                        elif j == 4:

                            gmsh_str += self._gmsh_surface([self._variables_bempp["i"][1] - 13,
                                                            -(self._variables_bempp["i"][1] - 1),
                                                            -(self._variables_bempp["i"][1] - 3),
                                                            self._variables_bempp["i"][1] - 4],
                                                           'Ruled', invert)

                            if i == ly - 1:
                                gmsh_str += self._gmsh_surface([(self._variables_bempp["i"][1] - 12),
                                                                (self._variables_bempp["i"][1] - 1),
                                                                -(self._variables_bempp["i"][1] - 2),
                                                                -(self._variables_bempp["i"][1] - 10)],
                                                               'Ruled', invert)

        self._variables_bempp["objects"][electrode_type] = {"gmsh_str": gmsh_str, "voltage": voltage}

        if self._debug:
            print(gmsh_str)

        return 0

    def get_parameter(self, key):

        if key in self._params_analytic.keys():

            return self._params_analytic[key]

        elif key in self._variables_analytic.keys():

            return self._variables_analytic[key]

        else:
            return None

    def _gmsh_point(self, pos, h=None):
        """
        returns the proper string for GMSH *.geo files point generation
        :param pos:
        :param h:
        :return:
        """

        if h is None:

            h = self._params_bempp["h"]

        return "Point({}) = {{ {}, {}, {}, {} }};\n".format(self._variables_bempp["i"][0], pos[0], pos[1], pos[2], h)

    def _gmsh_line(self, point_a, point_b):
        """
        returns the proper string for GMSH *.geo files line generation
        :param point_a:
        :param point_b:
        :return:
        """

        return "Line({}) = {{ {}, {} }};\n".format(self._variables_bempp["i"][1], point_a, point_b)

    def _gmsh_surface(self, lines, surface_type="Plane", invert=False):
        """
        returns the proper string for GMSH *.geo files surface generation
        :param surface_type: either 'Plane' or 'Ruled'
        :return:
        """

        if invert:
            c = 1
        else:
            c = -1

        return_str = "Line Loop({}) = {{".format(self._variables_bempp["i"][2])

        for line in lines[:-1]:
            return_str += " {},".format(c * line)

        return_str += " {} }};\n".format(c * lines[-1])
        return_str += "{} Surface({}) = {{ {} }};\n".format(surface_type,
                                                            self._variables_bempp["i"][2],
                                                            self._variables_bempp["i"][2])
        self._variables_bempp["i"][2] += 1

        return return_str

    def initialize(self, **kwargs):

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

        if si._initialized:
            si.initialize()

    def optimize_fringe(self, initial_guess=(None, None), maxiter=10, tol=1e-1, res=0.002):
        """
        This function optimizes the length of the spiral inflector to adjust for fringe fields.
        :param initial_guess: tuple, list or ndarray containing angle adjustment for entrance and exit
                              units: degrees
        :param maxiter: Maximum number of iterations to find the optimum entrance and exit adjustment before giving up
        :param tol: maximum deviation tolerated for a successful solution (degrees)
        :param res:
        :return:
        """

        print("Starting the optimization process...")

        # Start with an initial guess for bmin, bmax. If None, use 0.0
        # TODO: Update to use current internal adjustment, if present -DW
        _db_init = np.array(initial_guess)
        _db_init[np.where(_db_init is None)] = 0.0
        _db = _db_init[:]

        # Half the length of the cube for electric field calculation
        _hcl = (self._params_analytic["gap"] + 2.0 * self._params_analytic["dx"])
        # _ion = self._params_analytic["ion"]  # type: IonSpecies
        x_axis = Vector([1.0, 0.0, 0.0])
        y_axis = Vector([0.0, 1.0, 0.0])
        z_axis = Vector([0.0, 0.0, 1.0])

        _r = np.zeros([1, 3])

        # --- Optimize entrance fringe --- #
        deviation_x = 1e20  # Some large number as initialization
        it = 0
        while abs(deviation_x) > tol:

            # Apply the angle correction (first time use initial guess)
            self._set_blim(b_min=0.0 + _db[0])
            # (Re-)calculate the new geometry and BEM++ solution
            self.generate_meshed_model()
            self.solve_bempp()

            # Trajectory starting point:
            _trj = self._variables_analytic["trj_design"]  # type: np.ndarray
            _v_des = self._variables_analytic["v_design"]  # type: np.ndarray

            # Calculate E-Field upwards of starting point
            xs, ys, zs = _trj[0]
            self.calculate_efield(limits=((xs - 0.5 * _hcl, xs + 0.5 * _hcl),
                                          (ys - 0.5 * _hcl, ys + 0.5 * _hcl),
                                          (zs - 1.9 * _hcl, zs + 0.1 * _hcl)),
                                  res=res,
                                  domain_decomp=(4, 4, 4),
                                  overlap=0)

            _r, _v = self.track(r_start=_trj[0],  # Use starting point of design particle
                                v_start=-_v_des[0],  # Reverse direction of design particle
                                nsteps=5000,
                                dt=1e-11,
                                omit_b=True)

            trj_dir = Vector(_r[-1] - _r[-100])
            deviation_x = 90.0 - np.rad2deg(trj_dir.angle_with(x_axis))
            deviation_y = 90.0 - np.rad2deg(trj_dir.angle_with(y_axis))

            print("Current entrance adjustment: {:.4f}, "
                  "deviation from z-axis in x-dir: {:.4f} degrees, "
                  "deviation from z-axis in y-dir: {:.4f} degrees".format(_db[0], deviation_x, deviation_y))

            if it == 1:
                _db[0] += deviation_x
            else:
                _db[0] += 0.5 * deviation_x  # Dampen the oscillations a bit

            it += 1

            if it == maxiter:
                print("Entrance Fringe: Maximum number of iterations has been reached. Breaking.")
                break

        # Save the coordinates of trj endpoint --> need to shift x and y later
        shift = _r[-1]

        # --- Optimize exit fringe --- #
        deviation = 1e20  # Some large number as initialization
        it = 0
        while abs(deviation) > tol:

            # Apply the angle correction (first time use initial guess)
            self._set_blim(b_max=90.0 + _db[1])
            # (Re-)calculate the new geometry and BEM++ solution
            self.generate_meshed_model()
            self.solve_bempp()

            # Trajectory starting point:
            _trj = self._variables_analytic["trj_design"]  # type: np.ndarray
            _v_des = self._variables_analytic["v_design"]  # type: np.ndarray

            _ns = self._params_analytic["ns"]  # type: int
            start_idx = int(0.9 * _ns)  # start the tracking "10 percent" into the spiral inflector exit
            # xs, ys, zs = _trj[start_idx]
            rs = shift
            rs[2] = -0.15
            xs, ys, zs = rs
            vs = _v_des[0]
            # Calculate E-Field
            # TODO: Better way to determine the fringe field region
            self.calculate_efield(limits=((xs - 2.0 * _hcl, xs + 2.0 * _hcl),
                                          (ys - 2.0 * _hcl, ys + 2.0 * _hcl),
                                          (zs - _hcl, zs + _hcl)),
                                  res=res,
                                  domain_decomp=(4, 4, 4),
                                  overlap=0)

            _r, _v = self.track(r_start=rs,  # Use point close to exit along design particle
                                v_start=vs,  # Regular direction of design particle
                                nsteps=12500,
                                dt=1e-11,
                                omit_b=False)

            trj_dir = Vector(_r[-1] - _r[-100])
            deviation = 90.0 - np.rad2deg(trj_dir.angle_with(z_axis))

            print("Current exit adjustment: {:.4f}, "
                  "deviation from xy-plane: {:.4f} degrees".format(_db[1], deviation))

            if it == 1:
                _db[1] += deviation
            else:
                _db[1] += 0.65 * deviation  # Dampen the oscillations a bit

            it += 1

            if it == maxiter:
                print("Exit Fringe: Maximum number of iterations has been reached. Breaking.")
                break

        shift[2] = _r[-1, 2]
        self._variables_track["shift"] = -shift

        print("Done optimizing!")

        return shift

    def plot_potential(self, lims=[-0.1, 0.1, -0.1, 0.1], orientation="xy", **kwargs):

        if self._variables_bempp["solution"] is None:
            print("No BEM++ solution in memory to plot from.")
            return 1

        if "n_grid_points" in kwargs.keys():
            n_grid_points = kwargs["n_grid_points"]
        else:
            n_grid_points = 200

        if "offset" in kwargs.keys():
            offset = kwargs["offset"]
        else:
            offset = 0.0

        sol = self._bem_sol
        f_space = self._space

        plot_grid = np.mgrid[limits[0]:limits[1]:n_grid_points * 1j, \
                             limits[2]:limits[3]:n_grid_points * 1j]

        e1 = plot_grid[0].ravel()
        e2 = plot_grid[1].ravel()
        e3 = offset * np.ones(plot_grid[0].size)
        if orientation == "xy":
            points = np.vstack((e1, e2, e3))
        elif orientation == "yz":
            points = np.vstack((e3, e1, e2))
        elif orientation == "xz":
            points = np.vstack((e1, e3, e2))

        slp_pot = bempp.api.operators.potential.laplace.single_layer(f_space, points)
        u_evaluated = slp_pot * sol
        u_evaluated = u_evaluated.reshape((n_grid_points, n_grid_points))

        if vlims in kwargs.keys():
            vlims = kwargs["vlims"]
        else:
            vlims = [np.min(u_evaluated), np.max(u_evaluated)]

        plt.imshow(u_evaluated.T,
                   extent=(limits[0], limits[1], limits[2], limits[3]),
                   vmin=vlims[0],
                   vmax=vlims[1])

        if "colorbar" in kwargs.keys():
            if kwargs["colorbar"]:
                plt.colorbar()

        plt.show()

        if "save_fig" in kwargs.keys() and "filename" in kwargs.keys():
            if kwargs["save_fig"]:
                plt.savefig(filename)

        return 0

    def save_geo_files(self, filename=None):

        if filename is None:
            fd = FileDialog()
            folder, _ = fd.get_filename("folder")
            if folder is None:
                return 0
        else:
            folder = os.path.split(filename)[0]

        for name, electrode in self._variables_bempp["objects"].items():
            if self._debug:
                print("Working on {}".format(name))

            filename = os.path.join(folder, "{}.geo".format(name))
            with open(filename, 'w') as outfile:
                outfile.write(electrode["gmsh_str"])

        return 0

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

        elif key in self._params_bempp.keys():

            self._params_bempp[key] = value

        elif group in self._params_bempp.keys():

            if key in self._params_bempp[group].keys():

                self._params_bempp[group][key] = value

    def solve_bempp(self):

        if self._variables_bempp["full mesh"] is None:
            print("Please generate a mesh before solving with BEM++!")
            return 1

        print("Generating necessary BEM++ operators, function spaces. Solving... ", end="")

        dp0_space = bempp.api.function_space(self._variables_bempp["full mesh"], "DP", 0)
        slp = bempp.api.operators.boundary.laplace.single_layer(dp0_space, dp0_space, dp0_space)

        domain_mapping = {}
        for name, electrode in self._variables_bempp["objects"].items():
            domain_mapping[electrode["domain"]] = electrode["voltage"]

        def f(*args):

            domain_index = args[2]
            result = args[3]

            result[0] = domain_mapping[domain_index]

        dirichlet_fun = bempp.api.GridFunction(dp0_space, fun=f)

        if self._debug:
            dirichlet_fun.plot()

        # Solve
        sol, info = bempp.api.linalg.gmres(slp, dirichlet_fun, tol=1e-5, use_strong_form=True)

        print("Done!")

        # Save results
        self._variables_bempp["solution"] = sol
        self._variables_bempp["f_space"] = dp0_space
        self._variables_bempp["operator"] = slp
        self._variables_bempp["grid_fun"] = dirichlet_fun

        return 0

    def track(self, r_start=None, v_start=None, nsteps=10000, dt=1e-12, omit_b=False, omit_e=False):

        # TODO: For now break if r_start or v_start are not given, later get from class properties?
        assert(r_start is not None and v_start is not None), "Have to specify r_start and v_start for now!"

        if self._variables_track["ef_itp"] is None:
            print("No E-Field has been generated. Cannot track!")
            return 1

        pusher = ParticlePusher(self._params_analytic["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

        if omit_e:
            efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
        else:
            efield1 = self._variables_track["ef_itp"]  # type: Field

        if omit_b:
            bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
        else:
            bfield1 = self._params_analytic["bf_itp"]  # type: Field

        r = np.zeros([nsteps + 1, 3])
        v = np.zeros([nsteps + 1, 3])
        r[0, :] = r_start[:]
        v[0, :] = v_start[:]

        # initialize the velocity half a step back:
        ef = efield1(r[0])
        bf = bfield1(r[0])
        _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

        # Track for n steps
        for i in range(nsteps):
            ef = efield1(r[i])
            bf = bfield1(r[i])

            r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)

        self._variables_track["trj_tracker"] = r

        return r, v

    def __generate_numerical_trajectory(self, bf=None, nsteps=100000, dt=1e-12):

        pusher = ParticlePusher(self._params_analytic["ion"], "boris")  # Note: leapfrog is inaccurate above dt = 1e-12

        tilt = self._params_analytic["tilt"]  # type: float
        self._variables_analytic["kp"] = np.tan(np.deg2rad(tilt))  # Tilt parameter

        r_start = np.array([0.0, 0.0, -self._variables_analytic["height"]])
        v_start = np.array([0.0, 0.0, self._params_analytic["ion"].v_m_per_s()])

        _r = np.zeros([nsteps + 1, 3])
        _v = np.zeros([nsteps + 1, 3])
        _b = np.zeros([nsteps + 1])  # Store the "b" angle for geometry generation
        _r[0, :] = r_start[:]
        _v[0, :] = v_start[:]

        # Create a new electric field, which will be repeatedly re-defined
        field_val = self._variables_analytic["ef_design"]
        efield1 = Field(dim=0, field={"x": field_val, "y": 0.0, "z": 0.0})

        if bf is not None:
            bfield1 = bf
        else:
            bfield1 = self._params_analytic["bf_itp"]

        # initialize the velocity half a step back:
        ef = efield1(_r[0])
        bf = bfield1(_r[0])
        _, _v[0] = pusher.push(_r[0], _v[0], ef, bf, -0.5 * dt)

        # Track for n steps
        for i in range(nsteps):
            ef = efield1(_r[i])
            bf = bfield1(_r[i])

            _r[i + 1], _v[i + 1] = pusher.push(_r[i], _v[i], ef, bf, dt)

            vx, vy, vz = _v[i + 1]
            vo = np.sqrt(vx**2.0 + vy**2.0 + vz**2.0)
            _b[i + 1] = i * dt * vo / self._variables_analytic["height"]

            # Toprek theory with surgery
            Eh = field_val * self._variables_analytic["kp"] * np.sin(_b[i + 1])
            Ehx = -Eh * vy / (np.sqrt(vo**2.0 - vz**2.0))
            Ehy = Eh * vx / (np.sqrt(vo**2.0 - vz**2.0))

            ex = field_val * vx * np.abs(vz) / (vo * np.sqrt(vo**2.0 - vz**2.0)) + Ehx
            ey = field_val * vy * np.abs(vz) / (vo * np.sqrt(vo**2.0 - vz**2.0)) + Ehy
            ez = -field_val * (vo**2.0 - vz**2.0) / (vo * np.sqrt(vo**2.0 - vz**2.0))
            efield1 = Field(dim=0, field={"x": ex, "y": ey, "z": ez})
            if vz < 0:  # Stop when the z-component of the velocity is zero
                if self._debug:
                    print(_r[i + 1, :])  # Print the final position
                break

        ns = i
        try:
            i_init = np.where(_b >= self._params_analytic["b_lim"][0])[0][0]
            if i_init == 0:
              i_init = 1
        except IndexError:
            i_init = 1

        try:
            i_final = np.where(_b >= self._params_analytic["b_lim"][1])[0][0]
        except IndexError:
            i_final = i

        if self._debug:
            print("Design Trajectory: Initial index: %i, Final index: %i.".format(i_init, i_final))

        # The arrays cannot be as large as they're made initially.
        # This would cause the BEMPP routines to perform the computations
        # with incredibly high resolution and would never finish. -PW
        step = int(np.floor((ns / 50)))
        interval = [j for j in range(i_init, i_final, step)]
        interval.append(i_final)

        b = _b[interval]
        r = _r[interval, :]
        v = _v[interval, :]

        # Redefine the analytical variables (will change name eventually)
        self._variables_analytic["trj_design"] = r
        self._variables_analytic["trj_vel"] = v
        self._variables_analytic["b"] = b
        self._params_analytic["ns"] = len(r[:, 0])

        return r, v

    def __generate_numerical_geometry(self):
        # This is a slightly modified version of the normal analytical method
        if self._variables_analytic["trj_design"] is None:
            print("No numerical design trajectory yet, generating...")
            self.generate_design_trajectory()

        print("Generating numerical geometry... ", end="")

        geos = []
        _ion = self._params_analytic["ion"]  # type: IonSpecies
        ns = self._params_analytic["ns"]  # type: int
        gap = self._params_analytic["gap"]  # type: float
        sigma = self._params_analytic["sigma"]  # type: float
        kp = self._variables_analytic["kp"]  # type: float

        b = self._variables_analytic["b"]  # type: np.ndarray
        trj_design = self._variables_analytic["trj_design"]  # type: np.ndarray
        trj_vel = self._variables_analytic["trj_vel"]
        vx, vy, vz = trj_vel[:, 0], trj_vel[:, 1], trj_vel[:, 2]  # Unload the v components

        for thickness in [0.0, self._params_analytic["dx"]]:

            end_distance = 2.0 * thickness + gap  # End to end distance of electrodes in (m)
            aspect_ratio = 2.5 * (gap / end_distance)

            # Distance between electrodes at inflection angle theta
            d = end_distance * np.ones(ns) / (np.sqrt(1.0 + ((kp ** 2.0) * (np.sin(b)) ** 2.0)))

            # Save inner gap size vs deflection angle as class variable
            if thickness == 0.0:
                self._variables_analytic["d"] = d

            # Rotation/flip
            if not ((self._variables_analytic["bf_design"] > 0.0) ^ (_ion.q() > 0.0)):
                if self._debug:
                    print("Flipping direction of cyclotron motion...", end="")
                vy = -vy

            v2 = np.sqrt((vx ** 2.0) + (vy ** 2.0))  # xy-magnitude of the velocity
            v3 = np.sqrt((vx ** 2.0) + (vy ** 2.0) + (vz ** 2.0))  # 3-d magnitude of the velocity vector (Should = v)

            # Save vx, vy, vz in as class variable in same format as trj_design
            self._variables_analytic["v_design"] = np.array([vx, vy, vz]).T

            # Construction of the normalized vectors of the optical coordinate system
            v_path = np.transpose(np.array([vx, vy, vz]))
            v_optical = np.zeros((3, ns, 3))  # [h, u, v]

            for j in range(ns):

                for k in range(3):
                    v_optical[2, j, k] = v_path[j, k] / v3[j]

                if v2[j] == 0:
                    v_optical[0, j, :] = np.array([0, -1, 0])
                    v_optical[1, j, :] = np.array([-1, 0, 0])

                elif v2[j] > 0:
                    v_optical[0, j, :] = np.array([vy[j] / v2[j], -vx[j] / v2[j], 0])
                    v_optical[1, j, :] = np.array(
                        [-(vx[j] * vz[j]) / (v3[j] * v2[j]), -(vy[j] * vz[j]) / (v3[j] * v2[j]),
                         ((v3[j] ** 2) - (vz[j] ** 2)) / (v3[j] * v2[j])])

            # Rotation to the tilt angle if atilt > 0 || atilt < 0
            v_rh = np.copy(v_optical)  # "Rotated" or "right-handed" optical coordinate system vectors [hr, ur, v]

            # Note: the theta in the commented eq. below is not the same theta as in the code
            nemo = np.arctan(kp * np.sin(b))

            for i in range(ns):
                for j in range(3):
                    v_rh[0, i, j] = ((np.cos(nemo[i])) * v_optical[0, i, j]) - (np.sin(nemo[i])) * v_optical[1, i, j]
                    v_rh[1, i, j] = ((np.cos(nemo[i])) * v_optical[1, i, j]) + (np.sin(nemo[i])) * v_optical[0, i, j]

            # Turn track of unit vectors
            t1 = np.arange(0, 5, 0.01)
            v_er = np.zeros((3, 3, np.size(t1)))

            for i in range(3):
                for j in range(3):
                    for k in range(np.size(t1)):
                        v_er[i, j, k] = trj_design[ns - 1, j] + t1[k] * v_rh[i, ns - 1, j]

            # Construction of the electrodes
            edge_lines = np.zeros((5, ns, 3))

            xi = 0.5 * aspect_ratio * end_distance

            for i in range(ns):
                for j in range(3):
                    edge_lines[0, i, j] = trj_design[i, j] + 0.5 * d[i] * v_rh[1, i, j] + xi * v_rh[0, i, j]
                    edge_lines[1, i, j] = trj_design[i, j] + 0.5 * d[i] * v_rh[1, i, j] - xi * v_rh[0, i, j]
                    edge_lines[2, i, j] = trj_design[i, j] - 0.5 * d[i] * v_rh[1, i, j] + xi * v_rh[0, i, j]
                    edge_lines[3, i, j] = trj_design[i, j] - 0.5 * d[i] * v_rh[1, i, j] - xi * v_rh[0, i, j]

            geos.append(edge_lines)

        # Apply the v-shape modification:
        diff = (geos[0][0, :, :] - geos[1][0, :, :])  # Difference between ext and int
        diff_norm = (np.sqrt(diff[:, 0] ** 2 + diff[:, 1] ** 2 + diff[:, 2] ** 2))  # Magnitude of diff
        diff_hat = np.zeros(np.shape(diff))  # Initialize normalized vector

        for i in range(ns):
            diff_hat[i, :] = diff[i, :] / diff_norm[i]  # Calculate each normalized vector

        geo = np.zeros([10, ns, 3])  # Initialize geo array

        # Upper Electrode
        geo[0, :, :] = geos[0][0, :, :] + diff_hat * sigma
        geo[1, :, :] = geos[0][1, :, :] + diff_hat * sigma
        geo[2, :, :] = geos[1][0, :, :] + diff_hat * sigma  # Should reduce the thickness
        geo[3, :, :] = geos[1][1, :, :] + diff_hat * sigma  # Should reduce the thickness
        geo[4, :, :] = 0.5 * (geos[0][0, :, :] + geos[0][1, :, :])

        # Lower Electrode
        geo[5, :, :] = geos[0][2, :, :] + diff_hat * sigma
        geo[6, :, :] = geos[0][3, :, :] + diff_hat * sigma
        geo[7, :, :] = geos[1][2, :, :]
        geo[8, :, :] = geos[1][3, :, :]
        geo[9, :, :] = 0.5 * (geos[0][2, :, :] + geos[0][3, :, :])

        self._variables_analytic["geo"] = geo

        print("Done!")

        return self._variables_analytic["geo"]

    def plot_trajectories(self, trajectories=None):
        # Quick method for plotting trajectories
        fig = plt.figure()
        ax = Axes3D(fig)

        proj3d.persp_transformation = orthogonal_proj

        # Plot the beam trajectory
        for trj in trajectories:
          ax.plot(trj[:, 0], trj[:, 1], trj[:, 2])

        plt.show()

    def plot_bfield(self, lims=[-0.2, 0.2], num=5000):
        zpts = np.linspace(lims[0], lims[1], num)
        b = np.zeros(num)
        plt.rc('text', usetex=True)
        plt.rc('font', family="serif")
        matplotlib.rcParams.update({'font.size': 16})
        fig, ax = plt.subplots()
        for i in range(num):
          zpt = zpts[i]
          b[i] = self._params_analytic["bf_itp"](np.array([0.0, 0.0, zpt]))[2]
        ax.plot(np.array(zpts), -np.array(b), colors[0], linewidth=3.0)
        plt.xlim([-0.2, 0.2])
        plt.xlabel('z (m)')
        plt.ylabel('B (T)')
        plt.grid(True)
        plt.savefig("bfield.png")

if __name__ == "__main__":
    h2p = IonSpecies("H2_1+", 0.035)
    h2p.calculate_from_energy_mev(0.07 / h2p.a())

    si = SpiralInflector(ion=h2p,
                         method="numerical",
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
    si.draw_geometry(freq=10, show=True)

    si.set_parameter(key="h", value=0.005)  # Mesh characteristic length
    si.set_parameter(key="make_aperture", value=True)
    si.set_parameter(key="aperture_params", value={"thickness": 4e-3,
                                                   "radius": 50e-3,
                                                   "length": 45e-3,
                                                   "width": 18e-3,
                                                   "distance": 10e-3,
                                                   "voltage": 0.0})
    si.set_parameter(key="make_cylinder", value=False)
    si.set_parameter(key="cylinder_params", value={"radius": 120e-3,
                                                   "zmin": -150e-3,
                                                   "zmax": 80e-3,
                                                   "voltage": 0.0})

    si.generate_meshed_model()

    ts = time.time()
    myshift = si.optimize_fringe(initial_guess=(2.8683, -4.9849), maxiter=5, tol=0.01, res=0.002)
    print("Optimizing took {:.4f} s".format(time.time() - ts))

    rstart = myshift
    rstart[2] = -0.15

    ts = time.time()
    print("Calculating electric field...")
    si.calculate_efield(res=0.002,
                        limits=((-0.08, 0.08), (-0.08, 0.08), (-0.12, 0.05)),
                        domain_decomp=(7, 7, 7))
    print("Calculating field took {:.4f} s".format(time.time() - ts))

    with open('timing.txt', 'a') as outfile:
        outfile.write("Generating electric field took {:.4f} s\n".format(time.time() - ts))

    ts = time.time()
    si.track(r_start=rstart,
             v_start=np.array([0.0, 0.0, h2p.v_m_per_s()]),
             nsteps=15000,
             dt=1e-11)
    print("Tracking took {:.4f} s".format(time.time() - ts))

    with open('timing.txt', 'a') as outfile:
        outfile.write("Tracking took {:.4f} s\n".format(time.time() - ts))

    si.draw_geometry(show=True, filename='auto')
    si.export_electrode_geometry(fname='electrode_macro.ivb')
    si.export_aperture_geometry(fname='aperture_macro.ivb')
    si.save_geo_files()
