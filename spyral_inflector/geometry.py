from dans_pymodules import *
import bempp.api
# noinspection PyUnresolvedReferences
from bempp.api.shapes.shapes import __generate_grid_from_geo_string as generate_from_string
# noinspection PyUnresolvedReferences
from py_electrodes.py_electrodes import PyElectrode, PyElectrodeAssembly

X_AXIS = np.array([1, 0, 0], float)
Y_AXIS = np.array([0, 1, 0], float)
Z_AXIS = np.array([0, 0, 1], float)


class SIAperture(PyElectrode):
    def __init__(self, parent=None, name="New Aperture", voltage=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture

    def create_geo_str(self, r, dz, a, b, hole_type="ellipse", h=0.005, load=True):
        """

        Creates the geo string for a circular aperture plate with a elliptical or rectangular hole
        For circular or square holes set a = b
        This plate is centered around the origin (local coordinate system) with surface normal in z direction
        and needs to be shifted/rotated.

        :param r: plate radius
        :param dz: plate thickness
        :param a: ellipse/rectangle long half-axis
        :param b: ellipse/rectangle short half-axis
        :param hole_type: "ellipse", "rectangle"
        :param h: desired mesh resolution
        :param load: Flag whether to also load from geo string directly.
                     Cave: If False, geo str will not be saved internally!
        :return gmsh_str: the string object for gmsh
        """

        geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)

        geo_str += "// Base plate\n"
        geo_str += "Cylinder(1) = {{ 0, 0, {}, 0, 0, {}, {}, 2 * Pi }};\n\n".format(-0.5 * dz, dz, r)

        geo_str += "// Tool to subtract\n"
        if hole_type == "rectangle":
            geo_str += "Box(2) = {{ {}, {}, {}, {}, {}, {} }};\n\n".format(-0.5 * a, -0.5 * b, -dz, a, b, 2 * dz)
        elif hole_type == "ellipse":
            geo_str += "Disk(100) = {{ 0, 0, {}, {}, {} }};\n".format(-dz, 0.5 * a, 0.5 * b)
            geo_str += "Extrude {{ 0, 0, {} }} {{ Surface{{ 100 }}; }}\n".format(2 * dz)
        else:
            print("Don't understand hole type {}!".format(hole_type))
            return 1

        geo_str += "\nBooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };\n"

        # Call function in PyElectrode module we inherit from if 'load' is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


class SIPointSphere(PyElectrode):
    def __init__(self, parent=None, name="New Aperture", voltage=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture

    def create_geo_str(self, center, r=0.001, h=0.005, load=True):
        """

        Creates the geo string for a small sphere to show important vertices

        :param center: numpy array, list or tuple containing 3 floats forthe center of the sphere
        :param r: sphere radius - default is 1 mm
        :param h: desired mesh resolution
        :param load: Flag whether to also load from geo string directly.
                     Cave: If False, geo str will not be saved internally!
        :return gmsh_str: the string object for gmsh
        """

        center = np.asarray(center)
        assert center.shape == (3,), "Got wrong dimension of {} for center, should be (3, )".format(center.shape)

        geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)

        geo_str += "// Base plate\n"
        geo_str += "Sphere(1) = {{ {}, {}, {}, {} }};\n".format(center[0], center[1], center[2], r)

        # Call function in PyElectrode module we inherit from if 'load' is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


class SICylinder(PyElectrode):
    def __init__(self, parent=None, name="New Cylinder", voltage=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture

    def create_geo_str(self, r, dz, h=0.005, load=True):
        """

        Creates the geo string for a circular aperture plate with a elliptical or rectangular hole
        For circular or square holes set a = b
        This plate is centered around the origin (local coordinate system) with surface normal in z direction
        and needs to be shifted/rotated.

        :param r: cylinder radius
        :param dz: height
        :param h: desired mesh resolution
        :param load: Flag whether to also load from geo string directly.
                     Cave: If False, geo str will not be saved internally!
        :return gmsh_str: the string object for gmsh
        """

        geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)

        geo_str += "// Cylinder\n"
        geo_str += "Cylinder(1) = {{ 0, 0, {}, 0, 0, {}, {}, 2 * Pi }};\n\n".format(-0.5 * dz, dz, r)

        # Call function in PyElectrode module we inherit from if load is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


class SIElectrode(PyElectrode):
    def __init__(self, parent=None, name="New Spiral Electrode", voltage=10000):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture

    def create_geo_str(self, raw_geo, elec_type, h=0.005, load=True):
        """

        Creates the geo string for a circular aperture plate with a elliptical or rectangular hole
        For circular or square holes set a = b
        This plate is centered around the origin (local coordinate system) with surface normal in z direction
        and needs to be shifted/rotated.

        :param raw_geo: ndim=3 numpy array containing the guide rails of the spiral electrodes
        :param elec_type: 'anode' or 'cathode'
        :param h: desired mesh resolution
        :param load: Flag whether to also load from geo string directly.
                     Cave: If False, geo str will not be saved internally!
        :return gmsh_str: the string object for gmsh
        """

        if elec_type not in ["anode", "cathode"]:
            print("SIElectrode could not understand electrode type {}. Must be 'Anode' or 'Cathode'".format(type))
            return 1

        geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)

        new_pt = 1
        new_ln = 1
        new_loop = 1
        new_vol = 1

        # Shift the index in geo object for anode or cathode...
        if elec_type == "anode":
            k = 0
        elif elec_type == "cathode":
            k = 5

        num_sections = len(raw_geo[0, :, 0])

        for j in range(num_sections):
            for i in range(5):
                geo_str += "Point({}) = {{ {}, {}, {} }};\n".format(new_pt,
                                                                    raw_geo[i + k, j, 0],
                                                                    raw_geo[i + k, j, 1],
                                                                    raw_geo[i + k, j, 2])
                new_pt += 1

            # For each section, add the lines
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 0, (j * 5) + 4, (j * 5) + 2)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 1, (j * 5) + 3, (j * 5) + 1)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 2, (j * 5) + 3, (j * 5) + 4)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 3, (j * 5) + 5, (j * 5) + 2)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 4, (j * 5) + 1, (j * 5) + 5)

            new_ln += 5

            geo_str += "Wire({}) = {{ {}, {}, {}, {}, {} }};\n\n".format(new_loop,
                                                                         (j * 5) + 3,
                                                                         (j * 5) + 2,
                                                                         (j * 5) + 5,
                                                                         (j * 5) + 4,
                                                                         (j * 5) + 1)

            new_loop += 1

        geo_str += "Ruled ThruSections({}) = {{ 1:{} }};".format(new_vol, new_loop - 1)

        new_vol += 1

        # Call function in PyElectrode module we inherit from if load is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


# Geometrically, trajectories have much in common with electrodes...
class SITrajectory(PyElectrode):

    def __init__(self, parent=None, name="New Spiral Electrode", voltage=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture

    # --- Override some of the PyElectrode functions that don't make sense for a wire --- #
    @staticmethod
    def points_inside(_points):
        # By definition no points are 'inside' a PyWire
        return np.zeros(_points.shape, bool)

    @staticmethod
    def generate_mesh(brep_h=0.0):
        print("Can't generate a mesh from a PyWire")
        return 1

    def create_geo_str(self, points, max_points, load=True):
        """
        Create a geo string for gmsh
        :param points: np array of points along the trajectory
        :param max_points: maximum number of points to use if max is larger than num points all are used
        :param load: immediately load the geo str as an occ object
        :return:
        """

        points = np.asarray(points)

        assert points.ndim == 2 and points[0, :].shape == (3,), "points have wrong shape = {}".format(points.shape)

        # Reduce number of points to use in spline to max_points
        points = points[::int(np.ceil(len(points)) / max_points), :]

        geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
"""
        new_pt = 1
        new_ln = 1

        geo_str += "// Center Spline:\n"
        for _x, _y, _z in points:
            geo_str += "Point({}) = {{ {}, {}, {} }};\n".format(new_pt, _x, _y, _z)
            new_pt += 1

        geo_str += """
Spline({}) = {{ {}:{} }}; 
""".format(new_ln, 1, new_pt - 1)

        # Immediately delete the points used up in the spline
        geo_str += "Recursive Delete {{ Point{{ {}:{} }}; }}\n".format(1, new_pt - 1)

        # Call function in PyElectrode module we inherit from if load is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


class SIHousing(PyElectrode):

    def __init__(self, parent=None, name="Spiral Inflector Housing", voltage=0, experimental=False):
        super().__init__(name=name, voltage=voltage)

        assert parent is not None, "This class requires a parent."

        self._parent = parent
        self._aperture_params = None

        self._experimental = experimental

    def set_aperture_params(self, parameters):
        self._aperture_params = parameters

    def set_aperture_rot_angles(self, angles):
        self._tilt_angle, self._face_angle = angles

    def gen_convex_hull(self, geo, gap, thickness):
        from scipy.spatial import ConvexHull

        geo_list = []
        for i in range(9):
            geo_list.append(geo[i, :, :2])
        points = np.concatenate(geo_list)

        hull = ConvexHull(points)
        hull_pts = points[hull.vertices, :]

        hull_pts = np.concatenate((hull_pts, hull_pts[0, :][np.newaxis, :]), axis=0)

        if self._debug:
            plt.plot(hull_pts[:, 0], hull_pts[:, 1][:, np.newaxis], 'r--')

        # Idea: use the convex hull, generate points in a circle around each point of the hull and perform
        # another convex hull on that new set of points.
        def method_two():
            total_pts = []
            circle_pts = []
            circle_res = 48
            for i in range(circle_res):
                circle_pts.append(gap * np.array([np.cos(i * 2 * np.pi / circle_res),
                                                  np.sin(i * 2 * np.pi / circle_res)]))

            for pt in hull_pts:
                for cpt in circle_pts:
                    total_pts.append(pt + cpt)

            total_pts = np.array(total_pts)

            new_hull_inner = ConvexHull(total_pts)
            new_hull_pts_inner = total_pts[new_hull_inner.vertices, :]
            if self._debug:
                plt.plot(new_hull_pts_inner[:, 0], new_hull_pts_inner[:, 1], 'g')

            total_pts = []
            circle_pts = []
            circle_res = 48
            for i in range(circle_res):
                circle_pts.append((gap + thickness) * np.array([np.cos(i * 2 * np.pi / circle_res),
                                                                np.sin(i * 2 * np.pi / circle_res)]))

            for pt in hull_pts:
                for cpt in circle_pts:
                    total_pts.append(pt + cpt)

            total_pts = np.array(total_pts)

            new_hull_outer = ConvexHull(total_pts)
            new_hull_pts_outer = total_pts[new_hull_outer.vertices, :]

            if self._experimental:
                tilt_angle, face_angle = get_norm_vec_and_angles_from_geo(geo)
                face_vector = np.array([np.cos(face_angle), np.sin(face_angle), 0.0])
                norm_vector = np.cross(face_vector, np.array([0.0, 0.0, 1.0]))

                new_point_a_in = geo[2, -1, :2] + norm_vector[:2] * gap
                new_point_b_in = geo[8, -1, :2] + norm_vector[:2] * gap

                new_point_a_out = geo[2, -1, :2] + norm_vector[:2] * (gap + thickness)
                new_point_b_out = geo[8, -1, :2] + norm_vector[:2] * (gap + thickness)

                plt.scatter(new_point_a_in[0], new_point_a_in[1], color='r')
                plt.scatter(new_point_b_in[0], new_point_b_in[1], color='r')

                new_hull_pts_inner = np.vstack([new_hull_pts_inner, new_point_a_in])
                new_hull_pts_inner = np.vstack([new_hull_pts_inner, new_point_b_in])

                new_hull_pts_outer = np.vstack([new_hull_pts_outer, new_point_a_out])
                new_hull_pts_outer = np.vstack([new_hull_pts_outer, new_point_b_out])

            if self._debug:
                plt.plot(new_hull_pts_outer[:, 0], new_hull_pts_outer[:, 1], 'g')

            return new_hull_pts_inner, new_hull_pts_outer

        pts_in, pts_out = method_two()

        if self._debug:
            plt.show()

        return pts_in, pts_out

    def create_geo_str(self, geo, trj, zmin, zmax, gap, thickness, h=0.005, load=True):

        pts_in, pts_out = self.gen_convex_hull(geo, gap, thickness)

        dz = self._aperture_params["thickness"]
        r = self._aperture_params["radius"]
        a = self._aperture_params["length"]
        b = self._aperture_params["width"]
        t_gap = self._aperture_params["top_distance"]
        b_gap = self._aperture_params["bottom_distance"]

        norm_vec = Vector(trj[-1] - trj[-2]).normalized()

        translate = np.array([trj[-1][0] + norm_vec[0] * b_gap,
                              trj[-1][1] + norm_vec[1] * b_gap,
                              0.0])

        hole_type = self._aperture_params["hole_type"]

        # TODO: ToleranceBoolean is important for the aperture hole. This was tested from 1T to about 2.4T and
        # TODO: seems to work... 3T did not. -PW
        geo_str = """SetFactory("OpenCASCADE");
// Geometry.NumSubEdges = 100; // nicer display of curve
Geometry.ToleranceBoolean = 1E-5;
Geometry.Tolerance = 1E-10;
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
                """.format(h)

        geo_str += "// Outside points\n"
        n_pts_out = np.shape(pts_out)[0]

        for i, pt in enumerate(pts_out):
            geo_str += "Point({}) = {{ {}, {}, 0}};\n".format(i, pt[0], pt[1])

        geo_str += "// Outside lines\n"
        n_lines_out = n_pts_out  # Should be the same number

        for i in range(n_lines_out - 1):
            geo_str += "Line({}) = {{ {}, {} }};\n".format(i, i, i + 1)
        geo_str += "Line({}) = {{ {}, {} }};\n".format(i + 1, n_pts_out - 1, 0)  # Connect last point to first point

        geo_str += "// Inside points\n"
        n_pts_in = np.shape(pts_in)[0]

        for i, pt in enumerate(pts_in):
            geo_str += "Point({}) = {{ {}, {}, 0}};\n".format(i + n_pts_out, pt[0], pt[1])

        geo_str += "// Inside lines\n"
        n_lines_in = n_pts_in  # Should be the same number

        for i in range(n_lines_in - 1):
            geo_str += "Line({}) = {{ {}, {} }};\n".format(i + n_lines_out, i + n_pts_out, i + 1 + n_pts_out)
        # Connect last point to first point
        geo_str += "Line({}) = {{ {}, {} }};\n".format(i + 1 + n_lines_out, i + 1 + n_pts_out, n_pts_out)

        geo_str += "Line Loop(1) = {{ {}:{} }};\n".format(0, n_lines_out - 1)
        geo_str += "Line Loop(2) = {{ {}:{} }};\n".format(n_lines_out, n_lines_out + n_lines_in - 1)

        geo_str += "Plane Surface(1) = {1, 2};\n"
        geo_str += "// Plane Surface(2) = {2};\n"

        geo_str += "Extrude {{ 0, 0, {} }} {{ Surface{{ 1 }}; }}\n".format(zmax - zmin)
        geo_str += "// Extrude {{ 0, 0, {} }} {{ Surface{{ 2 }}; }}\n".format(zmax - zmin)

        # geo_str = 'SetFactory("OpenCASCADE");\n'

        geo_str += "// Tool to subtract\n"
        if hole_type == "rectangle":
            if self._experimental:
                geo_str += "Box(2) = {{ {}, {}, {}, {}, {}, {} }};\n\n".format(-0.5 * a, -0.5 * b,
                                                                               -0.025, a,
                                                                               b, 0.05)
            else:
                geo_str += "Box(2) = {{ {}, {}, {}, {}, {}, {} }};\n\n".format(-0.5 * a, -0.5 * b,
                                                                               0.0, a,
                                                                               b, 0.1)
        elif hole_type == "ellipse":
            geo_str += "Disk(5000) = {{ 0, 0, 0, {}, {} }};\n".format(0.5 * a + 5E-9, 0.5 * b + 5E-9)
            geo_str += "Extrude {{ 0, 0, {} }} {{ Surface{{ 5000 }}; }}\n".format(0.1)  # To be robust

        geo_str += "Rotate {{ {{ 1.0, 0.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ 2 }}; }}\n".format(
            np.pi / 2.0)
        geo_str += "Rotate {{ {{ 0.0, 1.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ 2 }}; }}\n".format(
            self._tilt_angle)
        geo_str += "Rotate {{ {{ 0.0, 0.0, 1.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ 2 }}; }}\n".format(
            self._face_angle)

        geo_str += "Translate {{ {}, {}, {} }} {{ Volume{{ 2 }}; }}\n".format(translate[0],
                                                                              translate[1],
                                                                              -zmin)

        geo_str += "BooleanDifference(50) = { Volume{1}; Delete; }{ Volume{2}; Delete; };\n"
        geo_str += "// BooleanDifference(75) = { Volume{50}; Delete; }{ Volume{3}; Delete; };\n"

        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


# def generate_aperture_geometry(si, electrode_type):
#
#     # Check if all parameters are set
#     abort_flag = False
#     for key, item in si._params_bempp["aperture_params"].items():
#         if item is None:
#             print("Item 'aperture_params/{}' is not set in BEM++ parameters!".format(key))
#             abort_flag = True
#     if abort_flag:
#         return 1
#
#     h = si._params_bempp["h"]
#     geo_str = """SetFactory("OpenCASCADE");
#     Geometry.NumSubEdges = 100; // nicer display of curve
#     Mesh.CharacteristicLengthMax = {};
#             """.format(h * 1.5)
#
#     geo_str += "// Create Plate \n"
#
#     geo = si._variables_analytic["geo"]
#     trj = si._variables_analytic["trj_design"]  # type: np.ndarray
#     thickness = si._params_bempp["aperture_params"]["thickness"]
#     radius = si._params_bempp["aperture_params"]["radius"]
#     length = si._params_bempp["aperture_params"]["length"]
#     width = si._params_bempp["aperture_params"]["width"]
#     top_gap = si._params_bempp["aperture_params"]["top_distance"]
#     bottom_gap = si._params_bempp["aperture_params"]["bottom_distance"]
#     voltage = si._params_bempp["aperture_params"]["voltage"]
#
#     if electrode_type == "top_aperture":
#         zmin = trj[0, 2] + top_gap
#         geo_str += "Cylinder(1) = {{ 0, 0, {}, 0, 0, {}, {}, 2 * Pi }};\n".format(zmin,
#                                                                                   thickness,
#                                                                                   radius)
#     elif electrode_type == "bottom_aperture":
#         i = -1
#     else:
#         print("You must specify either 'top_aperture' or 'bottom_aperture'!")
#         return 1
#
#     geo_str += "BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }\n"
#
#     geo_str += """
# s() = Surface "*";
# """.format(self._domain_idx)
#
#     if reverse_normals:
#         geo_str += """
# ReverseMesh Surface { s() };
# """
#
#     gmsh_str = """SetFactory("OpenCASCADE");
# Geometry.NumSubEdges = 100; // nicer display of curve
# h = {};
# """.format(h * 1.5)
#
#     # Note: These mid-vectors are swapped compared to the ones in the VB macro!
#     mid_vec_a = (geo[4, i, :] - geo[9, i, :])
#     mid_vec_b = (geo[8, i, :] - geo[7, i, :])
#
#     mid_vec_a /= np.linalg.norm(mid_vec_a)
#     mid_vec_b /= np.linalg.norm(mid_vec_b)
#
#     norm_vec = np.cross(mid_vec_b, mid_vec_a)
#     norm_vec /= np.linalg.norm(norm_vec)
#
#     if i == 0:
#         offset_a = -norm_vec * aperture_distance_top
#         offset_b = -norm_vec * (aperture_distance_top + thickness)
#     else:
#         offset_a = norm_vec * aperture_distance_bottom
#         offset_b = norm_vec * (aperture_distance_bottom + thickness)
#
#     gmsh_str += """
# Point(1) = {{ {}, {}, {}, h }};
# Point(2) = {{ {}, {}, {}, h }};
# Point(3) = {{ {}, {}, {}, h }};
# Point(4) = {{ {}, {}, {}, h }};
# Point(5) = {{ {}, {}, {}, h }};""".format(
#         trj[i, 0] + offset_a[0],
#         trj[i, 1] + offset_a[1],
#         trj[i, 2] + offset_a[2],
#         trj[i, 0] + offset_a[0] + mid_vec_a[0] * radius,
#         trj[i, 1] + offset_a[1] + mid_vec_a[1] * radius,
#         trj[i, 2] + offset_a[2] + mid_vec_a[2] * radius,
#         trj[i, 0] + offset_a[0] + mid_vec_b[0] * radius,
#         trj[i, 1] + offset_a[1] + mid_vec_b[1] * radius,
#         trj[i, 2] + offset_a[2] + mid_vec_b[2] * radius,
#         trj[i, 0] + offset_a[0] - mid_vec_a[0] * radius,
#         trj[i, 1] + offset_a[1] - mid_vec_a[1] * radius,
#         trj[i, 2] + offset_a[2] - mid_vec_a[2] * radius,
#         trj[i, 0] + offset_a[0] - mid_vec_b[0] * radius,
#         trj[i, 1] + offset_a[1] - mid_vec_b[1] * radius,
#         trj[i, 2] + offset_a[2] - mid_vec_b[2] * radius)
#
#     gmsh_str += """
# Point(6) = {{ {}, {}, {}, h }};
# Point(7) = {{ {}, {}, {}, h }};
# Point(8) = {{ {}, {}, {}, h }};
# Point(9) = {{ {}, {}, {}, h }};
# Point(10) = {{ {}, {}, {}, h }};""".format(
#         trj[i, 0] + offset_b[0],
#         trj[i, 1] + offset_b[1],
#         trj[i, 2] + offset_b[2],
#         trj[i, 0] + offset_b[0] + mid_vec_a[0] * radius,
#         trj[i, 1] + offset_b[1] + mid_vec_a[1] * radius,
#         trj[i, 2] + offset_b[2] + mid_vec_a[2] * radius,
#         trj[i, 0] + offset_b[0] + mid_vec_b[0] * radius,
#         trj[i, 1] + offset_b[1] + mid_vec_b[1] * radius,
#         trj[i, 2] + offset_b[2] + mid_vec_b[2] * radius,
#         trj[i, 0] + offset_b[0] - mid_vec_a[0] * radius,
#         trj[i, 1] + offset_b[1] - mid_vec_a[1] * radius,
#         trj[i, 2] + offset_b[2] - mid_vec_a[2] * radius,
#         trj[i, 0] + offset_b[0] - mid_vec_b[0] * radius,
#         trj[i, 1] + offset_b[1] - mid_vec_b[1] * radius,
#         trj[i, 2] + offset_b[2] - mid_vec_b[2] * radius)
#
#     gmsh_str += """
# Point(11) = {{ {}, {}, {}, h }};
# Point(12) = {{ {}, {}, {}, h }};
# Point(13) = {{ {}, {}, {}, h }};
# Point(14) = {{ {}, {}, {}, h }};""".format(
#         trj[i, 0] + offset_a[0] + mid_vec_a[0] * width / 2.0 + mid_vec_b[0] * length / 2.0,
#         trj[i, 1] + offset_a[1] + mid_vec_a[1] * width / 2.0 + mid_vec_b[1] * length / 2.0,
#         trj[i, 2] + offset_a[2] + mid_vec_a[2] * width / 2.0 + mid_vec_b[2] * length / 2.0,
#         trj[i, 0] + offset_a[0] - mid_vec_a[0] * width / 2.0 + mid_vec_b[0] * length / 2.0,
#         trj[i, 1] + offset_a[1] - mid_vec_a[1] * width / 2.0 + mid_vec_b[1] * length / 2.0,
#         trj[i, 2] + offset_a[2] - mid_vec_a[2] * width / 2.0 + mid_vec_b[2] * length / 2.0,
#         trj[i, 0] + offset_a[0] - mid_vec_a[0] * width / 2.0 - mid_vec_b[0] * length / 2.0,
#         trj[i, 1] + offset_a[1] - mid_vec_a[1] * width / 2.0 - mid_vec_b[1] * length / 2.0,
#         trj[i, 2] + offset_a[2] - mid_vec_a[2] * width / 2.0 - mid_vec_b[2] * length / 2.0,
#         trj[i, 0] + offset_a[0] + mid_vec_a[0] * width / 2.0 - mid_vec_b[0] * length / 2.0,
#         trj[i, 1] + offset_a[1] + mid_vec_a[1] * width / 2.0 - mid_vec_b[1] * length / 2.0,
#         trj[i, 2] + offset_a[2] + mid_vec_a[2] * width / 2.0 - mid_vec_b[2] * length / 2.0)
#
#     gmsh_str += """
# Point(15) = {{ {}, {}, {}, h }};
# Point(16) = {{ {}, {}, {}, h }};
# Point(17) = {{ {}, {}, {}, h }};
# Point(18) = {{ {}, {}, {}, h }};""".format(
#         trj[i, 0] + offset_b[0] + mid_vec_a[0] * width / 2.0 + mid_vec_b[0] * length / 2.0,
#         trj[i, 1] + offset_b[1] + mid_vec_a[1] * width / 2.0 + mid_vec_b[1] * length / 2.0,
#         trj[i, 2] + offset_b[2] + mid_vec_a[2] * width / 2.0 + mid_vec_b[2] * length / 2.0,
#         trj[i, 0] + offset_b[0] - mid_vec_a[0] * width / 2.0 + mid_vec_b[0] * length / 2.0,
#         trj[i, 1] + offset_b[1] - mid_vec_a[1] * width / 2.0 + mid_vec_b[1] * length / 2.0,
#         trj[i, 2] + offset_b[2] - mid_vec_a[2] * width / 2.0 + mid_vec_b[2] * length / 2.0,
#         trj[i, 0] + offset_b[0] - mid_vec_a[0] * width / 2.0 - mid_vec_b[0] * length / 2.0,
#         trj[i, 1] + offset_b[1] - mid_vec_a[1] * width / 2.0 - mid_vec_b[1] * length / 2.0,
#         trj[i, 2] + offset_b[2] - mid_vec_a[2] * width / 2.0 - mid_vec_b[2] * length / 2.0,
#         trj[i, 0] + offset_b[0] + mid_vec_a[0] * width / 2.0 - mid_vec_b[0] * length / 2.0,
#         trj[i, 1] + offset_b[1] + mid_vec_a[1] * width / 2.0 - mid_vec_b[1] * length / 2.0,
#         trj[i, 2] + offset_b[2] + mid_vec_a[2] * width / 2.0 - mid_vec_b[2] * length / 2.0)
#
#     gmsh_str += """
# Circle(1) = { 2, 1, 3 };
# Circle(2) = { 3, 1, 4 };
# Circle(3) = { 4, 1, 5 };
# Circle(4) = { 5, 1, 2 };
#
# Circle(5) = { 7, 6, 8 };
# Circle(6) = { 8, 6, 9 };
# Circle(7) = { 9, 6, 10 };
# Circle(8) = { 10, 6, 7 };
#
# // Skip Lines 9, 10, 11 because I messed up -PW
#
# Line(11) = { 11, 12 };
# Line(12) = { 12, 13 };
# Line(13) = { 13, 14 };
# Line(14) = { 14, 11 };
#
# Line(15) = { 15, 16 };
# Line(16) = { 16, 17 };
# Line(17) = { 17, 18 };
# Line(18) = { 18, 15 };
#
# Line(19) = { 11, 15 };
# Line(20) = { 12, 16 };
# Line(21) = { 13, 17 };
# Line(22) = { 14, 18 };
#
# Line(23) = { 2, 7 };
# Line(24) = { 3, 8 };
# Line(25) = { 4, 9 };
# Line(26) = { 5, 10 };
#
# Line Loop(1) = { 1, 2, 3, 4 };
# Line Loop(2) = { 5, 6, 7, 8 };
#
# Line Loop(3) = { 11, 12, 13, 14 };
# Line Loop(4) = { 15, 16, 17, 18 };
#
# Line Loop(5) = { 1, 24, -5, -23 };
# Line Loop(6) = { 2, 25, -6, -24 };
# Line Loop(7) = { 3, 26, -7, -25 };
# Line Loop(8) = { 4, 23, -8, -26 };
#
# Line Loop(9) = { 11, 20, -15, -19 };
# Line Loop(10) = { 12, 21, -16, -20 };
# Line Loop(11) = { 13, 22, -17, -21 };
# Line Loop(12) = { 14, 19, -18, -22 };
#
# Plane Surface(1) = { 1, 3 };
# Plane Surface(2) = { 2, 4 };
#
# Ruled Surface(3) = { 5 };
# Ruled Surface(4) = { 6 };
# Ruled Surface(5) = { 7 };
# Ruled Surface(6) = { 8 };
#
# Plane Surface(7) = { 9 };
# Plane Surface(8) = { 10 };
# Plane Surface(9) = { 11 };
# Plane Surface(10) = { 12 };
# """
#
#     si._variables_bempp["objects"][electrode_type] = {"gmsh_str": gmsh_str, "voltage": voltage}
#
#     if si._debug:
#         print(gmsh_str)
#
#     return 0


# def generate_cylinder_geometry(si):
#     # Check if all parameters are set
#     abort_flag = False
#     for key, item in si._params_bempp["cylinder_params"].items():
#         if item is None:
#             print("Item 'aperture_params/{}' is not set in BEM++ parameters!".format(key))
#             abort_flag = True
#     if abort_flag:
#         return 1
#
#     h = si._params_bempp["h"]
#     radius = si._params_bempp["cylinder_params"]["radius"]
#     zmin = si._params_bempp["cylinder_params"]["zmin"]
#     zmax = si._params_bempp["cylinder_params"]["zmax"]
#     voltage = si._params_bempp["cylinder_params"]["voltage"]
#     electrode_type = "cylinder"
#
#     gmsh_str = """SetFactory("OpenCASCADE");
# Geometry.NumSubEdges = 100; // nicer display of curve
# h = {};
# """.format(h * 5)
#
#     gmsh_str += """
# Point(1) = {{ {}, {}, {}, h }};
# Point(2) = {{ {}, {}, {}, h }};
# Point(3) = {{ {}, {}, {}, h }};
# Point(4) = {{ {}, {}, {}, h }};
# Point(5) = {{ {}, {}, {}, h }};""".format(0.0, 0.0, zmin,
#                                           radius, 0.0, zmin,
#                                           0.0, radius, zmin,
#                                           -radius, 0.0, zmin,
#                                           0.0, -radius, zmin)
#
#     gmsh_str += """
# Point(6) = {{ {}, {}, {}, h }};
# Point(7) = {{ {}, {}, {}, h }};
# Point(8) = {{ {}, {}, {}, h }};
# Point(9) = {{ {}, {}, {}, h }};
# Point(10) = {{ {}, {}, {}, h }};""".format(0.0, 0.0, zmax,
#                                            radius, 0.0, zmax,
#                                            0.0, radius, zmax,
#                                            -radius, 0.0, zmax,
#                                            0.0, -radius, zmax)
#
#     gmsh_str += """
# Circle(1) = { 2, 1, 3 };
# Circle(2) = { 3, 1, 4 };
# Circle(3) = { 4, 1, 5 };
# Circle(4) = { 5, 1, 2 };
#
# Circle(5) = { 7, 6, 8 };
# Circle(6) = { 8, 6, 9 };
# Circle(7) = { 9, 6, 10 };
# Circle(8) = { 10, 6, 7 };
#
# Line(23) = { 2, 7 };
# Line(24) = { 3, 8 };
# Line(25) = { 4, 9 };
# Line(26) = { 5, 10 };
#
# Line Loop(1) = { 1, 2, 3, 4 };
# Line Loop(2) = { 5, 6, 7, 8 };
#
# Line Loop(5) = { 1, 24, -5, -23 };
# Line Loop(6) = { 2, 25, -6, -24 };
# Line Loop(7) = { 3, 26, -7, -25 };
# Line Loop(8) = { 4, 23, -8, -26 };
#
# Plane Surface(1) = { 1 };
# Plane Surface(2) = { 2 };
#
# Ruled Surface(3) = { 5 };
# Ruled Surface(4) = { 6 };
# Ruled Surface(5) = { 7 };
# Ruled Surface(6) = { 8 };
# """
#
#     si._variables_bempp["objects"][electrode_type] = {"gmsh_str": gmsh_str, "voltage": voltage}
#
#     if si._debug:
#         print(gmsh_str)
#
#     return 0


def generate_analytical_geometry(si):
    """
    The process of generating the geometry is as follows:
    Create the inner and outer surface of the spiral electrodes, shift the inside edges according to sigma,
    put everything together in one array.
    :return:
    """

    if si._variables_analytic["trj_design"] is None:
        print("No analytical design trajectory yet, generating...")
        si.generate_design_trajectory()

    print("Generating analytical geometry... ", end="")

    geos = []
    _ion = si._params_analytic["ion"]  # type: IonSpecies
    ns = si._params_analytic["ns"]  # type: int
    gap = si._params_analytic["gap"]  # type: float
    sigma = si._params_analytic["sigma"]  # type: float
    kp = si._variables_analytic["kp"]  # type: float
    b = si._variables_analytic["b"]  # type: np.ndarray
    cp = si._variables_analytic["c+"]  # type: float
    cm = si._variables_analytic["c-"]  # type: float
    trj_design = si._variables_analytic["trj_design"]  # type: np.ndarray

    for thickness in [0.0, si._params_analytic["dx"]]:

        end_distance = 2.0 * thickness + gap  # End to end distance of electrodes in (m)
        # TODO: Aspect ratio should be a parameter that can be set -PWCalculation of Spiral loflector Orbits
        aspect_ratio = 2.5 * (gap / end_distance)

        # Distance between electrodes at inflection angle theta
        d = end_distance * np.ones(ns) / (np.sqrt(1.0 + ((kp ** 2.0) * (np.sin(b)) ** 2.0)))

        # Save inner gap size vs deflection angle as class variable
        if thickness == 0.0:
            si._variables_analytic["d"] = d

        # x-component of the velocity vector
        vx = np.array(0.5 * _ion.v_m_per_s() * (np.sin(cp * b) + np.sin(cm * b)))

        # y-component of the velocity vector
        vy = np.array(-0.5 * _ion.v_m_per_s() * (np.cos(cp * b) - np.cos(cm * b)))

        # Rotation/flip
        if not ((si._variables_analytic["bf_design"] > 0.0) ^ (_ion.q() > 0.0)):
            if si.debug:
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
        si._variables_analytic["v_design"] = np.array([vx, vy, vz]).T

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

        # # Turn track of unit vectors
        # t1 = np.arange(0, 5, 0.01)
        # v_er = np.zeros((3, 3, np.size(t1)))

        # for i in range(3):
        #     for j in range(3):
        #         for k in range(np.size(t1)):
        #             v_er[i, j, k] = trj_design[ns - 1, j] + t1[k] * v_rh[i, ns - 1, j]

        # Construction of the electrodes
        edge_lines = np.zeros((5, ns, 3))

        xi = 0.5 * aspect_ratio * end_distance

        if si._params_analytic["rotation"] != 0.0:
            for i in range(si._params_analytic["ns"]):
                v_rh[0, i, :] = np.matmul(si._variables_analytic["rot"], v_rh[0, i, :])
                v_rh[1, i, :] = np.matmul(si._variables_analytic["rot"], v_rh[1, i, :])

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

    si._variables_analytic["geo"] = geo

    print("Done!")

    return si._variables_analytic["geo"]


def generate_numerical_geometry(si):
    # This is a slightly modified version of the normal analytical method
    if si._variables_analytic["trj_design"] is None:
        print("No numerical design trajectory yet, generating...")
        si.generate_design_trajectory()

    print("Generating numerical geometry... ", end="")

    geos = []
    _ion = si._params_analytic["ion"]  # type: IonSpecies
    ns = si._params_analytic["ns"]  # type: int
    gap = si._params_analytic["gap"]  # type: float
    sigma = si._params_analytic["sigma"]  # type: float
    kp = si._variables_analytic["kp"]  # type: float

    b = si._variables_analytic["b"]  # type: np.ndarray
    trj_design = si._variables_analytic["trj_design"]  # type: np.ndarray
    trj_vel = si._variables_analytic["trj_vel"]
    vx, vy, vz = trj_vel[:, 0], trj_vel[:, 1], trj_vel[:, 2]  # Unload the v components

    for thickness in [0.0, si._params_analytic["dx"]]:

        end_distance = 2.0 * thickness + gap  # End to end distance of electrodes in (m)
        aspect_ratio = 2.5 * (gap / end_distance)

        # Distance between electrodes at inflection angle theta
        d = end_distance * np.ones(ns) / (np.sqrt(1.0 + ((kp ** 2.0) * (np.sin(b)) ** 2.0)))

        # Save inner gap size vs deflection angle as class variable
        if thickness == 0.0:
            si._variables_analytic["d"] = d

        # Rotation/flip
        if not ((si._variables_analytic["bf_design"] > 0.0) ^ (_ion.q() > 0.0)):
            if si._debug:
                print("Flipping direction of cyclotron motion...", end="")
            vy = -vy

        v2 = np.sqrt((vx ** 2.0) + (vy ** 2.0))  # xy-magnitude of the velocity
        v3 = np.sqrt((vx ** 2.0) + (vy ** 2.0) + (vz ** 2.0))  # 3-d magnitude of the velocity vector (Should = v)

        # Save vx, vy, vz in as class variable in same format as trj_design
        si._variables_analytic["v_design"] = np.array([vx, vy, vz]).T

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

        # # Turn track of unit vectors
        # t1 = np.arange(0, 5, 0.01)
        # v_er = np.zeros((3, 3, np.size(t1)))

        # for i in range(3):
        #     for j in range(3):
        #         for k in range(np.size(t1)):
        #             v_er[i, j, k] = trj_design[ns - 1, j] + t1[k] * v_rh[i, ns - 1, j]

        # Construction of the electrodes
        edge_lines = np.zeros((5, ns, 3))

        xi = 0.5 * aspect_ratio * end_distance

        # TODO: This is a work in progress
        # if si._variables_optimization["x_rot"] is not None and si._params_exp["y_opt"]:
        #     xrot = np.array([[1.0, 0.0, 0.0],
        #                     [0.0, np.cos(si._variables_optimization["x_rot"]),
        #                      -np.sin(si._variables_optimization["x_rot"])],
        #                     [0.0, np.sin(si._variables_optimization["x_rot"]),
        #                      np.cos(si._variables_optimization["x_rot"])]])
        #     for i in range(si._params_analytic["ns"]):
        #         v_rh[0, i, :] = np.matmul(xrot, v_rh[0, i, :])
        #         v_rh[1, i, :] = np.matmul(xrot, v_rh[1, i, :])
        #     # print("Applied a {:.4f} rad x rotation.".format(si._variables_optimization["x_rot"]))

        if si._params_analytic["rotation"] != 0.0:
            for i in range(si._params_analytic["ns"]):
                v_rh[0, i, :] = np.matmul(si._variables_analytic["rot"], v_rh[0, i, :])
                v_rh[1, i, :] = np.matmul(si._variables_analytic["rot"], v_rh[1, i, :])

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

    si._variables_analytic["geo"] = geo

    print("Done!")

    return si._variables_analytic["geo"]


# def generate_spiral_electrode_geometry(si, electrode_type):
#     abort_flag = False
#     for key, item in si._params_bempp.items():
#         if item is None:
#             print("Item {} is not set in BEM++ parameters!".format(key))
#             abort_flag = True
#     if abort_flag:
#         return 1
#
#     if si._variables_analytic["geo"] is None:
#         print("No analytic geometry generated yet... starting now...")
#         si.generate_geometry()
#
#     si._variables_bempp["i"] = [1, 1, 1]
#
#     voltage = 0.0
#
#     if electrode_type is "anode":
#         k = 0
#         invert = True
#         voltage = si._params_analytic["volt"]
#     elif electrode_type is "cathode":
#         k = 5
#         invert = False
#         voltage = -si._params_analytic["volt"]
#
#     geo = si._variables_analytic["geo"]  # type: np.ndarray
#
#     h = si._params_bempp["h"]
#     gmsh_str = """
# Geometry.NumSubEdges = 100; // nicer display of curve
# """
#
#     ly = geo.shape[1]
#
#     for i in range(ly):
#
#         for j in [0, 2, 3, 1, 4]:
#
#             if i < 3 or i > ly - 4:
#
#                 gmsh_str += _gmsh_point(si, geo[j + k, i, :], h=h)
#
#             else:
#
#                 gmsh_str += _gmsh_point(si, geo[j + k, i, :])
#
#             si._variables_bempp["i"][0] += 1
#
#             if j != 0:
#                 gmsh_str += _gmsh_line(si, si._variables_bempp["i"][0] - 2,
#                                        si._variables_bempp["i"][0] - 1)
#                 si._variables_bempp["i"][1] += 1
#
#             if j == 4:
#
#                 gmsh_str += _gmsh_line(si, si._variables_bempp["i"][0] - 1,
#                                        si._variables_bempp["i"][0] - 5)
#                 si._variables_bempp["i"][1] += 1
#
#                 if i == 0:  # Electrode end #1
#
#                     gmsh_str += _gmsh_surface(si, [-(si._variables_bempp["i"][1] - 5),
#                                                    -(si._variables_bempp["i"][1] - 4),
#                                                    -(si._variables_bempp["i"][1] - 3),
#                                                    -(si._variables_bempp["i"][1] - 2),
#                                                    -(si._variables_bempp["i"][1] - 1)],
#                                               'Plane', invert)
#                 if i == ly - 1:  # Electrode end #2
#
#                     gmsh_str += _gmsh_surface(si, [(si._variables_bempp["i"][1] - 8),
#                                                    (si._variables_bempp["i"][1] - 6),
#                                                    (si._variables_bempp["i"][1] - 4),
#                                                    (si._variables_bempp["i"][1] - 2),
#                                                    (si._variables_bempp["i"][1] - 1)],
#                                               'Plane', invert)
#             if i > 0:
#
#                 gmsh_str += _gmsh_line(si, si._variables_bempp["i"][0] - 1,
#                                        si._variables_bempp["i"][0] - 6)
#                 si._variables_bempp["i"][1] += 1
#
#                 if i == 1:
#
#                     if j == 2:
#
#                         gmsh_str += _gmsh_surface(si, [si._variables_bempp["i"][1] - 8,
#                                                        -(si._variables_bempp["i"][1] - 1),
#                                                        -(si._variables_bempp["i"][1] - 2),
#                                                        si._variables_bempp["i"][1] - 3],
#                                                   'Ruled', invert)
#                     elif j == 3:
#
#                         gmsh_str += _gmsh_surface(si, [si._variables_bempp["i"][1] - 9,
#                                                        -(si._variables_bempp["i"][1] - 1),
#                                                        -(si._variables_bempp["i"][1] - 2),
#                                                        si._variables_bempp["i"][1] - 3],
#                                                   'Ruled', invert)
#                     elif j == 1:
#
#                         gmsh_str += _gmsh_surface(si, [si._variables_bempp["i"][1] - 10,
#                                                        -(si._variables_bempp["i"][1] - 1),
#                                                        -(si._variables_bempp["i"][1] - 2),
#                                                        si._variables_bempp["i"][1] - 3],
#                                                   'Ruled', invert)
#                     elif j == 4:
#
#                         gmsh_str += _gmsh_surface(si, [si._variables_bempp["i"][1] - 12,
#                                                        -(si._variables_bempp["i"][1] - 1),
#                                                        -(si._variables_bempp["i"][1] - 3),
#                                                        si._variables_bempp["i"][1] - 4],
#                                                   'Ruled', invert)
#                 elif i > 1:
#
#                     if j == 0:
#
#                         if i == 2:
#                             c = 0
#                         else:
#                             c = 1
#
#                         gmsh_str += _gmsh_surface(si, [(si._variables_bempp["i"][1] - 12 - c),
#                                                        (si._variables_bempp["i"][1] - 2),
#                                                        -(si._variables_bempp["i"][1] - 3),
#                                                        -(si._variables_bempp["i"][1] - 11)],
#                                                   'Ruled', invert)
#                     elif j != 0 and j != 4:
#
#                         gmsh_str += _gmsh_surface(si, [si._variables_bempp["i"][1] - 12,
#                                                        -(si._variables_bempp["i"][1] - 1),
#                                                        -(si._variables_bempp["i"][1] - 2),
#                                                        si._variables_bempp["i"][1] - 3],
#                                                   'Ruled', invert)
#                     elif j == 4:
#
#                         gmsh_str += _gmsh_surface(si, [si._variables_bempp["i"][1] - 13,
#                                                        -(si._variables_bempp["i"][1] - 1),
#                                                        -(si._variables_bempp["i"][1] - 3),
#                                                        si._variables_bempp["i"][1] - 4],
#                                                   'Ruled', invert)
#
#                         if i == ly - 1:
#                             gmsh_str += _gmsh_surface(si, [(si._variables_bempp["i"][1] - 12),
#                                                            (si._variables_bempp["i"][1] - 1),
#                                                            -(si._variables_bempp["i"][1] - 2),
#                                                            -(si._variables_bempp["i"][1] - 10)],
#                                                       'Ruled', invert)
#
#     si._variables_bempp["objects"][electrode_type] = {"gmsh_str": gmsh_str, "voltage": voltage}
#
#     if si._debug:
#         print(gmsh_str)
#
#     return 0

def get_norm_vec_and_angles_from_geo(geo):
    mid_vec_b = Vector(geo[8, -1, :] - geo[7, -1, :]).normalized()

    # tilt_angle is the angle of mid_vec_b with x/y plane
    tilt_angle = 0.5 * np.pi - mid_vec_b.angle_with(Vector(-Z_AXIS))

    # face angle is the angle of mid_vec_b projected into x/y plane with x/z plane
    temp_vec = Vector([mid_vec_b[0], mid_vec_b[1], 0.0])
    face_angle = 0.5 * np.pi - temp_vec.angle_with(Vector(Y_AXIS))

    return tilt_angle, face_angle


def generate_solid_assembly(si, apertures=None, cylinder=None):
    analytic_pars = si.analytic_parameters
    analytic_vars = si.analytic_variables
    bempp_pars = si.bempp_parameters
    bempp_vars = si.bempp_variables

    if apertures is not None:
        bempp_pars["make_aperture"] = apertures

    if cylinder is not None:
        bempp_pars["make_cylinder"] = cylinder

    if analytic_vars["geo"] is None:
        print("No geometry generated yet... starting now...")
        si.generate_geometry()

    #  --- Create Electrode objects
    abort_flag = False
    for key, item in bempp_pars.items():
        if item is None:
            print("Item {} is not set in BEM++ parameters!".format(key))
            abort_flag = True
    if abort_flag:
        return 1

    geo = analytic_vars["geo"]
    trj = analytic_vars["trj_design"]
    voltage = analytic_pars["volt"]
    h = bempp_pars["h"]

    anode = SIElectrode(name="SI Anode", voltage=voltage)
    anode.create_geo_str(raw_geo=geo, elec_type="anode", h=h, load=True)
    anode.color = "RED"

    cathode = SIElectrode(name="SI Cathode", voltage=-voltage)
    cathode.create_geo_str(raw_geo=geo, elec_type="cathode", h=h, load=True)
    anode.color = "BLUE"

    # Create an assembly holding all the electrodes
    assy = PyElectrodeAssembly("Spiral Inflector Assembly")
    assy.add_electrode(anode)
    assy.add_electrode(cathode)

    tilt_angle, face_angle = get_norm_vec_and_angles_from_geo(geo)

    if bempp_pars["make_housing"]:
        zmin = bempp_pars["housing_params"]["zmin"]
        zmax = bempp_pars["housing_params"]["zmax"]
        gap = bempp_pars["housing_params"]["gap"]
        thickness = bempp_pars["housing_params"]["thickness"]
        voltage = bempp_pars["housing_params"]["voltage"]

        housing = SIHousing(parent=si, name="Housing", voltage=voltage)

        angles = (tilt_angle, face_angle)
        housing.set_aperture_params(bempp_pars["aperture_params"])
        housing.set_aperture_rot_angles(angles)

        translate = np.array([0.0, 0.0, zmin])
        housing.set_translation(translate, absolute=True)
        s = housing.create_geo_str(geo=geo,
                                   trj=trj,
                                   zmin=zmin,
                                   zmax=zmax,
                                   gap=gap,
                                   thickness=thickness,
                                   h=h,
                                   load=True)

        housing.color = "GREEN"

        with open('testing.geo', 'w') as outfile:
            outfile.write(s)

        assy.add_electrode(housing)

    if si.bempp_parameters["make_aperture"]:
        # Base aperture parameters:
        voltage = bempp_pars["aperture_params"]["voltage"]
        dz = bempp_pars["aperture_params"]["thickness"]
        r = bempp_pars["aperture_params"]["radius"]
        a = bempp_pars["aperture_params"]["length"]
        b = bempp_pars["aperture_params"]["width"]
        t_gap = bempp_pars["aperture_params"]["top_distance"]
        b_gap = bempp_pars["aperture_params"]["bottom_distance"]
        hole_type = bempp_pars["aperture_params"]["hole_type"]

        # --- Entrance aperture --- #
        # TODO: May have to be rotated more with entrance of SI
        entrance_aperture = SIAperture(name="Entrance Aperture", voltage=voltage)

        # Calculate correct translation
        entrance_aperture.set_translation(np.array([0, 0, trj[0][2] - t_gap - 0.5 * dz]), absolute=True)
        entrance_aperture.set_rotation_angle_axis(angle=np.deg2rad(90.0), axis=Z_AXIS, absolute=True)

        # Create geo string and load
        entrance_aperture.create_geo_str(r=r, dz=dz, a=a, b=b, hole_type=hole_type, h=h, load=True)
        entrance_aperture.color = "GREEN"

        assy.add_electrode(entrance_aperture)

        # --- Exit aperture (rotated and shifted) --- #
        if not bempp_pars["make_housing"]:
            exit_aperture = SIAperture(name="Exit Aperture", voltage=0)

            # DEBUG: Display points used for angles
            # p1 = SIPointSphere(name="P1")
            # p1.create_geo_str(geo[4, -1, :])
            # p1.color = "GREEN"
            # p2 = SIPointSphere(name="P2")
            # p2.create_geo_str(geo[7, -1, :])
            # p2.color = "BLUE"
            # p3 = SIPointSphere(name="P3")
            # p3.create_geo_str(geo[8, -1, :])
            # p3.color = "RED"
            # p4 = SIPointSphere(name="P4")
            # p4.create_geo_str(geo[9, -1, :])
            # p4.color = "BLACK"
            #
            # assy.add_electrode(p1)
            # assy.add_electrode(p2)
            # assy.add_electrode(p3)
            # assy.add_electrode(p4)

            # Calculate correct rotation and translation
            tilt_angle, face_angle = get_norm_vec_and_angles_from_geo(geo)
            norm_vec = Vector(trj[-1] - trj[-2]).normalized()

            translate = np.array([trj[-1][0] + norm_vec[0] * b_gap,
                                  trj[-1][1] + norm_vec[1] * b_gap,
                                  0.0])
            exit_aperture.set_translation(translate, absolute=True)

            # Calculate correct rotation
            exit_aperture.set_rotation_angle_axis(angle=np.deg2rad(90.0), axis=X_AXIS, absolute=True)  # upright
            exit_aperture.set_rotation_angle_axis(angle=tilt_angle, axis=Y_AXIS, absolute=False)  # match tilt
            exit_aperture.set_rotation_angle_axis(angle=face_angle, axis=Z_AXIS, absolute=False)  # match exit

            # Create geo string and load
            exit_aperture.create_geo_str(r=r, dz=dz, a=a, b=b, hole_type=hole_type, h=h, load=True)
            exit_aperture.color = "GREEN"

            assy.add_electrode(exit_aperture)

    if bempp_pars["make_cylinder"]:
        # Base cylinder parameters:
        r = bempp_pars["cylinder_params"]["radius"]
        zmin = bempp_pars["cylinder_params"]["zmin"]
        zmax = bempp_pars["cylinder_params"]["zmax"]
        voltage = bempp_pars["cylinder_params"]["voltage"]

        outer_cylinder = SICylinder(name="Outer Cylinder", voltage=voltage)
        translate = np.array([0.0, 0.0, zmin])
        outer_cylinder.set_translation(translate, absolute=True)
        outer_cylinder.create_geo_str(r=r, dz=zmax - zmin, h=h, load=True)

        assy.add_electrode(outer_cylinder)

    if si.debug:
        assy.show(show_screen=True)

    bempp_vars["objects"] = assy

    si.analytic_parameters = analytic_pars
    si.analytic_variables = analytic_vars
    si.bempp_parameters = bempp_pars
    si.bempp_variables = bempp_vars

    return assy


def generate_meshed_model(si, apertures=None, cylinder=None):
    # TODO: Think about surface normals for electrodes and outer cylinder!

    generate_solid_assembly(si, apertures, cylinder)

    bempp_vars = si.bempp_variables

    assy = bempp_vars["objects"]

    leaf_view = assy.get_bempp_mesh()

    bempp_vars["full mesh"] = {"verts": leaf_view["verts"],
                               "elems": leaf_view["elems"],
                               "domns": leaf_view["domns"]}

    if si.debug:
        _full_mesh = bempp.api.grid_from_element_data(leaf_view["verts"],
                                                      leaf_view["elems"],
                                                      leaf_view["domns"])
        _full_mesh.plot()

    si.bempp_variables = bempp_vars

    return bempp_vars["full mesh"]


def export_aperture_geometry(si, fname="aperture_macro.ivb"):
    geo = si._variables_analytic["geo"] * 100.0  # Scaling for inventor
    trj = si._variables_analytic["trj_design"] * 100.0  # type: np.ndarray
    thickness = si._params_bempp["aperture_params"]["thickness"] * 100.0
    radius = si._params_bempp["aperture_params"]["radius"] * 100.0
    length = si._params_bempp["aperture_params"]["length"] * 100.0
    width = si._params_bempp["aperture_params"]["width"] * 100.0
    aperture_distance_top = si._params_bempp["aperture_params"]["top_distance"] * 100.0
    aperture_distance_bottom = si._params_bempp["aperture_params"]["bottom_distance"] * 100.0
    # voltage = si._params_bempp["aperture_params"]["voltage"]

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

        offset = norm_vec * aperture_distance_top

        if i == -1:
            offset = -norm_vec * aperture_distance_bottom

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

        extrude_dir = "kNegativeExtentDirection"
        if i == -1:
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
    import os
    with open(os.path.join(si._outp_folder, fname), "w") as outfile:
        outfile.writelines(aperture_string)


def export_electrode_geometry(si, fname="electrode_macro.ivb"):
    geo = si._variables_analytic["geo"] * 100.0  # Fix scaling in inventor

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

            for i in range(si._params_analytic["ns"]):
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
    import os
    with open(os.path.join(si._outp_folder, fname), "w") as of:

        of.write(header_text + electrode_texts[0] + electrode_texts[1] + footer_text)

    print("Done!")


def save_geo_files(si, filename=None):
    import os
    from dans_pymodules import FileDialog
    if filename is None:
        fd = FileDialog()
        folder, _ = fd.get_filename("folder")
        if folder is None:
            return 0
    else:
        folder = os.path.split(filename)[0]

    for name, electrode in si._variables_bempp["objects"].items():
        if si._debug:
            print("Working on {}".format(name))

        filename = os.path.join(folder, "{}.geo".format(name))
        with open(filename, 'w') as of:
            of.write(electrode["gmsh_str"])

    return 0
