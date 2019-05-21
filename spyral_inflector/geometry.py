from dans_pymodules import *
from py_electrodes.py_electrodes import *

X_AXIS = np.array([1, 0, 0], float)
Y_AXIS = np.array([0, 1, 0], float)
Z_AXIS = np.array([0, 0, 1], float)

HAVE_BEMPP = False
try:
    import bempp.api
    from bempp.api.shapes.shapes import __generate_grid_from_geo_string as generate_from_string

    HAVE_BEMPP = True
except ImportError:
    bempp = None

HAVE_FENICS = False
try:
    import fenics as fn

    HAVE_FENICS = True
except:
    fn = None

HAVE_MESHIO = False
try:
    import meshio

    HAVE_MESHIO = True
except ImportError:
    meshio = None


class SIAperture(PyElectrode):
    def __init__(self, parent=None, name="New Aperture", voltage=0, offset=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture
        self._offset = offset

    def create_geo_str(self, r, dz, a, b, translation=None, rotation=None, hole_type="ellipse", h=0.005, load=True,
                       header=True):
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
        :param header: Flag whether to include the header for the geo string.
        :return gmsh_str: the string object for gmsh
        """

        offset = self._offset

        if translation is None:
            translation = np.array([0.0, 0.0, 0.0])

        if rotation is None:
            rotation = np.array([0.0, 0.0, 0.0])

        if header:
            geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
// Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)
        else:
            geo_str = ""

        geo_str += "// Base plate\n"
        geo_str += "Cylinder({}) = {{ 0, 0, {}, 0, 0, {}, {}, 2 * Pi }};\n\n".format(1 + offset, -0.5 * dz, dz, r)

        geo_str += "// Tool to subtract\n"
        if hole_type == "rectangle":
            geo_str += "Box({}) = {{ {}, {}, {}, {}, {}, {} }};\n\n".format(2 + offset, -0.5 * a, -0.5 * b, -dz, a, b,
                                                                            2 * dz)
        elif hole_type == "ellipse":
            geo_str += "Disk({}) = {{ 0, 0, {}, {}, {} }};\n".format(100 + offset, -dz, 0.5 * a, 0.5 * b)
            geo_str += "Extrude {{ 0, 0, {} }} {{ Surface{{ {} }}; }}\n".format(100 + offset, 2 * dz)
        else:
            print("Don't understand hole type {}!".format(hole_type))
            return 1

        geo_str += "\nBooleanDifference({}) = {{ Volume{{ {} }}; Delete; }}{{ Volume{{ {} }}; Delete; }};\n".format(
            3 + offset, 1 + offset, 2 + offset)

        geo_str += "Rotate {{ {{ 1.0, 0.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            rotation[0], 3 + offset)
        geo_str += "Rotate {{ {{ 0.0, 1.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            rotation[1], 3 + offset)
        geo_str += "Rotate {{ {{ 0.0, 0.0, 1.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            rotation[2], 3 + offset)

        geo_str += "Translate {{ {}, {}, {} }} {{ Volume{{ {} }}; }}\n".format(translation[0],
                                                                               translation[1],
                                                                               translation[2],
                                                                               3 + offset)

        # Call function in PyElectrode module we inherit from if 'load' is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


class SIPointSphere(PyElectrode):
    def __init__(self, parent=None, name="New Aperture", voltage=0, offset=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture
        self._offset = offset

    def create_geo_str(self, center, r=0.001, h=0.005, load=True, header=True):
        """

        Creates the geo string for a small sphere to show important vertices

        :param center: numpy array, list or tuple containing 3 floats forthe center of the sphere
        :param r: sphere radius - default is 1 mm
        :param h: desired mesh resolution
        :param load: Flag whether to also load from geo string directly.
                     Cave: If False, geo str will not be saved internally!
        :param header: Flag whether to include the header for the geo string.
        :return gmsh_str: the string object for gmsh
        """

        offset = self._offset

        center = np.asarray(center)
        assert center.shape == (3,), "Got wrong dimension of {} for center, should be (3, )".format(center.shape)

        if header:
            geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
// Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)
        else:
            geo_str = ""

        geo_str += "// Base plate\n"
        geo_str += "Sphere({}) = {{ {}, {}, {}, {} }};\n".format(1 + offset, center[0], center[1], center[2], r)

        # Call function in PyElectrode module we inherit from if 'load' is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


class SICylinder(PyElectrode):
    def __init__(self, parent=None, name="New Cylinder", voltage=0, offset=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture
        self._offset = offset

    def create_geo_str(self, r, zmin, zmax, h=0.0075, load=True, header=True):
        # TODO: This docstring is incorrect -PW
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
        :param header: Flag whether to include the header for the geo string.
        :return gmsh_str: the string object for gmsh
        """

        offset = self._offset

        if header:
            geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
// Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)
        else:
            geo_str = ""

        geo_str += "// Cylinder\n"
        geo_str += "Cylinder({}) = {{ 0, 0, {}, 0, 0, {}, {}, 2 * Pi }};\n\n".format(1 + offset, zmin, zmax - zmin, r)

        # Call function in PyElectrode module we inherit from if load is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


class SIElectrode(PyElectrode):
    def __init__(self, parent=None, name="New Spiral Electrode", voltage=10000, offset=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture
        self._offset = offset

    def create_geo_str(self, raw_geo, elec_type, h=0.005, load=True, header=True):
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
        :param header: Flag whether to include the header for the geo string.
        :return gmsh_str: the string object for gmsh
        """

        f = self._offset  # Normally I use offset = self._offset, but there are a lot of uses of it so to keep
        # it short, I will use f here. -PW

        if elec_type not in ["anode", "cathode"]:
            print("SIElectrode could not understand electrode type {}. Must be 'Anode' or 'Cathode'".format(type))
            return 1

        if header:
            geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
// Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)
        else:
            geo_str = ""

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
                geo_str += "Point({}) = {{ {}, {}, {}, {} }};\n".format(new_pt + f,
                                                                        raw_geo[i + k, j, 0],
                                                                        raw_geo[i + k, j, 1],
                                                                        raw_geo[i + k, j, 2],
                                                                        h)
                new_pt += 1

            # For each section, add the lines
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 0 + f, (j * 5) + 4 + f, (j * 5) + 2 + f)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 1 + f, (j * 5) + 3 + f, (j * 5) + 1 + f)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 2 + f, (j * 5) + 3 + f, (j * 5) + 4 + f)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 3 + f, (j * 5) + 5 + f, (j * 5) + 2 + f)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 4 + f, (j * 5) + 1 + f, (j * 5) + 5 + f)

            new_ln += 5

            geo_str += "Wire({}) = {{ {}, {}, {}, {}, {} }};\n\n".format(new_loop + f,
                                                                         (j * 5) + 3 + f,
                                                                         (j * 5) + 2 + f,
                                                                         (j * 5) + 5 + f,
                                                                         (j * 5) + 4 + f,
                                                                         (j * 5) + 1 + f)

            new_loop += 1

        geo_str += "Ruled ThruSections({}) = {{ {}:{} }};".format(new_vol + f, 1 + f, new_loop - 1 + f)

        new_vol += 1

        # Call function in PyElectrode module we inherit from if load is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


# Geometrically, trajectories have much in common with electrodes...
class SITrajectory(PyElectrode):

    def __init__(self, parent=None, name="New Spiral Electrode", voltage=0, offset=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture
        self._offset = offset

    # --- Override some of the PyElectrode functions that don't make sense for a wire --- #
    @staticmethod
    def points_inside(_points):
        # By definition no points are 'inside' a PyWire
        return np.zeros(_points.shape, bool)

    @staticmethod
    def generate_mesh(brep_h=0.0):
        print("Can't generate a mesh from a PyWire")
        return 1

    def create_geo_str(self, points, max_points, load=True, header=True):
        """
        Create a geo string for gmsh
        :param points: np array of points along the trajectory
        :param max_points: maximum number of points to use if max is larger than num points all are used
        :param load: immediately load the geo str as an occ object
        :param header: Flag whether to include the header for the geo string.
        :return:
        """

        offset = self._offset

        points = np.asarray(points)

        assert points.ndim == 2 and points[0, :].shape == (3,), "points have wrong shape = {}".format(points.shape)

        # Reduce number of points to use in spline to max_points
        if max_points is not None:
            points = points[::int(np.ceil(len(points)) / max_points), :]

        if header:
            geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
"""
        else:
            geo_str = ""

        new_pt = 1
        new_ln = 1

        geo_str += "// Center Spline:\n"
        for _x, _y, _z in points:
            geo_str += "Point({}) = {{ {}, {}, {} }};\n".format(new_pt + offset, _x, _y, _z)
            new_pt += 1

        geo_str += """
Spline({}) = {{ {}:{} }}; 
""".format(new_ln + offset, 1 + offset, new_pt - 1 + offset)

        # Immediately delete the points used up in the spline
        geo_str += "Recursive Delete {{ Point{{ {}:{} }}; }}\n".format(1 + offset, new_pt - 1 + offset)

        # Call function in PyElectrode module we inherit from if load is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


class SIHousing(PyElectrode):

    def __init__(self, parent=None, name="Spiral Inflector Housing", voltage=0, offset=0, experimental=False):
        super().__init__(name=name, voltage=voltage)

        assert parent is not None, "This class requires a parent."

        self._parent = parent
        self._offset = offset

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

        total_pts = []
        circle_pts = []
        circle_res = 8
        for i in range(circle_res):
            circle_pts.append(gap * np.array([np.cos(i * 2 * np.pi / circle_res),
                                              np.sin(i * 2 * np.pi / circle_res)]))

        for pt in hull_pts:
            for cpt in circle_pts:
                total_pts.append(pt + cpt)

        total_pts = np.array(total_pts)

        new_hull_inner = ConvexHull(total_pts)
        new_hull_pts_inner = total_pts[new_hull_inner.vertices, :]

        total_pts = []
        circle_pts = []
        circle_res = 8
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
            # TODO: This may only work for a positive tilt angle, if it's negative then
            # TODO: the 2, 8 indices for geo will be different (3 and 7?) -PW
            tilt_angle, face_angle = get_norm_vec_and_angles_from_geo(geo)
            face_vector = np.array([np.cos(face_angle), np.sin(face_angle), 0.0])
            norm_vector = np.cross(face_vector, np.array([0.0, 0.0, 1.0]))

            # new_point_a_in = geo[2, -1, :2] + norm_vector[:2] * gap
            new_point_b_in = geo[8, -1, :2] + norm_vector[:2] * gap

            # new_point_a_out = geo[2, -1, :2] + norm_vector[:2] * (gap + thickness)
            new_point_b_out = geo[8, -1, :2] + norm_vector[:2] * (gap + thickness)

            # new_hull_pts_inner = np.vstack([new_hull_pts_inner, new_point_a_in])
            new_hull_pts_inner = np.vstack([new_hull_pts_inner, new_point_b_in])

            # new_hull_pts_outer = np.vstack([new_hull_pts_outer, new_point_a_out])
            new_hull_pts_outer = np.vstack([new_hull_pts_outer, new_point_b_out])

        pts_in = self.sort_points_by_angle(new_hull_pts_inner)
        pts_out = self.sort_points_by_angle(new_hull_pts_outer)

        return pts_in, pts_out

    def sort_points_by_angle(self, points):
        angles = []
        for i, point in enumerate(points):
            theta = np.arctan2(point[1], point[0])
            angles.append((i, theta))

        angles.sort(key=lambda tup: tup[1])

        new_points = []
        for tup in angles:
            pt = points[tup[0]]
            new_points.append(pt)

        return np.array(new_points)

    def create_geo_str(self, geo, trj, zmin, zmax, span, gap, thickness, h=0.005, load=True, header=True):
        # TODO: Doc string -PW

        offset = self._offset

        pts_in, pts_out = self.gen_convex_hull(geo, gap, thickness)

        dz = self._aperture_params["thickness"]
        r = self._aperture_params["radius"]
        a = self._aperture_params["length"]
        b = self._aperture_params["width"]
        t_gap = self._aperture_params["top_distance"]
        b_gap = self._aperture_params["bottom_distance"]

        if span:
            zmin = np.min(geo[:, :, 2])
            zmax = np.max(geo[:, :, 2])

        norm_vec = Vector(trj[-1] - trj[-2]).normalized()

        translate = np.array([trj[-1][0] + norm_vec[0] * b_gap,
                              trj[-1][1] + norm_vec[1] * b_gap,
                              0.0])

        hole_type = self._aperture_params["hole_type"]

        # TODO: ToleranceBoolean is important for the aperture hole. This was tested from 1T to about 2.4T and
        # TODO: seems to work... 3T did not. -PW

        if header:
            geo_str = """SetFactory("OpenCASCADE");
// Geometry.NumSubEdges = 100; // nicer display of curve
Geometry.ToleranceBoolean = 1E-5;
Geometry.Tolerance = 1E-10;
// Mesh.CharacteristicLengthMax = {};  // maximum mesh size""".format(h)
        else:
            geo_str = "Geometry.ToleranceBoolean = 1E-5;\n"

        geo_str += "// Outside points\n"
        n_pts_out = np.shape(pts_out)[0]

        for i, pt in enumerate(pts_out):
            geo_str += "Point({}) = {{ {}, {}, 0, {}}};\n".format(i + offset, pt[0], pt[1], h)

        geo_str += "// Outside lines\n"
        n_lines_out = n_pts_out  # Should be the same number

        for i in range(n_lines_out - 1):
            geo_str += "Line({}) = {{ {}, {} }};\n".format(i + offset, i + offset, i + 1 + offset)
        geo_str += "Line({}) = {{ {}, {} }};\n".format(i + 1 + offset, n_pts_out - 1 + offset,
                                                       0 + offset)  # Connect last point to first point

        geo_str += "// Inside points\n"
        n_pts_in = np.shape(pts_in)[0]

        for i, pt in enumerate(pts_in):
            geo_str += "Point({}) = {{ {}, {}, 0}};\n".format(i + n_pts_out + offset, pt[0], pt[1])

        geo_str += "// Inside lines\n"
        n_lines_in = n_pts_in  # Should be the same number

        for i in range(n_lines_in - 1):
            geo_str += "Line({}) = {{ {}, {} }};\n".format(i + n_lines_out + offset, i + n_pts_out + offset,
                                                           i + 1 + n_pts_out + offset)
        # Connect last point to first point
        geo_str += "Line({}) = {{ {}, {} }};\n".format(n_lines_in + n_lines_out + offset - 1,
                                                       i + 1 + n_pts_out + offset,
                                                       n_pts_out + offset)

        geo_str += "Line Loop({}) = {{ {}:{} }};\n".format(1 + offset, 0 + offset, n_lines_out - 1 + offset)
        geo_str += "Line Loop({}) = {{ {}:{} }};\n".format(2 + offset, n_lines_out + offset,
                                                           n_lines_out + n_lines_in - 1 + offset)

        geo_str += "Plane Surface({}) = {{ {}, {} }};\n".format(1 + offset, 1 + offset, 2 + offset)

        num_wire_points = 12
        # This will be extruded from 0 to zmin, then translated.
        z_values = np.linspace(0.0, zmax - zmin, num_wire_points)
        point_index = n_pts_in + n_pts_out + offset
        for k in range(num_wire_points):
            geo_str += "Point({}) = {{ 0.0, 0.0, {} }};\n".format(point_index + k, z_values[k])

        # Create the lines that connect these new points
        line_index = n_lines_in + n_lines_out + offset + 1
        for m in range(num_wire_points - 1):
            geo_str += "Line({}) = {{ {}, {} }};\n".format(line_index + m,
                                                           point_index + m, point_index + m + 1)

        geo_str += "Wire({}) = {{ {}:{} }};\n".format(3 + offset,
                                                      n_lines_in + n_lines_out + offset + 1,
                                                      n_lines_in + n_lines_out + offset + num_wire_points - 1)

        geo_str += "housing_out[] = Extrude {{ Surface{{ {} }}; }} Using Wire {{ {} }};\n".format(1 + offset,
                                                                                                  3 + offset)
        geo_str += "Translate {{ 0, 0, {} }} {{ Volume{{ housing_out[] }}; }}\n".format(zmin)

        geo_str += "Recursive Delete {{ Point{{ {}:{} }}; }}\n".format(point_index, point_index + num_wire_points - 1)
        geo_str += "Recursive Delete {{ Line{{ {}:{} }}; }}\n".format(line_index, line_index + num_wire_points - 2)

        geo_str += "// Tool to subtract\n"
        if hole_type == "rectangle":
            if self._experimental:
                geo_str += "Box({}) = {{ {}, {}, {}, {}, {}, {} }};\n\n".format(2 + offset,
                                                                                -0.5 * a, -0.5 * b,
                                                                                -0.025, a,
                                                                                b, 0.05)
            else:
                geo_str += "Box({}) = {{ {}, {}, {}, {}, {}, {} }};\n\n".format(2 + offset,
                                                                                -0.5 * a, -0.5 * b,
                                                                                0.0, a,
                                                                                b, 0.1)
        elif hole_type == "ellipse":  # TODO: The ellipse doesn't work -PW
            geo_str += "Disk ({}) = {{ 0, 0, 0, {}, {} }};\n".format(500 + offset, 0.5 * a + 5E-9, 0.5 * b + 5E-9)
            geo_str += "disk_out[] = Extrude {{ 0, 0, {} }} {{ Surface{{ {} }}; }};\n".format(2 + offset, 0.1,
                                                                                              500 + offset)

        geo_str += "Rotate {{ {{ 1.0, 0.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            np.pi / 2.0, 2 + offset)
        geo_str += "Rotate {{ {{ 0.0, 1.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            self._tilt_angle, 2 + offset)
        geo_str += "Rotate {{ {{ 0.0, 0.0, 1.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            self._face_angle, 2 + offset)

        geo_str += "Translate {{ {}, {}, {} }} {{ Volume{{ {} }}; }}\n".format(translate[0],
                                                                               translate[1],
                                                                               0.0,
                                                                               2 + offset)

        geo_str += "BooleanDifference({}) = {{ Volume {{ housing_out[] }}; Delete; }}{{ Volume {{ {} }}; Delete; }};\n".format(
            50 + offset, 2 + offset)

        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str


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
    aspect_ratio = si._params_analytic["aspect_ratio"]  # type: float
    kp = si._variables_analytic["kp"]  # type: float
    b = si._variables_analytic["b"]  # type: np.ndarray
    cp = si._variables_analytic["c+"]  # type: float
    cm = si._variables_analytic["c-"]  # type: float
    trj_design = si._variables_analytic["trj_design"]  # type: np.ndarray

    for thickness in [0.0, si._params_analytic["dx"]]:

        end_distance = 2.0 * thickness + gap  # End to end distance of electrodes in (m)

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

        xi = 0.5 * aspect_ratio * gap

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
    aspect_ratio = si._params_analytic["aspect_ratio"]  # type: float
    kp = si._variables_analytic["kp"]  # type: float

    b = si._variables_analytic["b"]  # type: np.ndarray
    trj_design = si._variables_analytic["trj_design"]  # type: np.ndarray
    trj_vel = si._variables_analytic["trj_vel"]
    vx, vy, vz = trj_vel[:, 0], trj_vel[:, 1], trj_vel[:, 2]  # Unload the v components

    for thickness in [0.0, si._params_analytic["dx"]]:

        end_distance = 2.0 * thickness + gap  # End to end distance of electrodes in (m)

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

        xi = 0.5 * aspect_ratio * gap

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


def get_norm_vec_and_angles_from_geo(geo):
    mid_vec_b = Vector(geo[8, -1, :] - geo[7, -1, :]).normalized()
    # tilt_angle is the angle of mid_vec_b with x/y plane

    tilt_angle = 0.5 * np.pi - mid_vec_b.angle_with(Vector(-Z_AXIS))

    # face angle is the angle of mid_vec_b projected into x/y plane with x/z plane
    temp_vec = Vector([mid_vec_b[0], mid_vec_b[1], 0.0])
    face_angle = 0.5 * np.pi - temp_vec.angle_with(Vector(Y_AXIS))

    return tilt_angle, face_angle


def generate_vacuum_space(si):
    # TODO: Clean up
    assert si.numerical_parameters["make_cylinder"], "You need a cylinder/boundary to create the vacuum space!"

    numerical_vars = si.numerical_variables
    numerical_pars = si.numerical_parameters

    assy = numerical_vars["objects"]

    master_geo_str = "// Full .geo file for fenics mesh generation\n"

    if numerical_pars["make_aperture"]:
        make_top_aperture = True
        if numerical_pars["make_housing"]:
            make_housing = True
            make_bottom_aperture = False
        else:
            make_housing = False
            make_bottom_aperture = True
    else:
        make_top_aperture = False
        make_housing = False
        make_bottom_aperture = False

    for _, electrode in assy.electrodes.items():
        master_geo_str += electrode._geo_str

    master_geo_str += """
//    anode_offset = 0
//    cathode_offset = 1000
//    housing_offset = 5000
//    exit_offset = 5000
//    entrance_offset = 3000

// Anode
anode_boundary[] = Boundary { Volume{ 1 }; };
N_anode = #anode_boundary[];
For i In {0:N_anode-1}
    Physical Surface (i + 1) = { anode_boundary[i] };
EndFor

// Cathode
cathode_boundary[] = Boundary { Volume{ 1001 }; };
N_cathode = #cathode_boundary[];
For k In {0:N_cathode-1}
    Physical Surface (1000 + k) = { cathode_boundary[k] };
EndFor
"""

    master_geo_str += """
// Vacuum Cylinder
vacuum_boundary[] = Boundary { Volume{ 4001 }; };
N_vacuum = #vacuum_boundary[];
For j In {0:N_vacuum-1}
    Physical Surface (4000 + j) = { vacuum_boundary[j] };
EndFor

// Surface Loop (2) = { vacuum_boundary };
// Surface Loop (1) = { anode_boundary };
// Surface Loop (3) = { cathode_boundary };
"""

    if make_bottom_aperture or make_housing:
        if make_housing:
            ap_id = 5050
        else:
            ap_id = 5003

        master_geo_str += """
// Bottom/Exit Aperture or Housing
exit_boundary[] = Boundary {{ Volume{{ {} }}; }};""".format(ap_id)

        master_geo_str += """
N_exit = #exit_boundary[];
For k In {0:N_exit-1}
    Physical Surface (5000 + k) = { exit_boundary[k] };
EndFor
"""
    if make_top_aperture:
        master_geo_str += """
// Top/Entrance Aperture
entrance_boundary[] = Boundary { Volume{ 3003 }; };
N_entrance = #entrance_boundary[];
For k In {0:N_entrance-1}
    Physical Surface (3000 + k) = { entrance_boundary[k] };
EndFor
"""

    master_geo_str += "Delete{ Volume{1, 1001, 4001"

    if make_bottom_aperture or make_housing:
        master_geo_str += ", {}".format(ap_id)
    if make_top_aperture:
        master_geo_str += ", 3003"

    master_geo_str += "}; }\n"

    master_geo_str += "Volume (1) = {3, 1, 2"

    if make_bottom_aperture or make_housing:
        master_geo_str += ", 4"
    if make_top_aperture:
        if make_bottom_aperture or make_housing:
            master_geo_str += ", 5"
        else:
            master_geo_str += ", 4"

    master_geo_str += "};\n"

    master_geo_str += """

Physical Volume(1) = { 1 };

edge_list[] = Boundary { Volume{ 1 }; };
N_edges = #edge_list[];
Field[1] = Distance;
Field[1].NNodesByEdge = 100;
Field[1].EdgesList = {1:N_edges};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.0015;
Field[2].LcMax = 0.01;
Field[2].DistMin = 0.003;
Field[2].DistMax = 0.01;
Background Field = 2;

// Mesh.OptimizeNetgen = 1;

"""

    with open('master_geometry.geo', 'w') as outfile:
        outfile.write(master_geo_str)


def generate_solid_assembly(si, apertures=None, cylinder=None):
    analytic_pars = si.analytic_parameters
    analytic_vars = si.analytic_variables
    numerical_pars = si.numerical_parameters
    numerical_vars = si.numerical_variables
    solver = si.solver

    if apertures is not None:
        numerical_pars["make_aperture"] = apertures

    if cylinder is not None:
        numerical_pars["make_cylinder"] = cylinder

    if analytic_vars["geo"] is None:
        print("No geometry generated yet... starting now...")
        si.generate_geometry()

    #  --- Create Electrode objects
    abort_flag = False
    for key, item in numerical_pars.items():
        if item is None:
            print("Item {} is not set in BEM++ parameters!".format(key))
            abort_flag = True
    if abort_flag:
        return 1

    geo = analytic_vars["geo"]
    trj = analytic_vars["trj_design"]
    voltage = analytic_pars["volt"]
    h = numerical_pars["h"]

    # Variables for fenics solving, won't affect anything BEMPP related (ideally) -PW
    anode_offset = 0
    cathode_offset = 1000
    housing_offset = 5000
    exit_offset = 5000
    entrance_offset = 3000
    cylinder_offset = 4000

    anode = SIElectrode(name="SI Anode", voltage=voltage, offset=anode_offset)
    anode.create_geo_str(raw_geo=geo, elec_type="anode", h=h, load=True, header=True)
    anode.color = "RED"

    cathode = SIElectrode(name="SI Cathode", voltage=-voltage, offset=cathode_offset)
    cathode.create_geo_str(raw_geo=geo, elec_type="cathode", h=h, load=True, header=True)
    anode.color = "BLUE"

    # Create an assembly holding all the electrodes
    assy = PyElectrodeAssembly("Spiral Inflector Assembly")
    assy.add_electrode(anode)
    assy.add_electrode(cathode)


    if numerical_pars["make_housing"]:
        zmin = numerical_pars["housing_params"]["zmin"]
        zmax = numerical_pars["housing_params"]["zmax"]
        gap = numerical_pars["housing_params"]["gap"]
        thickness = numerical_pars["housing_params"]["thickness"]
        voltage = numerical_pars["housing_params"]["voltage"]
        experimental = numerical_pars["housing_params"]["experimental"]
        span = numerical_pars["housing_params"]["span"]

        housing = SIHousing(parent=si, name="Housing", voltage=voltage,
                            offset=housing_offset, experimental=experimental)

        tilt_angle, face_angle = get_norm_vec_and_angles_from_geo(geo)
        angles = (tilt_angle, face_angle)
        housing.set_aperture_params(numerical_pars["aperture_params"])
        housing.set_aperture_rot_angles(angles)

        # translate = np.array([0.0, 0.0, zmin])
        # housing.set_translation(translate, absolute=True)
        s = housing.create_geo_str(geo=geo,
                                   trj=trj,
                                   zmin=zmin,
                                   zmax=zmax,
                                   span=span,
                                   gap=gap,
                                   thickness=thickness,
                                   h=h,
                                   load=True,
                                   header=True)

        # with open('housing_geo_str.geo', 'w') as f:
        #     f.write(s)

        housing.color = "GREEN"

        assy.add_electrode(housing)

    if si.numerical_parameters["make_aperture"]:
        # Base aperture parameters:
        voltage = numerical_pars["aperture_params"]["voltage"]
        dz = numerical_pars["aperture_params"]["thickness"]
        r = numerical_pars["aperture_params"]["radius"]
        a = numerical_pars["aperture_params"]["length"]
        b = numerical_pars["aperture_params"]["width"]
        t_gap = numerical_pars["aperture_params"]["top_distance"]
        b_gap = numerical_pars["aperture_params"]["bottom_distance"]
        hole_type = numerical_pars["aperture_params"]["hole_type"]

        # --- Entrance aperture --- #
        # TODO: May have to be rotated more with entrance of SI
        entrance_aperture = SIAperture(name="Entrance Aperture", voltage=voltage, offset=entrance_offset)

        # Calculate correct translation and rotation
        translation = np.array([0, 0, trj[0][2] - t_gap - 0.5 * dz])
        rotation = np.array([0.0, 0.0, np.deg2rad(90.0)])

        electrode_angle_vec = geo[7, 0, :2] - geo[8, 0, :2]
        electrode_angle_vec /= np.linalg.norm(electrode_angle_vec)
        electrode_angle = np.arctan(electrode_angle_vec[1] / electrode_angle_vec[0]) - np.pi / 2.0
        rotation += np.array([0.0, 0.0, electrode_angle])

        if solver == "bempp":
            entrance_aperture.set_translation(translation, absolute=True)
            entrance_aperture.set_rotation_angle_axis(angle=rotation[2], axis=Z_AXIS, absolute=True)

        # Create geo string and load
        if solver == "bempp":
            translation = np.array([0.0, 0.0, 0.0])
            rotation = np.array([0.0, 0.0, 0.0])
        entrance_aperture.create_geo_str(r=r, dz=dz, a=a, b=b, translation=translation, rotation=rotation,
                                         hole_type=hole_type,
                                         h=h, load=True, header=True)
        entrance_aperture.color = "GREEN"

        assy.add_electrode(entrance_aperture)

        # --- Exit aperture (rotated and shifted) --- #
        if not numerical_pars["make_housing"]:
            exit_aperture = SIAperture(name="Exit Aperture", voltage=0, offset=exit_offset)

            # DEBUG: Display points used for angles
            # p1 = SIPointSphere(name="P1")
            # p1.create_geo_str(geo[4, 0, :])
            # p1.color = "GREEN"
            # p2 = SIPointSphere(name="P2")
            # p2.create_geo_str(geo[7, 0, :])
            # p2.color = "BLUE"
            # p3 = SIPointSphere(name="P3")
            # p3.create_geo_str(geo[8, 0, :])
            # p3.color = "RED"
            # p4 = SIPointSphere(name="P4")
            # p4.create_geo_str(geo[9, 0, :])
            # p4.color = "BLACK"
            #
            # assy.add_electrode(p1)
            # assy.add_electrode(p2)
            # assy.add_electrode(p3)
            # assy.add_electrode(p4)

            # Calculate correct rotation and translation
            tilt_angle, face_angle = get_norm_vec_and_angles_from_geo(geo)

            tvec = geo[7, -1, :] - geo[8, -1, :]
            tvec /= np.linalg.norm(tvec)
            norm_vec = np.cross(tvec, np.array([0.0, 0.0, -1.0]))

            # norm_vec = Vector(trj[-1] - trj[-2]).normalized()

            translation = np.array([trj[-1][0] + norm_vec[0] * b_gap,
                                    trj[-1][1] + norm_vec[1] * b_gap,
                                    0.0])

            rotation = np.array([np.deg2rad(90.0),
                                 tilt_angle,
                                 face_angle])

            # Calculate correct rotation
            if solver == "bempp":
                exit_aperture.set_translation(translation, absolute=True)
                exit_aperture.set_rotation_angle_axis(angle=rotation[0], axis=X_AXIS, absolute=True)  # upright
                exit_aperture.set_rotation_angle_axis(angle=rotation[1], axis=Y_AXIS, absolute=False)  # match tilt
                exit_aperture.set_rotation_angle_axis(angle=rotation[2], axis=Z_AXIS, absolute=False)  # match exit
                translation = np.array([0.0, 0.0, 0.0])
                rotation = np.array([0.0, 0.0, 0.0])
            # Create geo string and load
            exit_aperture.create_geo_str(r=r, dz=dz, a=a, b=b, translation=translation, rotation=rotation,
                                         hole_type=hole_type, h=h, load=True, header=True)
            exit_aperture.color = "GREEN"

            assy.add_electrode(exit_aperture)

    if numerical_pars["make_cylinder"]:
        # Base cylinder parameters:
        r = numerical_pars["cylinder_params"]["radius"]
        zmin = numerical_pars["cylinder_params"]["zmin"]
        zmax = numerical_pars["cylinder_params"]["zmax"]
        voltage = numerical_pars["cylinder_params"]["voltage"]

        outer_cylinder = SICylinder(name="Outer Cylinder", voltage=voltage, offset=cylinder_offset)
        # translate = np.array([0.0, 0.0, 0,0])
        # outer_cylinder.set_translation(translate, absolute=True)
        # outer_cylinder.create_geo_str(r=r, dz=zmax - zmin, h=h, load=True, header=False)
        outer_cylinder.create_geo_str(r=r, zmin=zmin, zmax=zmax, h=0.0025, load=True, header=True)

        assy.add_electrode(outer_cylinder)

    if si.debug:
        assy.show(show_screen=True)

    numerical_vars["objects"] = assy

    si.analytic_parameters = analytic_pars
    si.analytic_variables = analytic_vars
    si.numerical_parameters = numerical_pars
    si.numerical_variables = numerical_vars

    return assy


def generate_meshed_model(si, apertures=None, cylinder=None):
    # TODO: Think about surface normals for electrodes and outer cylinder!

    generate_solid_assembly(si, apertures, cylinder)

    numerical_vars = si.numerical_variables

    if si._solver == "bempp":

        assert HAVE_BEMPP, "BEMPP not found. Aborting!"

        assy = numerical_vars["objects"]

        leaf_view = assy.get_bempp_mesh()

        numerical_vars["full mesh"] = {"verts": leaf_view["verts"],
                                       "elems": leaf_view["elems"],
                                       "domns": leaf_view["domns"]}

        if si.debug:
            _full_mesh = bempp.api.grid_from_element_data(leaf_view["verts"],
                                                          leaf_view["elems"],
                                                          leaf_view["domns"])
            _full_mesh.plot()

    elif si._solver == "fenics":

        assert HAVE_FENICS, "Fenics not found. Aborting!"
        assert HAVE_MESHIO, "Meshio not found. Aborting!"

        si.generate_vacuum_space()

        import fenics as fn
        import meshio
        import os

        os.system('gmsh -3 master_geometry.geo -format msh2 -v 0 -o master_geometry.msh')
        # os.system('dolfin-convert master_geometry.msh master_geometry.xml')
        msh = meshio.read('master_geometry.msh')
        meshio.write("master_geometry.xdmf", msh)

        meshio.write_points_cells("master_geometry_markers.xdmf",
                                  msh.points,
                                  {"tetra": msh.cells["tetra"]},
                                  cell_data={"tetra": {"gmsh:physical": msh.cell_data["tetra"]["gmsh:physical"]}})
        meshio.write_points_cells("master_geometry_boundaries.xdmf",
                                  msh.points,
                                  {"triangle": msh.cells["triangle"]},
                                  cell_data={"triangle": {"gmsh:physical": msh.cell_data["triangle"]["gmsh:physical"]}})

        mesh = fn.Mesh()
        fn.XDMFFile("master_geometry_markers.xdmf").read(mesh)

        markers = fn.MeshFunction("size_t", mesh, mesh.topology().dim())
        fn.XDMFFile("master_geometry_markers.xdmf").read(markers, "gmsh:physical")

        boundaries = fn.MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
        fn.XDMFFile("master_geometry_boundaries.xdmf").read(boundaries, "gmsh:physical")
        boundaries = fn.MeshFunction("size_t", mesh, boundaries)

        full_mesh = [mesh, markers, boundaries]

        numerical_vars["full_mesh"] = full_mesh

    si.numerical_variables = numerical_vars

    return numerical_vars["full mesh"]


def aperture_geometry_macro(si, fname=None):
    # File type for inventor is .ivb
    geo = si._variables_analytic["geo"] * 100.0  # Scaling for inventor
    trj = si._variables_analytic["trj_design"] * 100.0  # type: np.ndarray
    thickness = si._params_numerical["aperture_params"]["thickness"] * 100.0
    radius = si._params_numerical["aperture_params"]["radius"] * 100.0
    length = si._params_numerical["aperture_params"]["length"] * 100.0
    width = si._params_numerical["aperture_params"]["width"] * 100.0
    aperture_distance_top = si._params_numerical["aperture_params"]["top_distance"] * 100.0
    aperture_distance_bottom = si._params_numerical["aperture_params"]["bottom_distance"] * 100.0
    # voltage = si._params_numerical["aperture_params"]["voltage"]

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
    if fname is not None:
        import os
        with open(os.path.join(si._outp_folder, fname), "w") as outfile:
            outfile.writelines(aperture_string)

    return aperture_string


def electrode_geometry_macro(si, fname=None):
    # File type for inventor is .ivb
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
    if fname is not None:
        import os
        with open(os.path.join(si._outp_folder, fname), "w") as of:
            of.write(header_text + electrode_texts[0] + electrode_texts[1] + footer_text)

    return header_text + electrode_texts[0] + electrode_texts[1] + footer_text


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

    for name, electrode in si._variables_numerical["objects"].items():
        if si._debug:
            print("Working on {}".format(name))

        filename = os.path.join(folder, "{}.geo".format(name))
        with open(filename, 'w') as of:
            of.write(electrode["gmsh_str"])

    return 0
