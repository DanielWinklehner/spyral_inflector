from py_electrodes.py_electrodes import *  # From py_electrodes we also get HAVE_GMSH and GMSH_EXE
import matplotlib.pyplot as plt
from .vector import Vector
from PyPATools.particles import ParticleDistribution
import numpy as np

X_AXIS = np.array([1, 0, 0], float)
Y_AXIS = np.array([0, 1, 0], float)
Z_AXIS = np.array([0, 0, 1], float)


def _longitudinal_profile(specification, num_sections, *, scalar_mode, name):
    """Evaluate a geometry modifier on the normalized inflector length.

    The normalized coordinate ``xi`` runs from -1 at array index zero to +1
    at the final section. Supported specifications are:

    * scalar: retains the legacy behavior selected by ``scalar_mode``;
    * sequence: polynomial coefficients in ascending order, i.e.
      ``[c0, c1, c2]`` evaluates ``c0 + c1*xi + c2*xi**2``;
    * callable: called once as ``specification(xi)``;
    * dictionary with ``kind='polynomial'`` and ``coefficients=[...]``;
    * dictionary with ``kind='samples'``, ``values=[...]`` and optional
      ``positions=[...]``. Samples are linearly interpolated.

    ``positions`` use the same [-1, +1] normalized coordinate. A callable
    may return either a scalar or one value per section.
    """

    if num_sections < 1:
        raise ValueError(f"{name} requires at least one geometry section")

    xi = np.linspace(-1.0, 1.0, num_sections)

    if callable(specification):
        profile = specification(xi)

    elif isinstance(specification, dict):
        kind = str(specification.get("kind", "polynomial")).lower()

        if kind in ("polynomial", "poly"):
            coefficients = specification.get(
                "coefficients", specification.get("coeffs")
            )
            if coefficients is None:
                raise ValueError(
                    f"{name} polynomial specification requires "
                    "'coefficients'"
                )
            coefficients = np.asarray(coefficients, dtype=float)
            if coefficients.ndim != 1 or coefficients.size == 0:
                raise ValueError(
                    f"{name} polynomial coefficients must be a nonempty "
                    "one-dimensional sequence"
                )
            profile = np.polynomial.polynomial.polyval(xi, coefficients)

        elif kind in ("samples", "sampled", "control_points"):
            values = np.asarray(specification.get("values"), dtype=float)
            if values.ndim != 1 or values.size == 0:
                raise ValueError(
                    f"{name} sampled specification requires a nonempty "
                    "one-dimensional 'values' sequence"
                )

            positions = specification.get("positions")
            if positions is None:
                positions = np.linspace(-1.0, 1.0, values.size)
            else:
                positions = np.asarray(positions, dtype=float)

            if positions.shape != values.shape:
                raise ValueError(
                    f"{name} sampled 'positions' and 'values' must have "
                    "the same shape"
                )
            if np.any(np.diff(positions) <= 0.0):
                raise ValueError(
                    f"{name} sampled 'positions' must be strictly increasing"
                )
            if positions[0] > -1.0 or positions[-1] < 1.0:
                raise ValueError(
                    f"{name} sampled 'positions' must span [-1, +1]"
                )

            profile = np.interp(xi, positions, values)

        else:
            raise ValueError(
                f"Unsupported {name} profile kind {kind!r}; expected "
                "'polynomial' or 'samples'"
            )

    else:
        values = np.asarray(specification, dtype=float)

        if values.ndim == 0:
            value = float(values)
            if scalar_mode == "legacy_angling":
                # Preserve the old sentinel/disable behavior for zero and
                # negative scalar values, including the historical -1.
                if value <= 0.0:
                    profile = np.zeros(num_sections, dtype=float)
                else:
                    profile = np.linspace(value, -value, num_sections)
            elif scalar_mode == "constant":
                profile = np.full(num_sections, value, dtype=float)
            else:
                raise ValueError(f"Unknown scalar mode {scalar_mode!r}")

        elif values.ndim == 1 and values.size > 0:
            # A bare sequence is concise polynomial notation.
            profile = np.polynomial.polynomial.polyval(xi, values)
        else:
            raise ValueError(
                f"{name} must be a scalar, callable, dictionary, or "
                "nonempty one-dimensional polynomial coefficient sequence"
            )

    profile = np.asarray(profile, dtype=float)
    if profile.ndim == 0:
        profile = np.full(num_sections, float(profile), dtype=float)
    if profile.shape != (num_sections,):
        raise ValueError(
            f"{name} profile returned shape {profile.shape}; expected "
            f"({num_sections},)"
        )
    if not np.all(np.isfinite(profile)):
        raise ValueError(f"{name} profile contains non-finite values")

    return profile

HAVE_BEMPP = False
try:
    import bempp_cl.api
    from bempp_cl.api.shapes.shapes import __generate_grid_from_geo_string as generate_from_string
    from bempp_cl.api.grid import Grid as BemppGrid
    HAVE_BEMPP = True
except ImportError:
    bempp = None
    BemppGrid = None

HAVE_FENICS = False
try:
    import fenics as fn
    HAVE_FENICS = True
except ImportError:
    fn = None

HAVE_MESHIO = False
try:
    import meshio
    HAVE_MESHIO = True
except ImportError:
    meshio = None



class SIHyperbolicDipole(PyElectrode):
    def __init__(self, parent=None, name="New Dipole", voltage=0, offset=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this dipole
        self._offset = offset    


    def create_geo_str_old(self, r, dz, a, b, translation=None, rotation=None, h=0.005, load=True,header=True):

        offset = self._offset

        if translation is None:
            translation = np.array([0.0, 0.0, 0.0])

        if rotation is None:
            rotation = np.array([0.0, 0.0, 0.0])
        
        if header:
            geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)
        else:
            geo_str = ""

        #Define the point at which the ENTRANCE to the dipole is centered
        #Not the geometric center, it will be extruded only in the +z direction from here
        geo_str+="Point(1000) = {%f,%f,%f};\n"%(translation[0],translation[1],translation[2])

        #Generate the points that will define the hyperbolic curve
        geo_str+="N = 100;\n"
        geo_str+="umin = -1;\n"
        geo_str+="umax = 1;\n"
        geo_str+="x0 = 0;\n"
        geo_str+="For i In {0:N-1}\n"
        geo_str+="u = umin + (umax - umin) / (N - 1) * i;\n"
        geo_str+="x = x0 + %f * Sqrt(1 + u^2);\n"%(a)
        geo_str+="y = x0 + %f * u;\n"%(b)
        geo_str+="x2= -1.0*x;\n"
        geo_str+="y2= -1.0*y;\n"
        geo_str+="Point(1 + i) = {x, y, %f, 1.0};\n"%(translation[2])
        geo_str+="Point(N+1+i) = {x2,y2,%f, 1.0};\n"%(translation[2])
        geo_str+="EndFor\n"

        #Create a hyperbolic curve by splining the points above
        geo_str+="Spline(1) = {1:N};\n"
        
        #Create a circle centered at the origin and connecting the end pts of the hyperbola
        geo_str+="Circle(2) = {1,1000,N};\n"
        
        #Combine the hyperbolic and circular curves into a sincle wire element
        geo_str+="Wire(3) = {1,2};\n"

        #Generate a surface defined by the closed wire created above
        geo_str+="Plane Surface(4) = {3};\n"
        
        # Extrude the plane surface we just made along the z-axis by an amount dz.
        # Layers has to follow dz/h: it used to be Layers{1}, i.e. a single element
        # spanning the whole electrode length however fine h was. Recombine is left
        # off so the extruded walls come out as triangles -- BEM++ needs triangles,
        # and a mixed triangle/quadrangle mesh breaks the element assembly in
        # py_electrodes.PyElectrode.generate_mesh().
        n_layers = int(max(1, np.ceil(abs(dz) / h)))
        geo_str+="Ex[] = Extrude {0,0,%f} {Surface{4}; Layers{%d};};\n"%(dz, n_layers)

        #Rotate if needed
        geo_str+="Rotate {{0,0,1},{0,0,0},%f}{Volume{Ex[1]};}\n"%(rotation[2])
        
        if load:
            self.generate_from_geo_str(geo_str=geo_str)


        return geo_str

    def create_geo_str(
            self,
            r,
            dz,
            a,
            b,
            translation=None,
            rotation=None,
            h=0.005,
            load=True,
            header=True
    ):
        """
        Create a hyperbolic quadrupole electrode cross-section and extrude it.

        Parameters
        ----------
        r : float
            Outer radius of the electrode cross-section.
        dz : float
            Extrusion length along +z.
        a : float
            Hyperbola vertex radius. For a symmetric ideal quadrupole,
            this is the clear-aperture radius.
        b : float
            Hyperbola transverse scale. Use b == a for ideal quadrupole symmetry.
        translation : array-like, optional
            Electrode center position [x, y, z].
        rotation : array-like, optional
            Rotation angles [rx, ry, rz] in radians. Only rz is used here.
        h : float
            Maximum mesh characteristic length.
        load : bool
            Load the generated geometry into the PyElectrode object.
        header : bool
            Include the Gmsh OpenCASCADE header.
        """

        if translation is None:
            translation = np.array([0.0, 0.0, 0.0], dtype=float)
        else:
            translation = np.asarray(translation, dtype=float)

        if rotation is None:
            rotation = np.array([0.0, 0.0, 0.0], dtype=float)
        else:
            rotation = np.asarray(rotation, dtype=float)

        if translation.shape != (3,):
            raise ValueError("translation must contain exactly three values")

        if rotation.shape != (3,):
            raise ValueError("rotation must contain exactly three values")

        if a <= 0.0:
            raise ValueError("a must be positive")

        if b <= 0.0:
            raise ValueError("b must be positive")

        if r <= a:
            raise ValueError(
                f"Outer radius r must be greater than the hyperbola "
                f"vertex radius a; received r={r}, a={a}"
            )

        if dz == 0.0:
            raise ValueError("dz must be nonzero")

        # For:
        #
        #     x_local = a * sqrt(1 + u^2)
        #     y_local = b * u
        #
        # require x_local^2 + y_local^2 = r^2 at the endpoints.
        u_max = np.sqrt(
            (r**2 - a**2) / (a**2 + b**2)
        )

        tx, ty, tz = translation

        if header:
            geo_str = (
                'SetFactory("OpenCASCADE");\n'
                'Geometry.NumSubEdges = 100;\n'
                f'Mesh.CharacteristicLengthMax = {h:.16g};\n'
            )
        else:
            geo_str = ""

        # Center used by the circular rear boundary.
        geo_str += (
            f"Point(1000) = "
            f"{{{tx:.16g}, {ty:.16g}, {tz:.16g}, {h:.16g}}};\n"
        )

        # Hyperbolic face:
        #
        #     (x - tx)^2/a^2 - (y - ty)^2/b^2 = 1
        #
        geo_str += "N = 100;\n"
        geo_str += f"umin = {-u_max:.16g};\n"
        geo_str += f"umax = {u_max:.16g};\n"

        geo_str += "For i In {0:N-1}\n"
        geo_str += "    u = umin + (umax - umin) / (N - 1) * i;\n"
        geo_str += f"    x = {tx:.16g} + {a:.16g} * Sqrt(1 + u^2);\n"
        geo_str += f"    y = {ty:.16g} + {b:.16g} * u;\n"
        geo_str += (
            f"    Point(1 + i) = "
            f"{{x, y, {tz:.16g}, {h:.16g}}};\n"
        )
        geo_str += "EndFor\n"

        # Hyperbolic front face.
        geo_str += "Spline(1) = {1:N};\n"

        # Circular rear face centered at the translated quadrupole axis.
        # Point 1 and Point N are both exactly radius r from Point 1000.
        geo_str += "Circle(2) = {N, 1000, 1};\n"

        # Closed cross-sectional boundary and surface.
        geo_str += "Wire(3) = {1, 2};\n"
        geo_str += "Plane Surface(4) = {3};\n"

        # Extrude along +z. The layer count has to follow dz/h: Layers{1} puts a
        # single element across the whole electrode length however fine h is,
        # which left the quadrupole dipoles at a few dozen triangles.
        n_layers = int(max(1, np.ceil(abs(dz) / h)))
        geo_str += (
            f"Ex[] = Extrude {{0, 0, {dz:.16g}}} "
            f"{{Surface{{4}}; Layers{{{n_layers}}};}};\n"
        )

        # Rotate around the electrode's translated z-axis.
        if abs(rotation[2]) > 0.0:
            geo_str += (
                f"Rotate {{{{0, 0, 1}}, "
                f"{{{tx:.16g}, {ty:.16g}, {tz:.16g}}}, "
                f"{rotation[2]:.16g}}} "
                "{Volume{Ex[1]};}\n"
            )

        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        return geo_str

        
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
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
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
            geo_str += "Extrude {{ 0, 0, {} }} {{ Surface{{ {} }}; }}\n".format(2 * dz, 100 + offset)
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
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
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

    def create_geo_str(self, r, zmin, zmax, h=0.0075, load=True, header=True,offsetXY=[0.0,0.0]):
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
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)
        else:
            geo_str = ""

        geo_str += "// Cylinder\n"
        geo_str += "Cylinder({}) = {{ 0, 0, {}, 0, 0, {}, {}, 2 * Pi }};\n\n".format(1 + offset, zmin, zmax - zmin, r)

        
        if load:
            self.generate_from_geo_str(geo_str=geo_str)
        return geo_str


class SIElectrode(PyElectrode):
    # Match the intent of the old "last 10 of 100 slices" guard without
    # making the physical cut depth depend on the geometry resolution.
    _GAMMA_TERMINAL_FRACTION = 0.10

    def __init__(self, parent=None, name="New Spiral Electrode", voltage=10000, offset=0):
        super().__init__(name=name, voltage=voltage)
        self._parent = parent  # the spiral inflector that contains this aperture
        self._offset = offset

    def rotate_point_around_axis(self, point, axis, theta, axis_point=(0, 0, 0)):

        theta = np.radians(theta)

        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)

        point = np.array(point)
        axis_point = np.array(axis_point)
        translated_point = point - axis_point

        # Rodrigues' rotation formula (rotation matrix)
        W1, W2, W3 = axis
        x, y, z = translated_point

        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        one_minus_cos_t = 1 - cos_t

        Rot_mat = np.array([
            [cos_t + W1**2 * one_minus_cos_t,       W1*W2 * one_minus_cos_t - W3*sin_t, W1*W3 * one_minus_cos_t + W2*sin_t],
            [W2*W1 * one_minus_cos_t + W3*sin_t, cos_t + W2**2 * one_minus_cos_t,       W2*W3 * one_minus_cos_t - W1*sin_t],
            [W3*W1 * one_minus_cos_t - W2*sin_t, W3*W2 * one_minus_cos_t + W1*sin_t, cos_t + W3**2 * one_minus_cos_t]
        ])
        
        rotated_translated_point = Rot_mat @ np.array([x, y, z])
        rotated_point = rotated_translated_point + axis_point

        return rotated_point

    def _transformed_section_points(self, section_points, angle_degrees):
        """Apply the existing plate-angle rotation to one five-point slice."""
        points = np.asarray(section_points, dtype=float).copy()
        W1, W2, W3, W4 = points[:4]
        segment_center = 0.25 * (W1 + W2 + W3 + W4)
        section_normal = np.cross(W2 - W1, W3 - W1)
        normal_size = np.linalg.norm(section_normal)
        if normal_size <= 1.0e-14:
            raise ValueError("Cannot rotate a degenerate electrode section")
        section_normal /= normal_size

        if abs(angle_degrees) > 1.0e-14:
            points = np.array([
                self.rotate_point_around_axis(
                    point,
                    section_normal,
                    angle_degrees,
                    axis_point=segment_center,
                )
                for point in points
            ])

        return points

    def _gamma_exit_frame(self, raw_geo, angle_at_segment):
        """Return an orthonormal frame anchored to the upper exit edge."""
        upper_end = self._transformed_section_points(
            raw_geo[0:5, -1, :], angle_at_segment[-1]
        )
        lower_end = self._transformed_section_points(
            raw_geo[5:10, -1, :], angle_at_segment[-1]
        )

        outer_edge_1 = upper_end[2]
        outer_edge_2 = upper_end[3]
        plane_point = 0.5 * (outer_edge_1 + outer_edge_2)

        width_axis = outer_edge_2 - outer_edge_1
        width_size = np.linalg.norm(width_axis)
        if width_size <= 1.0e-14:
            raise ValueError(
                "Cannot construct gamma plane from a zero-width upper "
                "electrode exit edge"
            )
        width_axis /= width_size

        zero_normal = np.cross(width_axis, upper_end[4] - plane_point)
        zero_size = np.linalg.norm(zero_normal)
        if zero_size <= 1.0e-14:
            raise ValueError("Cannot determine the upper-electrode exit plane")
        zero_normal /= zero_size

        body_index = -2 if raw_geo.shape[1] > 1 else -1
        upper_body_section = self._transformed_section_points(
            raw_geo[0:5, body_index, :], angle_at_segment[body_index]
        )
        upper_body_center = np.mean(upper_body_section[:4], axis=0)
        if np.dot(zero_normal, upper_body_center - plane_point) < 0.0:
            zero_normal *= -1.0

        upper_center = np.mean(upper_end[:4], axis=0)
        lower_center = np.mean(lower_end[:4], axis=0)
        upper_direction = upper_center - lower_center

        ninety_normal = (
            upper_direction
            - np.dot(upper_direction, width_axis) * width_axis
            - np.dot(upper_direction, zero_normal) * zero_normal
        )
        ninety_size = np.linalg.norm(ninety_normal)
        if ninety_size <= 1.0e-14:
            ninety_normal = np.cross(zero_normal, width_axis)
            ninety_size = np.linalg.norm(ninety_normal)
        ninety_normal /= ninety_size
        if np.dot(ninety_normal, upper_direction) < 0.0:
            ninety_normal *= -1.0

        # Choose the width sign so (width, vertical, backward) is a
        # right-handed frame. Flipping this axis does not alter the cut plane.
        if np.dot(np.cross(width_axis, ninety_normal), zero_normal) < 0.0:
            width_axis *= -1.0

        return width_axis, ninety_normal, zero_normal, plane_point

    def _make_gamma_plane(self, raw_geo, angle_at_segment, gamma_angle):
        """Define the shared exit cut plane from the final angled geometry.

        The zero-degree normal points back into the electrode body, ensuring
        that the negative half-space contains no electrode at gamma=0. The
        90-degree normal points from the lower electrode toward the upper one.
        The Boolean cutter is separately bounded to the terminal region.
        """
        _, ninety_normal, zero_normal, plane_point = self._gamma_exit_frame(
            raw_geo,
            angle_at_segment,
        )

        gamma_radians = np.radians(gamma_angle)
        gamma_normal = (
            np.cos(gamma_radians) * zero_normal
            + np.sin(gamma_radians) * ninety_normal
        )
        gamma_normal /= np.linalg.norm(gamma_normal)

        return [
            gamma_normal,
            plane_point,
            np.dot(gamma_normal, plane_point),
        ]

    def _gamma_boolean_geo_str(
            self, electrode_volume, cutter_volume, raw_geo,
            angle_at_segment, gamma_normal, plane_point):
        """Return a planar gamma cut confined to the electrode terminus."""
        bounds = np.ptp(np.asarray(raw_geo).reshape(-1, 3), axis=0)
        cut_size = max(1.0, 10.0 * np.linalg.norm(bounds))

        geo_str = "\n// True planar gamma cut\n"
        geo_str += "Box({}) = {{ {}, {}, {}, {}, {}, {} }};\n".format(
            cutter_volume,
            -cut_size,
            -cut_size,
            -cut_size,
            2.0 * cut_size,
            2.0 * cut_size,
            cut_size,
        )

        z_axis = np.array([0.0, 0.0, 1.0])
        rotation_axis = np.cross(z_axis, gamma_normal)
        rotation_size = np.linalg.norm(rotation_axis)
        z_alignment = np.clip(np.dot(z_axis, gamma_normal), -1.0, 1.0)

        if rotation_size > 1.0e-14:
            rotation_axis /= rotation_size
            rotation_angle = np.arccos(z_alignment)
            geo_str += (
                "Rotate {{ {{ {}, {}, {} }}, {{ 0, 0, 0 }}, {} }} "
                "{{ Volume{{ {} }}; }}\n".format(
                    rotation_axis[0],
                    rotation_axis[1],
                    rotation_axis[2],
                    rotation_angle,
                    cutter_volume,
                )
            )
        elif z_alignment < 0.0:
            geo_str += (
                "Rotate { { 1, 0, 0 }, { 0, 0, 0 }, Pi } "
                "{ Volume{ " + str(cutter_volume) + " }; }\n"
            )

        geo_str += "Translate {{ {}, {}, {} }} {{ Volume{{ {} }}; }}\n".format(
            plane_point[0],
            plane_point[1],
            plane_point[2],
            cutter_volume,
        )

        # Bound the half-space in all three spatial directions around the
        # physical terminus. A slab bounded only along the exit normal is not
        # sufficient for a curved inflector: the distant entrance can curl
        # back into that same slab.
        upper_centers = np.mean(raw_geo[0:4, :, :], axis=0)
        segment_lengths = np.linalg.norm(
            np.diff(upper_centers, axis=0),
            axis=1,
        )
        path_length = np.sum(segment_lengths)
        if path_length <= 1.0e-14:
            raise ValueError(
                "Cannot bound the gamma cut on a zero-length electrode"
            )

        target_tail_length = self._GAMMA_TERMINAL_FRACTION * path_length
        tail_start = raw_geo.shape[1] - 1
        accumulated_length = 0.0
        while tail_start > 0 and accumulated_length < target_tail_length:
            tail_start -= 1
            accumulated_length += segment_lengths[tail_start]

        terminal_sections = []
        for section_index in range(tail_start, raw_geo.shape[1]):
            terminal_sections.append(self._transformed_section_points(
                raw_geo[0:5, section_index, :],
                angle_at_segment[section_index],
            ))
            terminal_sections.append(self._transformed_section_points(
                raw_geo[5:10, section_index, :],
                angle_at_segment[section_index],
            ))
        terminal_points = np.concatenate(terminal_sections, axis=0)

        terminal_min = np.min(terminal_points, axis=0)
        terminal_max = np.max(terminal_points, axis=0)
        terminal_span = terminal_max - terminal_min
        terminal_padding = max(
            1.0e-6,
            0.05 * np.linalg.norm(terminal_span),
            segment_lengths[-1] if segment_lengths.size else 0.0,
        )
        terminal_min -= terminal_padding
        terminal_max += terminal_padding
        terminal_volume = cutter_volume + 100_000

        geo_str += "Box({}) = {{ {}, {}, {}, {}, {}, {} }};\n".format(
            terminal_volume,
            terminal_min[0],
            terminal_min[1],
            terminal_min[2],
            terminal_max[0] - terminal_min[0],
            terminal_max[1] - terminal_min[1],
            terminal_max[2] - terminal_min[2],
        )

        geo_str += (
            "gamma_cutter_{}() = BooleanIntersection "
            "{{ Volume{{ {} }}; Delete; }}"
            "{{ Volume{{ {} }}; Delete; }};\n".format(
                electrode_volume,
                cutter_volume,
                terminal_volume,
            )
        )

        # Do not explicitly assign the result back to electrode_volume: OCC
        # creates the Boolean result before deleting its inputs, so reusing
        # that live tag raises "OpenCASCADE entity ... already exists". Let
        # Gmsh return the result tags and preserve the original tag when it can.
        geo_str += (
            "gamma_result_{}() = BooleanDifference "
            "{{ Volume{{ {} }}; Delete; }}"
            "{{ Volume{{ gamma_cutter_{}() }}; Delete; }};\n".format(
                electrode_volume,
                electrode_volume,
                electrode_volume,
            )
        )
        return geo_str

            
    def create_geo_str(self, raw_geo, elec_type, h=0.005,
                       load=True, header=True, gammaAng = -1.0,
                       gamma_plane=None, angling = -1,
                       vee_shape="linear"):

        f = self._offset

        if elec_type not in ["anode", "cathode"]:
            print("SIElectrode could not understand electrode type {}. Must be 'Anode' or 'Cathode'".format(type))
            return 1

        if header:
            geo_str = """SetFactory("OpenCASCADE");
Geometry.NumSubEdges = 100; // nicer display of curve
Mesh.CharacteristicLengthMax = {};  // maximum mesh size
""".format(h)
        else:
            geo_str = ""

        new_pt   = 1
        new_ln   = 1
        new_loop = 1
        new_vol  = 1

        # Shift the index in geo object for anode or cathode...
        if elec_type   == "anode":
            k = 0
        elif elec_type == "cathode":
            k = 5

        num_sections    = len(raw_geo[0, :, 0])

        vee_shape = str(vee_shape).strip().lower()
        if vee_shape in ("v", "vee", "legacy"):
            vee_shape = "linear"
        if vee_shape not in ("linear", "parabolic"):
            raise ValueError(
                f"Unsupported vee_shape {vee_shape!r}; expected "
                "'linear' or 'parabolic'"
            )

        #### ------- Set up the angling per slice ------ #####
        # A positive scalar retains the historical +angle -> -angle ramp.
        # Polynomial, sampled, and callable specifications may contain any
        # signed profile. This isn't flipped for anode vs cathode because the
        # cross products defining each slice normal are antiparallel.
        angle_at_segment = _longitudinal_profile(
            angling,
            num_sections,
            scalar_mode="legacy_angling",
            name="angling",
        )
        apply_angling = np.any(np.abs(angle_at_segment) > 1.0e-14)

        # Apply gamma as a Boolean half-space cut to the completed solid. This
        # makes the cut a genuine plane and avoids special-casing only the last
        # few loft sections.
        requested_gamma_angle = float(gammaAng)
        if requested_gamma_angle > 90.0:
            raise ValueError(
                "gammaAng must not exceed 90 degrees; received "
                f"{requested_gamma_angle}"
            )
        apply_boolean_gamma = requested_gamma_angle > 0.0

        if apply_boolean_gamma and elec_type == "anode":
            boolean_gamma_plane = self._make_gamma_plane(
                raw_geo,
                angle_at_segment,
                requested_gamma_angle,
            )
        elif apply_boolean_gamma:
            if gamma_plane is None or len(gamma_plane) != 3:
                raise ValueError(
                    "The anode must be processed first to define the shared "
                    "gamma cut plane"
                )
            boolean_gamma_plane = gamma_plane
        else:
            boolean_gamma_plane = [None, None, None]

        for j in range(num_sections):
            section_points = self._transformed_section_points(
                raw_geo[k:k + 5, j, :],
                angle_at_segment[j],
            )
            W1, W2, W3, W4, W5 = section_points

            for point in section_points:
                geo_str += "Point({}) = {{ {}, {}, {}, {} }};\n".format(
                    new_pt + f,
                    point[0],
                    point[1],
                    point[2],
                    h,
                )
                new_pt += 1

            # The legacy inner face consists of the two straight segments
            # W1--W5 and W5--W2.  For the parabolic option, sample the unique
            # parabola through those same edge and center points and ask Gmsh
            # to construct a smooth spline through the samples.  This changes
            # only the facing surface: sigma retains exactly the same meaning
            # as the edge-to-center sag, including a varying sigma(s).
            smooth_face_point_tags = None
            if vee_shape == "parabolic":
                face_edge_1, face_edge_2, face_center = W1, W2, W5
                edge_midpoint = 0.5 * (face_edge_1 + face_edge_2)
                half_span = 0.5 * (face_edge_2 - face_edge_1)
                sag_vector = edge_midpoint - face_center

                # W1, W5 and W2 keep their normal point tags. Additional
                # spline points use a distant, electrode-specific tag range
                # so the established five-point-per-section numbering and all
                # downstream offsets remain untouched.
                interior_eta = (-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)
                extra_tag_base = (
                    10_000_000
                    + (0 if elec_type == "anode" else 1_000_000)
                    + j * len(interior_eta)
                )
                extra_tags = []
                for local_index, eta in enumerate(interior_eta):
                    point = (
                        face_center
                        + eta * half_span
                        + (eta ** 2.0) * sag_vector
                    )
                    point_tag = extra_tag_base + local_index
                    extra_tags.append(point_tag)
                    geo_str += (
                        "Point({}) = {{ {}, {}, {}, {} }};\n".format(
                            point_tag,
                            point[0],
                            point[1],
                            point[2],
                            h,
                        )
                    )

                base_point = j * 5
                smooth_face_point_tags = [
                    base_point + 1 + f,
                    extra_tags[0],
                    extra_tags[1],
                    extra_tags[2],
                    base_point + 5 + f,
                    extra_tags[3],
                    extra_tags[4],
                    extra_tags[5],
                    base_point + 2 + f,
                ]
                    
                    
            # For each section, add the lines
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 0 + f, (j * 5) + 4 + f, (j * 5) + 2 + f)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 1 + f, (j * 5) + 3 + f, (j * 5) + 1 + f)
            geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 2 + f, (j * 5) + 3 + f, (j * 5) + 4 + f)

            if vee_shape == "linear":
                geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 3 + f, (j * 5) + 5 + f, (j * 5) + 2 + f)
                geo_str += "Line({}) = {{ {}, {} }};\n".format(new_ln + 4 + f, (j * 5) + 1 + f, (j * 5) + 5 + f)
            else:
                geo_str += "Spline({}) = {{ {} }};\n".format(
                    new_ln + 3 + f,
                    ", ".join(str(tag) for tag in smooth_face_point_tags),
                )

            new_ln += 5

            if vee_shape == "linear":
                geo_str += "Wire({}) = {{ {}, {}, {}, {}, {} }};\n\n".format(new_loop + f,
                                                                             (j * 5) + 3 + f,
                                                                             (j * 5) + 2 + f,
                                                                             (j * 5) + 5 + f,
                                                                             (j * 5) + 4 + f,
                                                                             (j * 5) + 1 + f)
            else:
                geo_str += "Wire({}) = {{ {}, {}, {}, {} }};\n\n".format(new_loop + f,
                                                                         (j * 5) + 3 + f,
                                                                         (j * 5) + 2 + f,
                                                                         (j * 5) + 4 + f,
                                                                         (j * 5) + 1 + f)

            new_loop += 1

        electrode_volume = new_vol + f
        geo_str += "Ruled ThruSections({}) = {{ {}:{} }};\n".format(
            electrode_volume,
            1 + f,
            new_loop - 1 + f,
        )

        if apply_boolean_gamma:
            gamma_normal = np.asarray(boolean_gamma_plane[0], dtype=float)
            plane_point = np.asarray(boolean_gamma_plane[1], dtype=float)
            cutter_volume = 9_000_000 + f
            geo_str += self._gamma_boolean_geo_str(
                electrode_volume,
                cutter_volume,
                raw_geo,
                angle_at_segment,
                gamma_normal,
                plane_point,
            )

        
        new_vol += 1

        # Call function in PyElectrode module we inherit from if load is not False
        if load:
            self.generate_from_geo_str(geo_str=geo_str)

        if elec_type == "cathode":
            return geo_str
        else:
            return geo_str, boolean_gamma_plane


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
            if np.isnan(_x):
                break
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
        self._tilt_angle = None
        self._face_angle = None

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

    @staticmethod
    def sort_points_by_angle(points):
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
        # The housing exit opening defaults to the entrance-aperture hole (length x width), rotated
        # with the exit electrodes; "exit_length" / "exit_width" override it so the exit can be
        # opened wider than the entrance slot without changing the entrance aperture.
        a = self._aperture_params.get("exit_length") or self._aperture_params["length"]
        b = self._aperture_params.get("exit_width") or self._aperture_params["width"]
        t_gap = self._aperture_params["top_distance"]
        b_gap = self._aperture_params["bottom_distance"]

        if span:
            zmin = np.min(geo[:, :, 2])
            zmax = np.max(geo[:, :, 2])

        norm_vec = Vector(trj[-1] - trj[-2]).normalized()

        translate = np.array([trj[-1][0] + norm_vec[0] * b_gap * 0.99,
                              trj[-1][1] + norm_vec[1] * b_gap * 0.99,
                              0.0])

        hole_type = self._aperture_params["hole_type"]

        if header:
            geo_str = """SetFactory("OpenCASCADE");
// Geometry.NumSubEdges = 100; // nicer display of curve
// Geometry.ToleranceBoolean = 1E-5;
// Geometry.Tolerance = 1E-10;
Mesh.CharacteristicLengthMax = {};  // maximum mesh size""".format(h)
        else:
            geo_str = "Geometry.ToleranceBoolean = 1E-5;\n"

        geo_str += "// Outside points\n"
        n_pts_out = np.shape(pts_out)[0]

        for i, pt in enumerate(pts_out):
            geo_str += "Point({}) = {{ {}, {}, 0, {} }};\n".format(i + offset, pt[0], pt[1], h)

        geo_str += "// Outside lines\n"
        n_lines_out = n_pts_out  # Should be the same number

        for i in range(n_lines_out - 1):
            geo_str += "Line({}) = {{ {}, {} }};\n".format(i + offset, i + offset, i + 1 + offset)
        geo_str += "Line({}) = {{ {}, {} }};\n".format(i + 1 + offset, n_pts_out - 1 + offset,
                                                       0 + offset)  # Connect last point to first point

        geo_str += "// Inside points\n"
        n_pts_in = np.shape(pts_in)[0]

        for i, pt in enumerate(pts_in):
            geo_str += "Point({}) = {{ {}, {}, 0, {} }};\n".format(i + n_pts_out + offset, pt[0], pt[1], h)

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

        geo_str += "// Tool to subtract\n"

        sub_tool = None  # This is either 5002 or disk_out[1], depending on the type of hole
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
            sub_tool = 2 + offset

        elif hole_type == "ellipse":
            geo_str += "Disk ({}) = {{ 0, 0, 0, {}, {} }};\n".format(500 + offset, 0.5 * a + 5E-9, 0.5 * b + 5E-9)
            geo_str += "disk_out[] = Extrude {{ 0, 0, {} }} {{ Surface{{ {} }}; }};\n".format(0.05, 500 + offset)
            sub_tool = "disk_out[1]"

        geo_str += "housing_out[] = Extrude {{ Surface{{ {} }}; }} Using Wire {{ {} }};\n".format(1 + offset,
                                                                                                  3 + offset)
        geo_str += "Translate {{ 0, 0, {} }} {{ Volume{{ housing_out[] }}; }}\n".format(zmin)

        geo_str += "Recursive Delete {{ Point{{ {}:{} }}; }}\n".format(point_index, point_index + num_wire_points - 1)
        geo_str += "Recursive Delete {{ Line{{ {}:{} }}; }}\n".format(line_index, line_index + num_wire_points - 2)
        geo_str += "Recursive Delete {{ Surface{{ {} }}; }}\n".format(1 + offset)

        geo_str += "Rotate {{ {{ 1.0, 0.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            np.pi / 2.0, sub_tool)
        geo_str += "Rotate {{ {{ 0.0, 1.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            self._tilt_angle, sub_tool)
        geo_str += "Rotate {{ {{ 0.0, 0.0, 1.0 }}, {{ 0.0, 0.0, 0.0 }}, {} }} {{ Volume{{ {} }}; }}\n".format(
            self._face_angle, sub_tool)

        geo_str += "Translate {{ {}, {}, {} }} {{ Volume{{ {} }}; }}\n".format(translate[0],
                                                                               translate[1],
                                                                               0.0,
                                                                               sub_tool)

        geo_str += "BooleanDifference({}) = {{ Volume {{ housing_out[] }}; Delete; }}{{ Volume {{ {} }}; Delete; }};\n".format(
            50 + offset, sub_tool)

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
    
    analytic_params = si.analytic_parameters
    analytic_vars = si.analytic_variables

    if analytic_vars["trj_design"] is None:
        print("No analytical design trajectory yet, generating...")
        si.generate_design_trajectory()

    print("Generating analytical geometry... ", end="")

    geos = []
    _ion = analytic_params["ion"]  # type: ParticleDistribution
    ns = analytic_params["ns"]  # type: int
    gap = analytic_params["gap"]  # type: float
    sigma = analytic_params["sigma"]  # type: float
    aspect_ratio = analytic_params["aspect_ratio"]  # type: float
    kp = analytic_vars["kp"]  # type: float
    b = analytic_vars["b"]  # type: np.ndarray
    cp = analytic_vars["c+"]  # type: float
    cm = analytic_vars["c-"]  # type: float
    trj_design = analytic_vars["trj_design"]  # type: np.ndarray

    for thickness in [0.0, analytic_params["dx"]]:

        end_distance = 2.0 * thickness + gap  # End to end distance of electrodes in (m)

        # Distance between electrodes at inflection angle theta
        d = end_distance * np.ones(ns) / (np.sqrt(1.0 + ((kp ** 2.0) * (np.sin(b)) ** 2.0)))

        # Save inner gap size vs deflection angle as class variable
        if thickness == 0.0:
            analytic_vars["d"] = d

        # x-component of the velocity vector
        vx = np.array(0.5 * _ion.v_mean_m_per_s * (np.sin(cp * b) + np.sin(cm * b)))

        # y-component of the velocity vector
        vy = np.array(-0.5 * _ion.v_mean_m_per_s * (np.cos(cp * b) - np.cos(cm * b)))

        # Rotation/flip
        if not ((analytic_vars["bf_design"] > 0.0) ^ (_ion.species.q > 0.0)):
            if si.debug:
                print("Flipping direction of cyclotron motion...", end="")
            vy = -vy

        v2 = np.sqrt((vx ** 2.0) + (vy ** 2.0))  # xy-magnitude of the velocity
        vz = np.sqrt((_ion.v_mean_m_per_s ** 2.0) - (v2 ** 2.0) + 0.0j)  # z-component of the velocity

        # Checks for imaginary components of the z-component
        for i in range(ns):
            if np.imag(vz[i]) != 0.0:
                vz[i] = 0.0 + 0.0j  # Redefines that element as 0

        vz = np.real(vz)
        v3 = np.sqrt((vx ** 2.0) + (vy ** 2.0) + (vz ** 2.0))  # 3-d magnitude of the velocity vector (Should = v)

        # Save vx, vy, vz in as class variable in same format as trj_design
        analytic_vars["v_design"] = np.array([vx, vy, vz]).T

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

        if analytic_params["rotation"] != 0.0:
            for i in range(analytic_params["ns"]):
                v_rh[0, i, :] = np.matmul(analytic_vars["rot"], v_rh[0, i, :])
                v_rh[1, i, :] = np.matmul(analytic_vars["rot"], v_rh[1, i, :])

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

    sigma_at_segment = _longitudinal_profile(
        sigma,
        ns,
        scalar_mode="constant",
        name="sigma",
    )
    sigma_displacement = diff_hat * sigma_at_segment[:, np.newaxis]

    geo = np.zeros([10, ns, 3])  # Initialize geo array

    # Upper Electrode
    geo[0, :, :] = geos[0][0, :, :] + sigma_displacement
    geo[1, :, :] = geos[0][1, :, :] + sigma_displacement
    geo[2, :, :] = geos[1][0, :, :] + sigma_displacement  # Should reduce the thickness
    geo[3, :, :] = geos[1][1, :, :] + sigma_displacement  # Should reduce the thickness
    geo[4, :, :] = 0.5 * (geos[0][0, :, :] + geos[0][1, :, :])

    # Lower Electrode
    geo[5, :, :] = geos[0][2, :, :] + sigma_displacement
    geo[6, :, :] = geos[0][3, :, :] + sigma_displacement
    geo[7, :, :] = geos[1][2, :, :]
    geo[8, :, :] = geos[1][3, :, :]
    geo[9, :, :] = 0.5 * (geos[0][2, :, :] + geos[0][3, :, :])

    analytic_vars["geo"] = geo

    print("Done!")

    return analytic_vars["geo"]


def generate_numerical_geometry(si):
    # This is a slightly modified version of the normal analytical method
    
    analytic_vars = si.analytic_variables
    analytic_params = si.analytic_parameters
    
    if analytic_vars["trj_design"] is None:
        print("No numerical design trajectory yet, generating...")
        si.generate_design_trajectory()

    print("Generating numerical geometry... ", end="")

    geos   = []
    _ion   = analytic_params["ion"]  # type: ParticleDistribution
    ns     = analytic_params["ns"]  # type: int
    gap    = analytic_params["gap"]  # type: float
    sigma  = analytic_params["sigma"]  # type: float
    aspect_ratio = analytic_params["aspect_ratio"]  # type: float
    kp = analytic_vars["kp"]  # type: float
    b = analytic_vars["b"]  # type: np.ndarray
    trj_design = analytic_vars["trj_design"]  # type: np.ndarray
    trj_vel = analytic_vars["trj_vel"]
    vx, vy, vz = trj_vel[:, 0], trj_vel[:, 1], trj_vel[:, 2]  # Unload the v components

    for thickness in [0.0, analytic_params["dx"]]:

        end_distance = 2.0 * thickness + gap  # End to end distance of electrodes in (m)

        # Distance between electrodes at inflection angle theta
        d = end_distance * np.ones(ns) / (np.sqrt(1.0 + ((kp ** 2.0) * (np.sin(b)) ** 2.0)))

        # Save inner gap size vs deflection angle as class variable
        if thickness == 0.0:
            analytic_vars["d"] = d

        v2 = np.sqrt((vx ** 2.0) + (vy ** 2.0))  # xy-magnitude of the velocity
        v3 = np.sqrt((vx ** 2.0) + (vy ** 2.0) + (vz ** 2.0))  # 3-d magnitude of the velocity vector (Should = v)

        # Save vx, vy, vz in as class variable in same format as trj_design
        analytic_vars["v_design"] = np.array([vx, vy, vz]).T

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
        #     for i in range(analytic_params["ns"]):
        #         v_rh[0, i, :] = np.matmul(xrot, v_rh[0, i, :])
        #         v_rh[1, i, :] = np.matmul(xrot, v_rh[1, i, :])
        #     # print("Applied a {:.4f} rad x rotation.".format(si._variables_optimization["x_rot"]))

        # if analytic_params["rotation"] != 0.0:
        #     for i in range(analytic_params["ns"]):
        #         v_rh[0, i, :] = np.matmul(analytic_vars["rot"], v_rh[0, i, :])
        #         v_rh[1, i, :] = np.matmul(analytic_vars["rot"], v_rh[1, i, :])

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

    sigma_at_segment = _longitudinal_profile(
        sigma,
        ns,
        scalar_mode="constant",
        name="sigma",
    )
    sigma_displacement = diff_hat * sigma_at_segment[:, np.newaxis]

    geo = np.zeros([10, ns, 3])  # Initialize geo array

    # Upper Electrode
    geo[0, :, :] = geos[0][0, :, :] + sigma_displacement
    geo[1, :, :] = geos[0][1, :, :] + sigma_displacement
    geo[2, :, :] = geos[1][0, :, :] + sigma_displacement  # Should reduce the thickness
    geo[3, :, :] = geos[1][1, :, :] + sigma_displacement  # Should reduce the thickness
    geo[4, :, :] = 0.5 * (geos[0][0, :, :] + geos[0][1, :, :])

    # Lower Electrode
    geo[5, :, :] = geos[0][2, :, :] + sigma_displacement
    geo[6, :, :] = geos[0][3, :, :] + sigma_displacement
    geo[7, :, :] = geos[1][2, :, :]
    geo[8, :, :] = geos[1][3, :, :]
    geo[9, :, :] = 0.5 * (geos[0][2, :, :] + geos[0][3, :, :])

    analytic_vars["geo"] = geo
    si.analytic_variables = analytic_vars

    print("Done!")

    return analytic_vars["geo"]


def get_norm_vec_and_angles_from_geo(geo):
    mid_vec_b = Vector(geo[8, -1, :] - geo[7, -1, :]).normalized()
    # tilt_angle is the angle of mid_vec_b with x/y plane

    tilt_angle = 0.5 * np.pi - mid_vec_b.angle_with(Vector(-Z_AXIS))

    # face angle is the angle of mid_vec_b projected into x/y plane with x/z plane
    temp_vec = Vector([mid_vec_b[0], mid_vec_b[1], 0.0])

    #face_angle = 0.5 * np.pi - temp_vec.angle_with(Vector(Y_AXIS))
    face_angle = np.arctan2(mid_vec_b[1],mid_vec_b[0])

    
    return tilt_angle, face_angle


def generate_vacuum_space(si):
    # TODO: Clean up
    assert si.numerical_parameters["make_cylinder"], "You need a cylinder/boundary to create the vacuum space!"

    numerical_vars = si.numerical_variables
    numerical_pars = si.numerical_parameters

    assy = numerical_vars["objects"]

    master_geo_str = "// Full .geo file for fenics mesh generation\n"

    master_geo_str += assy.get_electrode_by_name('SI Anode')._geo_str
    master_geo_str += assy.get_electrode_by_name('SI Cathode')._geo_str
    master_geo_str += assy.get_electrode_by_name('Outer Cylinder')._geo_str

    if numerical_pars["make_aperture"]:
        make_top_aperture = True
        master_geo_str += assy.get_electrode_by_name('Entrance Aperture')._geo_str
        if numerical_pars["make_housing"]:
            make_housing = True
            make_bottom_aperture = False
            master_geo_str += assy.get_electrode_by_name('Housing')._geo_str
        else:
            make_housing = False
            make_bottom_aperture = True
            master_geo_str += assy.get_electrode_by_name('Exit Aperture')._geo_str
    else:
        make_top_aperture = False
        make_housing = False
        make_bottom_aperture = False

    # for _, electrode in assy.electrodes.items():
    #     master_geo_str += electrode._geo_str

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
vacuum_boundary[] = Boundary { Volume{ 2001 }; };
N_vacuum = #vacuum_boundary[];
For j In {0:N_vacuum-1}
    Physical Surface (2000 + j) = { vacuum_boundary[j] };
EndFor
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

    master_geo_str += "Delete{ Volume{1, 1001, 2001"

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

// edge_list[] = Boundary { Volume{ 1 }; };
// N_edges = #edge_list[];
// Field[1] = Distance;
// Field[1].NNodesByEdge = 100;
// Field[1].EdgesList = {1:N_edges};

// Field[2] = Threshold;
// Field[2].IField = 1;
// Field[2].LcMin = 0.0015;
// Field[2].LcMax = 0.01;
// Field[2].DistMin = 0.003;
// Field[2].DistMax = 0.01;
// Background Field = 2;

// Mesh.OptimizeNetgen = 1;

"""

    with open(TEMP_DIR+'/master_geometry.geo', 'w') as outfile:
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

    geo        = analytic_vars["geo"]
    trj        = analytic_vars["trj_design"]
    # Design voltage times the operating scale (.get: objects saved before the key
    # existed). The design orbit and the electrodes stay those of the design voltage.
    voltage    = analytic_pars["volt"] * numerical_pars.get("volt_scale", 1.0)
    h          = numerical_pars["h"]
    gamma      = analytic_pars["gammaAng"] 
    anglingAng = analytic_pars["anglingAng"]
    vee_shape  = analytic_pars.get("vee_shape", "linear")
    
    # Variables for fenics solving, won't affect anything BEMPP related (ideally) -PW
    anode_offset = 0
    cathode_offset = 1000
    cylinder_offset = 2000
    entrance_offset = 3000
    exit_offset = 5000
    housing_offset = 5000

    anode = SIElectrode(name="SI Anode", voltage=voltage, offset=anode_offset)
    _, anode_gamma_plane = anode.create_geo_str(raw_geo=geo, elec_type="anode", h=h, load=True, header=True, gammaAng = gamma, angling = anglingAng, vee_shape=vee_shape)
    anode.color = "RED"

    cathode = SIElectrode(name="SI Cathode", voltage=-voltage, offset=cathode_offset)
    cathode.create_geo_str(raw_geo=geo, elec_type="cathode", h=h, load=True, header=True, gammaAng = gamma, gamma_plane = anode_gamma_plane, angling = anglingAng, vee_shape=vee_shape)
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
                                   h=h * 3,
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
            translation = np.array([0.0, 0.0, 0.0])
            rotation = np.array([0.0, 0.0, 0.0])

        entrance_aperture.create_geo_str(r=r, dz=dz, a=a, b=b, translation=translation, rotation=rotation,
                                         hole_type=hole_type,
                                         h=h, load=True, header=True)
        entrance_aperture.color = "BLACK"

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

            #norm_vec = Vector(trj[-1] - trj[-2]).normalized()

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
            exit_aperture.color = "BLACK"

            assy.add_electrode(exit_aperture)

    if numerical_pars["make_cylinder"]:
        # Base cylinder parameters:
        r = numerical_pars["cylinder_params"]["radius"]
        zmin = numerical_pars["cylinder_params"]["zmin"]
        zmax = numerical_pars["cylinder_params"]["zmax"]
        voltage = numerical_pars["cylinder_params"]["voltage"]
        
        outer_cylinder = SICylinder(name="Outer Cylinder", voltage=voltage, offset=cylinder_offset)
        outer_cylinder.create_geo_str(r=r, zmin=zmin, zmax=zmax, h=0.01, load=True, header=True,offsetXY=offsetXY)
        outer_cylinder.color = "GREEN"

        assy.add_electrode(outer_cylinder)


    # The quadrupoles and their grounded apertures are fixed in the lab frame, so
    # the entrance-centering shift must not move them. Everything created before
    # this point belongs to the inflector itself.
    _inflector_ids = set(assy.electrodes.keys())

    if numerical_pars["make_quadrupoles"]:
        a          = numerical_pars["quadrupole_params"]["a"]
        b          = numerical_pars["quadrupole_params"]["b"]
        r          = numerical_pars["quadrupole_params"]["radius"]
        z_starts   = numerical_pars["quadrupole_params"]["z_starts"]
        quad_lens  = numerical_pars["quadrupole_params"]["lengths"]
        quad_volts = numerical_pars["quadrupole_params"]["voltages"]
        aper_rad   = numerical_pars["quadrupole_params"]["aper_rad"]

        pi     = 3.14159265358
        # Grounded aperture plates at both ends of each quadrupole: thickness and the axial gap
        # between a plate and the pole ends. The old fixed 1 mm gap put 10 kV across 1 mm and
        # shorted the fringe field; both are now quadrupole_params entries, defaults unchanged.
        aper_t = numerical_pars["quadrupole_params"].get("plate_thickness", 0.005)
        gap    = numerical_pars["quadrupole_params"].get("plate_gap", 0.001)
        # shared_plates: one grounded plate between consecutive quadrupoles instead of an exit
        # plate of one and an entrance plate of the next. The caller places quad k at
        # z_starts[k-1] + lengths[k-1] + 2 * gap + plate_thickness so the shared plate (the
        # previous quad's "ext" plate) sits one gap in front of it.
        shared = numerical_pars["quadrupole_params"].get("shared_plates", False)

        # Mesh size for the quadrupole electrodes. These used to be hardcoded
        # (0.005 for the apertures, 0.01 for the dipoles) and so ignored the "h"
        # parameter set on the SpiralInflector.
        quad_h = numerical_pars["h"]
        
        nquads = len(quad_volts)

        for pole in range(nquads):

            if pole == 0 or not shared:
                A1 = SIAperture(name="ent%i"%(4*pole),voltage=0.0)
                A1.create_geo_str(r=r, dz=aper_t, a=aper_rad, b=aper_rad, translation=[0,0,z_starts[pole]-aper_t/2.0-gap], hole_type="ellipse", h=quad_h, load=True,header=True)
                A1.color="BLACK"
                assy.add_electrode(A1)
            
            A2 = SIAperture(name="ext%i"%(4*pole),voltage=0.0)
            A2.create_geo_str(r=r, dz=aper_t, a=aper_rad, b=aper_rad, translation=[0,0,z_starts[pole]+quad_lens[pole]+aper_t/2.0+gap], hole_type="ellipse", h=quad_h, load=True,header=True)
            A2.color="BLACK"
            assy.add_electrode(A2)

            
            D1 = SIHyperbolicDipole(name="D%i"%(4*pole), voltage=quad_volts[pole])
            D1.create_geo_str(r=r,dz=quad_lens[pole],a=a,b=b,h=quad_h,translation=[0,0,z_starts[pole]],rotation=[0.0,0.0,0.0],load=True,header=True)
            D1.color="BLUE"            
            assy.add_electrode(D1)
            
            D2 = SIHyperbolicDipole(name="D%i"%(4*pole+1), voltage=quad_volts[pole])
            D2.create_geo_str(r=r,dz=quad_lens[pole],a=a,b=b,h=quad_h,translation=[0,0,z_starts[pole]],rotation=[0.0,0.0,pi],load=True,header=True)
            D2.color="BLUE"
            assy.add_electrode(D2)
         
            D3 = SIHyperbolicDipole(name="D%i"%(4*pole+2), voltage=-1.0*quad_volts[pole])
            D3.create_geo_str(r=r,dz=quad_lens[pole],a=a,b=b,h=quad_h,translation=[0,0,z_starts[pole]],rotation=[0.0,0.0,pi/2],load=True,header=True)
            D3.color="RED"
            assy.add_electrode(D3)
            
            D4 = SIHyperbolicDipole(name="D%i"%(4*pole+3), voltage=-1.0*quad_volts[pole])
            D4.create_geo_str(r=r,dz=quad_lens[pole],a=a,b=b,h=quad_h,translation=[0,0,z_starts[pole]],rotation=[0.0,0.0,3*pi/2],load=True,header=True)
            D4.color="RED"
            assy.add_electrode(D4)
        

    if si.debug:
        assy.show(show_screen=True)

        
    # Centering shift from optimize_fringe(apply_shift=True). Applied here so it
    # survives assembly regeneration; get_bempp_mesh() then bakes it into the BEM
    # mesh, and the STEP export picks it up too.
    _shift = si.track_variables.get("shift_applied")
    _shift_lab = si.track_variables.get("shift_applied_lab")

    if _shift is not None:
        _shift = np.asarray(_shift, dtype=float)
        # Lab-frame hardware (quadrupoles and their grounded apertures) takes only the
        # axial component, so it stays centred on the incoming beam while preserving
        # its spacing to the inflector.
        _shift_lab = (np.zeros(3) if _shift_lab is None
                      else np.asarray(_shift_lab, dtype=float))

        # Added to, not replacing, the translation an electrode already carries: the
        # entrance (and exit) aperture is positioned through set_translation above,
        # and an absolute shift here put it at the origin instead of just below the
        # electrode entrance. The electrodes are rebuilt on every call, so the shift
        # is added exactly once.
        for _eid, _elec in assy.electrodes.items():
            _elec.set_translation(_shift if _eid in _inflector_ids else _shift_lab,
                                  absolute=False)

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

    if si.solver == "bempp":

        assert HAVE_BEMPP, "BEMPP not found. Aborting!"

        assy = numerical_vars["objects"]

        leaf_view = assy.get_bempp_mesh()

        # print("len verts:", len(leaf_view["verts"]),
        #       "len elems:", len(leaf_view["elems"]),
        #       "len domns:", len(leaf_view["domns"]))
        # print()
        # print(leaf_view["verts"])
        # print(leaf_view["elems"])
        # print(leaf_view["domns"])

        numerical_vars["full mesh"] = {"verts": leaf_view["verts"],
                                       "elems": leaf_view["elems"],
                                       "domns": leaf_view["domns"]}

        if si.debug:
            _full_mesh = BemppGrid(leaf_view["verts"],
                                   leaf_view["elems"],
                                   leaf_view["domns"])

            # bempp.api.PLOT_BACKEND = "gmsh"
            # _full_mesh.plot()

    elif si.solver == "fenics":

        assert HAVE_FENICS, "Fenics not found. Aborting!"
        assert HAVE_MESHIO, "Meshio not found. Aborting!"
        assert HAVE_GMSH, "Gmsh was not properly installed with py_electrodes. Aborting!"

        si.generate_vacuum_space()

        import os

        os.system(GMSH_EXE+' -3 '+TEMP_DIR+'/master_geometry.geo -format msh2 -v 0 -o '+TEMP_DIR+'/master_geometry.msh')
        # os.system('dolfin-convert master_geometry.msh master_geometry.xml')
        msh = meshio.read(TEMP_DIR + "/master_geometry.msh")
        meshio.write(TEMP_DIR + "/master_geometry.xdmf", msh)

        meshio.write_points_cells(TEMP_DIR + "/master_geometry_markers.xdmf",
                                  msh.points,
                                  {"tetra": msh.cells["tetra"]},
                                  cell_data={"tetra": {"gmsh:physical": msh.cell_data["tetra"]["gmsh:physical"]}})
        meshio.write_points_cells(TEMP_DIR + "/master_geometry_boundaries.xdmf",
                                  msh.points,
                                  {"triangle": msh.cells["triangle"]},
                                  cell_data={"triangle": {"gmsh:physical": msh.cell_data["triangle"]["gmsh:physical"]}})

        mesh = fn.Mesh()
        fn.XDMFFile(TEMP_DIR + "/master_geometry_markers.xdmf").read(mesh)

        markers = fn.MeshFunction("size_t", mesh, mesh.topology().dim())
        fn.XDMFFile(TEMP_DIR + "/master_geometry_markers.xdmf").read(markers, "gmsh:physical")

        boundaries = fn.MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
        fn.XDMFFile(TEMP_DIR + "/master_geometry_boundaries.xdmf").read(boundaries, "gmsh:physical")
        boundaries = fn.MeshFunction("size_t", mesh, boundaries)

        full_mesh = [mesh, markers, boundaries]

        numerical_vars["full_mesh"] = full_mesh

    si.numerical_variables = numerical_vars

    return numerical_vars["full mesh"]


def aperture_geometry_macro(si, fname=None):
    # File type for inventor is .ivb
    
    analytic_vars = si.analytic_variables
    analytic_params = si.analytic_parameters
    numerical_params = si.numerical_parameters
    
    geo = analytic_vars["geo"] * 100.0  # Scaling for inventor
    trj = analytic_vars["trj_design"] * 100.0  # type: np.ndarray
    thickness = numerical_params["aperture_params"]["thickness"] * 100.0
    radius = numerical_params["aperture_params"]["radius"] * 100.0
    length = numerical_params["aperture_params"]["length"] * 100.0
    width = numerical_params["aperture_params"]["width"] * 100.0
    aperture_distance_top = numerical_params["aperture_params"]["top_distance"] * 100.0
    aperture_distance_bottom = numerical_params["aperture_params"]["bottom_distance"] * 100.0
    # voltage = numerical_params["aperture_params"]["voltage"]

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
    analytic_vars = si.analytic_variables
    analytic_params = si.analytic_parameters
    
    geo = analytic_vars["geo"] * 100.0  # Fix scaling in inventor

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

            for i in range(analytic_params["ns"]):
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
    from .tk_filedialog import FileDialog
    if filename is None:
        fd = FileDialog()
        folder, _ = fd.get_filename("folder")
        if folder is None:
            return 0
    else:
        folder = os.path.split(filename)[0]

    for name, electrode in si._variables_numerical["objects"].items():
        if si.debug:
            print("Working on {}".format(name))

        filename = os.path.join(folder, "{}.geo".format(name))
        with open(filename, 'w') as of:
            of.write(electrode["gmsh_str"])

    return 0
