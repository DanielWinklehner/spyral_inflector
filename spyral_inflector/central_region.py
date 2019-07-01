import numpy as np
from dans_pymodules import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from .field_solving import *
from .geometry import *
from .optimization import *
from .plotting import *
from .trajectories import *

colors = MyColors()

gmsh_macros = """

Macro MakeDeeSurface

    // Make the n'th dee
    // Assume there's a variable start_idx, which defines the first point
    // and another variable end_idx, which defines the last
    
    first_line_idx = newl; Line(first_line_idx) = {start_idx, start_idx + 1};
    
    For k In {1:end_idx-start_idx-1}
        Line (newl) = {start_idx + k, start_idx + k + 1};
    EndFor
    
    // Make the last line, close the loop with a wire, and create a surface
    last_line_idx = newl; Line (last_line_idx) = {start_idx, end_idx};
    wire_idx = newll; Wire (wire_idx) = { first_line_idx:last_line_idx };
    surface_idx = news; Surface (surface_idx) = { wire_idx };

Return

Macro MakeDeeVolume
    // Needs surface_idx, thickness, and gap defined
    
    // Extrude the top part of the dee
    top_extrude[] = Extrude {0.0, 0.0, thickness} { Surface { surface_idx }; };
    dee_top_vol = top_extrude[1];  // Volume tag stored in [1]
    
    // Extrude the bottom part of the dee
    bottom_extrude[] = Extrude {0.0, 0.0, -thickness} { Surface { surface_idx }; };
    dee_bottom_vol = bottom_extrude[1];
    
    // Translate to create a gap
    Translate {0.0, 0.0, gap/2} { Volume { dee_top_vol }; }
    Translate {0.0, 0.0, -gap/2} { Volume { dee_bottom_vol }; }
    
Return

Macro MakePost

    // Make an extrusion based on 4 points and boolean union with dee vols
    // Point tags are p1, p2, p3, p4
    
    first_line_idx = newl; Line (first_line_idx) = { p1, p2 };
    Line (newl) = { p2, p3 };
    Line (newl) = { p3, p4 };
    last_line_idx = newl; Line (last_line_idx) = { p4, p1 };
    
    wire_idx = newll; Wire (wire_idx) = { first_line_idx:last_line_idx };
    surface_idx = news; Surface (surface_idx) = { wire_idx };
    
    height = gap + 2*thickness;
    post_extrude[] = Extrude { 0.0, 0.0, height } { Surface { surface_idx }; };
    post_vol = post_extrude[1];
    
    Translate { 0.0, 0.0, -height/2 } { Volume { post_vol }; }
    dee_vol = newv;
    BooleanUnion(dee_vol) = { Volume { dee_top_vol, dee_bottom_vol }; Delete; } { Volume { post_vol }; Delete; };
    
Return
"""


class CentralRegion(PyElectrodeAssembly):
    def __init__(self, spiral_inflector=None,
                 r_cr=0.3,
                 dee_voltage=70e3,
                 dee_opening_angle=30.0,
                 rf_phase=0.0,
                 rf_freq=32.8e6,
                 harmonic=4,
                 ion=None):
        super().__init__(name="Central Region")

        self._debug = False
        self._initialized = False

        self._spiral_inflector = None

        self._params_analytic = {}
        self._variables_analytic = {}
        self._params_numerical = {}
        self._variables_numerical = {}
        self._variables_track = {}

        self._abstract_dees = []
        self._dees = []
        self._dummy_dees = []

        self._tracked_trjs = []
        self._track_trjs_vel = []

        self.rf_phase = rf_phase
        self.rf_freq = rf_freq
        self.harmonic = harmonic
        self.dee_voltage = dee_voltage
        self.r_cr = r_cr  # Radius at which we no longer consider the central region as a different entity
        self.dee_opening_angle = dee_opening_angle

        self._dee_z_func = None

        if spiral_inflector is not None:
            self.set_inflector(spiral_inflector)

        if ion is not None:
            self.analytic_parameters["ion"] = ion

        self._xi, self._vi = None, None

    @property
    def debug(self):
        return self._debug

    @property
    def analytic_parameters(self):
        return self._params_analytic

    @property
    def analytic_variables(self):
        return self._variables_analytic

    @property
    def numerical_parameters(self):
        return self._params_numerical

    @property
    def numerical_variables(self):
        return self._variables_numerical

    @property
    def track_variables(self):
        return self._variables_track

    @analytic_parameters.setter
    def analytic_parameters(self, analytic_parameters):
        self._params_analytic = analytic_parameters

    @analytic_variables.setter
    def analytic_variables(self, analytic_variables):
        self._variables_analytic = analytic_variables

    @numerical_parameters.setter
    def numerical_parameters(self, numerical_parameters):
        self._params_numerical = numerical_parameters

    @numerical_variables.setter
    def numerical_variables(self, numerical_variables):
        self._variables_numerical = numerical_variables

    @track_variables.setter
    def track_variables(self, track_variables):
        self._variables_track = track_variables

    # def dee_crossing(self, rc, rd, v, t):
    #     ion = self._params_analytic["ion"]
    #     initial_v = ion.v_m_per_s()
    #
    #     dE = 0.0
    #     ttf = 1.0
    #
    #     for segment in self._segments:
    #         if segment.check_intersection(rc[:2], rd[:2]):
    #             print('Intersection')
    #             dE += ion.q() * np.sin(2 * np.pi * self.rf_freq * t + self.rf_phase + segment.phase) * self.v_peak * ttf
    #             print(np.mod(np.rad2deg(2 * np.pi * self.rf_freq * t + self.rf_phase), 360))
    #
    #     # if not intersected:
    #     #     return v
    #     # dE = ion.q() * np.sin( 2 * np.pi * self.rf_freq * t + self.rf_phase) * self.v_peak * ttf
    #     # print("Gained energy {:.3f} keV".format(dE / 1E3))
    #     ion.calculate_from_energy_mev(ion._energy_mev + dE / 1E6)
    #     return ion.v_m_per_s() * v / initial_v

    def initialize(self, xi=None, vi=None, dee_z_func=None):
        if xi is None and vi is None:
            trj_design = self._spiral_inflector.analytic_variables["trj_design"]
            v_design = self._spiral_inflector.analytic_variables["v_design"]
            self._xi, self._vi = trj_design[-1, :], v_design[-1, :]
        else:
            self._xi, self._vi = xi, vi

        if dee_z_func is not None:
            self._dee_z_func = dee_z_func
        else:
            def default_z_func(r):
                return 0.0

            self._dee_z_func = default_z_func

    def set_inflector(self, spiral_inflector):
        self._spiral_inflector = spiral_inflector
        self._params_analytic["ion"] = self._spiral_inflector.analytic_parameters["ion"]
        self._params_analytic["bf_itp"] = self._spiral_inflector.analytic_parameters["bf_itp"]

        for electrode in spiral_inflector.numerical_variables["objects"].electrodes.values():
            self.add_electrode(electrode)


    def load_bfield(self, bfield=None, bf_scale=1.0, spatial_unit="cm", **kwargs):

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

            self._params_analytic["bf_itp"].load_field_from_file(filename=bfield, **kwargs)

            print("Successfully loaded B-Field from file")

        if self._initialized:
            self.initialize()

    def add_dee(self, dee):

        dee.generate_geometry()
        self._dees.append(dee)
        self.add_electrode(dee)

        return 0

    def split_dees(self):

        for dee in self._dees:
            dee.split()

        return 0

    def plot_bfield(self, nx=100, ny=100):
        from matplotlib import cm

        x, y = np.linspace(-0.1, 0.1, nx), np.linspace(-0.1, 0.1, ny)
        B = np.zeros([nx, ny])
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                B[i, j] = self._params_analytic["bf_itp"](np.array([xi, yi, 0.0]))[2]

        fig = plt.figure()
        xx, yy = np.meshgrid(x, y)
        ax = fig.add_subplot(111)
        sc = ax.contourf(xx, yy, B, levels=30, cmap=cm.coolwarm)
        ax.set_aspect(1)
        fig.colorbar(sc)
        # ax.plot_surface(xx, yy, B, cmap=cm.coolwarm)
        plt.show()

    def track(self, **kwargs):
        r, v = central_region_track(self, r_start=self._xi, v_start=self._vi, dt=1e-11, **kwargs)
        self._tracked_trjs.append(r)
        self._track_trjs_vel.append(v)
        return r, v

    def track_segment(self, **kwargs):

        r, v = track_segment(self, r_start=self._xi, v_start=self._vi, dt=1e-11, **kwargs)
        self._tracked_trjs.append(r)
        self._track_trjs_vel.append(v)

        return r, v

    def make_dees(self, n, gap, thickness, cl=0.1, **kwargs):

        # cl = 0.1
        my_dee = AbstractDee(opening_angle=self.dee_opening_angle, char_len=cl,
                             gap=gap, thickness=thickness)
        my_dee.initialize()

        my_second_dee = AbstractDee(opening_angle=self.dee_opening_angle, char_len=cl,
                                    gap=gap, thickness=thickness)
        my_second_dee.initialize()
        my_second_dee.rotate(90, angle_unit="deg")

        my_third_dee = AbstractDee(opening_angle=self.dee_opening_angle, char_len=cl,
                                   gap=gap, thickness=thickness)
        my_third_dee.initialize()
        my_third_dee.rotate(180, angle_unit="deg")

        my_fourth_dee = AbstractDee(opening_angle=self.dee_opening_angle, char_len=cl,
                                    gap=gap, thickness=thickness)
        my_fourth_dee.initialize()
        my_fourth_dee.rotate(270, angle_unit="deg")

        self._abstract_dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]

        for abs_dee in self._abstract_dees:
            for _ in range(n):
                abs_dee.next_top_segment(angle_offset=0)
                abs_dee.next_bottom_segment(angle_offset=0)
            abs_dee.split()
            dee = Dee(parent=abs_dee,
                      gap=gap,
                      thickness=thickness,
                      **kwargs)
            dee.generate_geometry()
            self.add_electrode(dee)
            self._dees.append(dee)

    def make_dummy_dees(self, gap, thickness, **kwargs):

        for i in range(len(self._abstract_dees) - 1):
            dummy_dee = DummyDee(self._abstract_dees[i], self._abstract_dees[i + 1],
                                 gap=gap,
                                 thickness=thickness,
                                 **kwargs)
            dummy_dee.generate_geometry()
            self.add_electrode(dummy_dee)
            self._dummy_dees.append(dummy_dee)
        dummy_dee = DummyDee(self._abstract_dees[-1], self._abstract_dees[0],
                             gap=gap,
                             thickness=thickness,
                             **kwargs)
        dummy_dee.generate_geometry()
        self.add_electrode(dummy_dee)
        self._dummy_dees.append(dummy_dee)

    def solve_bempp(self):
        self.numerical_variables["objects"] = self  # TODO: Make the solve_bempp more generic for electrode assms
        leaf_view = self.get_bempp_mesh()

        self.numerical_variables["full mesh"] = {"verts": leaf_view["verts"],
                                                 "elems": leaf_view["elems"],
                                                 "domns": leaf_view["domns"]}
        self.numerical_parameters["gmres_tol"] = 0.0001
        solve_bempp(self)

    def plot_dees(self, ax=None, show=False):
        for abs_dee in self._abstract_dees:
            abs_dee.plot_segments(show=show, ax=ax)

class AbstractDee(PyElectrode):
    def __init__(self,
                 char_len=0.03,
                 opening_angle=30.0,
                 gap=0.05,
                 thickness=0.025,
                 **kwargs):
        super().__init__(name="Dee", **kwargs)

        self.opening_angle = opening_angle
        self._opening_angle = np.deg2rad(self.opening_angle)
        self.angle = 0.0
        self._angle = 0.0

        self._top_segments, self._bottom_segments = [], []
        self._char_len = char_len  # Characteristic length of CRSegments
        self._gap = gap
        self._thickness = thickness

        self._is_split = False

        # Assume dee is centered on x-axis, "top" side is in +y, "bottom" side is in -y
        self.top_dummy_dee_segs = []
        self.top_dee_segs = []
        self.bottom_dummy_dee_segs = []
        self.bottom_dee_segs = []

        self.top_gap_vec = None
        self.bottom_gap_vec = None

        self._debug = True
        self.zfun = None

    def rotate(self, angle, angle_unit="deg"):
        # TODO: This should be done using the PyElectrodes methods
        """
        Rotates the entire Dee object by angle
        :param angle:
        :return:
        """
        if angle_unit == "deg":
            self.angle = angle
            self._angle = np.deg2rad(angle)
        elif angle_unit == "rad":
            self.angle = np.rad2deg(angle)
            self._angle = angle
        else:
            return 1

        for seg in self._top_segments + self._bottom_segments:
            seg.rotate(angle, angle_unit=angle_unit)

        return 0

    def initialize(self):

        # Initial top segment

        ra = self._char_len * np.array([np.cos(self._opening_angle / 2.0),
                                        np.sin(self._opening_angle / 2.0),
                                        0.0])
        rb = 2 * self._char_len * np.array([np.cos(self._opening_angle / 2.0),
                                            np.sin(self._opening_angle / 2.0),
                                            0.0])
        top_seg = CRSegment(ra, rb)

        # Initial bottom segment
        ra = self._char_len * np.array([np.cos(-self._opening_angle / 2.0),
                                        np.sin(-self._opening_angle / 2.0),
                                        0.0])
        rb = 2 * self._char_len * np.array([np.cos(-self._opening_angle / 2.0),
                                            np.sin(-self._opening_angle / 2.0),
                                            0.0])
        bottom_seg = CRSegment(ra, rb)

        self._top_segments.append(top_seg)
        self._bottom_segments.append(bottom_seg)

        pass

    def next_bottom_segment(self, angle_offset=0.0, length=None, angle_unit="deg"):

        if angle_unit == "deg":
            _angle = np.deg2rad(angle_offset)
        else:
            _angle = angle_offset

        if length is None:
            length = self._char_len

        ra = np.array([0.0, 0.0, 0.0])
        rb = length * np.array([np.cos(_angle), np.sin(_angle), 0.0])

        prev_seg = self._bottom_segments[-1]
        next_seg = CRSegment(ra, rb, color=1)

        next_seg.rotate(self.angle - self.opening_angle / 2.0, angle_unit="deg")
        next_seg.translate(prev_seg.rb)

        self._bottom_segments.append(next_seg)

        return 0

    def next_top_segment(self, angle_offset=0.0, length=None, angle_unit="deg"):

        if angle_unit == "deg":
            _angle = np.deg2rad(angle_offset)
        else:
            _angle = angle_offset

        if length is None:
            length = self._char_len

        ra = np.array([0.0, 0.0, 0.0])
        rb = length * np.array([np.cos(_angle), np.sin(_angle), 0.0])

        prev_seg = self._top_segments[-1]
        next_seg = CRSegment(ra, rb, color=1)

        next_seg.rotate(self.angle + self.opening_angle / 2.0, angle_unit="deg")
        next_seg.translate(prev_seg.rb)

        self._top_segments.append(next_seg)

        return 0

    def split(self, gap=0.012):

        # Split top
        for mid_seg in self._top_segments:
            theta = self._opening_angle / 2.0 + np.pi / 2.0 + self._angle
            gap_vec = (gap / 2.0) * np.array([np.cos(theta), np.sin(theta), 0.0])

            dd_seg = CRSegment(ra=mid_seg.ra + gap_vec, rb=mid_seg.rb + gap_vec, color=2)
            d_seg = CRSegment(ra=mid_seg.ra - gap_vec, rb=mid_seg.rb - gap_vec, color=1)

            self.top_dummy_dee_segs.append(dd_seg)
            self.top_dee_segs.append(d_seg)
            self.top_gap_vec = gap_vec

        # Split bottom
        for mid_seg in self._bottom_segments:
            theta = -self._opening_angle / 2.0 + np.pi / 2.0 + self._angle
            gap_vec = (gap / 2.0) * np.array([np.cos(theta), np.sin(theta), 0.0])

            dd_seg = CRSegment(ra=mid_seg.ra - gap_vec, rb=mid_seg.rb - gap_vec, color=2)
            d_seg = CRSegment(ra=mid_seg.ra + gap_vec, rb=mid_seg.rb + gap_vec, color=1)

            self.bottom_dummy_dee_segs.append(dd_seg)
            self.bottom_dee_segs.append(d_seg)
            self.bottom_gap_vec = gap_vec

        self._is_split = True

        return 0

    def set_meeting_point(self, point):

        self.meeting_point = point

        return 0

    def plot_segments(self, show=True, ax=None):

        if show:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if not self._is_split:
            for seg in self._top_segments:
                ax.plot([seg.ra[0], seg.rb[0]], [seg.ra[1], seg.rb[1]], color=colors[seg.color])
                ax.scatter(seg.ra[0], seg.ra[1], marker='o', color=colors[0])
                ax.scatter(seg.rb[0], seg.rb[1], marker='o', color=colors[0])

            for seg in self._bottom_segments:
                ax.plot([seg.ra[0], seg.rb[0]], [seg.ra[1], seg.rb[1]], color=colors[seg.color])
                ax.scatter(seg.ra[0], seg.ra[1], marker='o', color=colors[0])
                ax.scatter(seg.rb[0], seg.rb[1], marker='o', color=colors[0])
        else:
            all_segs = self.top_dee_segs + self.top_dummy_dee_segs + self.bottom_dee_segs + self.bottom_dummy_dee_segs
            for seg in all_segs:
                ax.plot([seg.ra[0], seg.rb[0]], [seg.ra[1], seg.rb[1]], color=colors[seg.color])
                ax.scatter(seg.ra[0], seg.ra[1], marker='o', color=colors[seg.color])
                ax.scatter(seg.rb[0], seg.rb[1], marker='o', color=colors[seg.color])

        if show:
            ax.grid(True)
            ax.set_xlim([-0.25, 0.25])
            ax.set_ylim([-0.25, 0.25])
            ax.set_aspect(1)
            ax.set_xlabel('x (m)')
            ax.set_ylabel('y (m)')
            plt.show()

        return 0


class DummyDee(PyElectrode):
    def __init__(self, parent1, parent2, gap=0.05, thickness=0.025, **kwargs):
        super().__init__(name="Dummy Dee", **kwargs)

        self.parent1 = parent1
        self.parent2 = parent2
        self.gap = gap
        self.thickness = thickness

        self.color = 'BLUE'

    def generate_geometry(self):
        geo_str, _ = generate_dee_geometry(self.parent1.top_dummy_dee_segs,
                                           self.parent2.bottom_dummy_dee_segs,
                                           gap=self.gap, thickness=self.thickness)

        self.generate_from_geo_str(geo_str)


class Dee(PyElectrode):
    def __init__(self, parent, gap=0.05, thickness=0.025, **kwargs):
        super().__init__(name="Dee", **kwargs)

        self.parent = parent
        self.gap = gap
        self.thickness = thickness

        self.color = 'GREEN'

    def generate_geometry(self):
        geo_str, _ = generate_dee_geometry(self.parent.top_dee_segs,
                                           self.parent.bottom_dee_segs,
                                           gap=self.gap, thickness=self.thickness)

        self.generate_from_geo_str(geo_str)


class CRSegment(object):
    def __init__(self, ra, rb, color=0):
        self.ra = ra  # Inner coordinate (x,y)
        self.rb = rb  # Outer coordinate (x,y)
        self.mid = 0.5 * (ra + rb)

        self.color = color

    def rotate(self, angle, angle_unit="deg"):

        if angle_unit == "deg":
            _angle = np.deg2rad(angle)
        elif angle_unit == "rad":
            _angle = angle
        else:
            return 1

        rot = np.array([[np.cos(_angle), -np.sin(_angle), 0.0],
                        [np.sin(_angle), np.cos(_angle), 0.0],
                        [0.0, 0.0, 1.0]])

        self.ra = np.matmul(rot, self.ra[np.newaxis].T)[:, 0]
        self.rb = np.matmul(rot, self.rb[np.newaxis].T)[:, 0]

        self.mid = 0.5 * (self.ra + self.rb)

        return 0

    def translate(self, dr):

        self.ra = self.ra + dr
        self.rb = self.rb + dr

        return 0

    def check_intersection(self, rc, rd):
        # 2D
        return check_intersection(self.ra[:2], self.rb[:2], rc[:2], rd[:2])


def check_intersection(q, b, p, d):
    s = b - q
    r = d - p

    rscross = np.cross(r, s)

    if rscross == 0:
        return False

    qmp = q - p

    t = np.cross(qmp, s) / rscross
    u = np.cross(qmp, r) / rscross

    if np.cross(r, s) != 0 and t >= 0 and t <= 1 and u >= 0 and u <= 1:
        return True
    else:
        return False


def generate_dee_geometry(top_segments, bottom_segments, h=0.025,
                          gap=0.05, thickness=0.025,
                          starting_indices=(0, 0, 0, 0, 1)):
    # indices[0] = Points
    # indices[1] = Lines
    # indices[2] = Line Loops
    # indices[3] = Surfaces
    # indices[4] = Volumes

    indices = list(starting_indices)

    geo_str = 'SetFactory("OpenCASCADE");\n'
    geo_str += "Mesh.CharacteristicLengthMax = {};\n".format(h)

    geo_str += gmsh_macros

    geo_str += "gap = {};\n".format(gap)
    geo_str += "thickness = {};\n".format(thickness)

    geo_str += "start_idx = {};\n".format(indices[0])

    # Points are 2D (x,y)
    # For N segments, there are N+1 points
    for i, seg in enumerate(top_segments):
        if i == 0:
            geo_str += "Point ({}) = {{ {}, {}, {} }};\n".format(indices[0], seg.ra[0], seg.ra[1], seg.ra[2])
            indices[0] += 1
        geo_str += "Point ({}) = {{ {}, {}, {} }};\n".format(indices[0], seg.rb[0], seg.rb[1], seg.rb[2])
        indices[0] += 1

    # Outer arc center point
    # r1, r2 = self.top_dee_segs[-1].rb, self.bottom_dee_segs[-1].rb
    # rc = 0.5*(r2 - r1)
    # ro = r1 + rc
    # R = 1.0  # TODO
    # a = np.linalg.norm(rc)
    # dr = np.array([rc[1], -rc[0]]) * np.sqrt((R/a)**2 - 1)
    # rp = ro + dr
    # geo_str += "Point ({}) = {{ {}, {}, 0.0 }};\n".format(indices[0], rp[0], rp[1])
    # outer_arc_center_idx = indices[0]
    # indices[0] += 1
    #
    # dtheta = np.arccos(np.dot(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2)))

    for i, seg in enumerate(bottom_segments[::-1]):
        geo_str += "Point ({}) = {{ {}, {}, {} }};\n".format(indices[0], seg.rb[0], seg.rb[1], seg.rb[2])
        indices[0] += 1

        if i == len(bottom_segments) - 1:
            geo_str += "Point ({}) = {{ {}, {}, {} }};\n".format(indices[0], seg.ra[0], seg.ra[1], seg.ra[2])
            indices[0] += 1

    geo_str += "\nend_idx = {};\n".format(indices[0] - 1)

    geo_str += "Call MakeDeeSurface;\n"
    geo_str += "Call MakeDeeVolume;\n"

    geo_str += "p1 = newp; Point(p1) = {{ {}, {}, 0.0 }};\n".format(top_segments[0].ra[0],
                                                                    top_segments[0].ra[1])
    geo_str += "p2 = newp; Point(p2) = {{ {}, {}, 0.0 }};\n".format(top_segments[0].mid[0],
                                                                    top_segments[0].mid[1])
    geo_str += "p3 = newp; Point(p3) = {{ {}, {}, 0.0 }};\n".format(bottom_segments[0].mid[0],
                                                                    bottom_segments[0].mid[1])
    geo_str += "p4 = newp; Point(p4) = {{ {}, {}, 0.0 }};\n".format(bottom_segments[0].ra[0],
                                                                    bottom_segments[0].ra[1])
    geo_str += "Call MakePost;\n"

    with open('_dee.geo', 'w') as f:
        f.write(geo_str)

    # self.generate_from_geo_str(geo_str=geo_str)
    # self.show()

    return geo_str, indices
