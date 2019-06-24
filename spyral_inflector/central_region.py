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


class CentralRegion(object):
    def __init__(self, spiral_inflector=None, r_cr=0.3):

        self._initialized = False

        self._spiral_inflector = None

        self._params_analytic = {}
        self._variables_analytic = {}
        self._variables_numerical = {}
        self._variables_track = {}

        self._dees = []

        self._tracked_trjs = []
        self._track_trjs_vel = []

        self.rf_phase = 0.0
        self.rf_freq = 32.8E6
        self.harmonic = 4
        self.v_peak = 70.0E3
        self.r_cr = r_cr  # Radius at which we no longer consider the central region as a different entity

        if spiral_inflector is not None:
            self.set_inflector(spiral_inflector)

        self._xi, self._vi = None, None

    @property
    def analytic_parameters(self):
        return self._params_analytic

    @property
    def analytic_variables(self):
        return self._variables_analytic

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

    def initialize(self, xi=None, vi=None):
        if xi is None and vi is None:
            trj_design = self._spiral_inflector.analytic_variables["trj_design"]
            v_design = self._spiral_inflector.analytic_variables["v_design"]
            self._xi, self._vi = trj_design[-1, :], v_design[-1, :]
        else:
            self._xi, self._vi = xi, vi

    def set_inflector(self, spiral_inflector):
        self._spiral_inflector = spiral_inflector
        self._params_analytic["ion"] = self._spiral_inflector.analytic_parameters["ion"]
        self._params_analytic["bf_itp"] = self._spiral_inflector.analytic_parameters["bf_itp"]

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
        self._dees.append(dee)

        return 1

    def plot_dees(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for dee in self._dees:
            dee.plot_segments(show=False, ax=ax)

        ax.grid(True)
        ax.set_xlim([-self.r_cr, self.r_cr])
        ax.set_ylim([-self.r_cr, self.r_cr])
        ax.set_aspect(1)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        plt.show()

    def plot_bfield(self, nx=100, ny=100):
        from matplotlib import cm

        x, y = np.linspace(-0.1, 0.1, nx), np.linspace(-0.1, 0.1, ny)
        B = np.zeros([nx, ny])
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                B[i, j] = self._params_analytic["bf_itp"](np.array([xi, yi, 0.0]))[2]

        fig = plt.figure()
        xx, yy = np.meshgrid(x, y)
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(xx, yy, B, cmap=cm.coolwarm)
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


class Dee(object):

    def __init__(self):

        self.meeting_point = [0.0, 0.0]
        self.opening_angle = 45.0
        self._opening_angle = np.deg2rad(self.opening_angle)
        self.angle = 0.0

        self._top_segments, self._bottom_segments = [], []
        self._char_len = 0.05  # Characteristic length of CRSegments

    def rotate(self, angle, angle_unit="deg"):
        """
        Rotates the entire Dee object by angle
        :param angle:
        :return:
        """
        self.angle = angle

        for seg in self._top_segments + self._bottom_segments:
            seg.rotate(angle, angle_unit=angle_unit)

        return 0

    def initialize(self):

        # Initial top segment
        ra = self._char_len * np.array([np.cos(self._opening_angle / 2.0),
                                        np.sin(self._opening_angle / 2.0)])

        rb = 2 * self._char_len * np.array([np.cos(self._opening_angle / 2.0),
                                            np.sin(self._opening_angle / 2.0)])

        top_seg = CRSegment(ra, rb)

        # Initial bottom segment
        ra = self._char_len * np.array([np.cos(-self._opening_angle / 2.0),
                                        np.sin(-self._opening_angle / 2.0)])

        rb = 2 * self._char_len * np.array([np.cos(-self._opening_angle / 2.0),
                                            np.sin(-self._opening_angle / 2.0)])

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

        ra = np.array([0.0, 0.0])
        rb = length * np.array([np.cos(_angle), np.sin(_angle)])

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

        ra = np.array([0.0, 0.0])
        rb = length * np.array([np.cos(_angle), np.sin(_angle)])

        prev_seg = self._top_segments[-1]
        next_seg = CRSegment(ra, rb, color=1)

        next_seg.rotate(self.angle + self.opening_angle / 2.0, angle_unit="deg")
        next_seg.translate(prev_seg.rb)

        self._top_segments.append(next_seg)

        return 0

    def split(self):
        pass

    # def add_seg(self, seg, dummy=False):
    #     if dummy:
    #         self.dummy_dee_segs.append(seg)
    #     else:
    #         self.dee_segs.append(seg)
    #
    #     return 0

    def set_meeting_point(self, point):
        self.meeting_point = point

        return 0

    def plot_segments(self, show=True, ax=None):

        if show:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for seg in self._top_segments:
            ax.plot([seg.ra[0], seg.rb[0]], [seg.ra[1], seg.rb[1]], color=colors[seg.color])
            ax.scatter(seg.ra[0], seg.ra[1], marker='o', color=colors[0])
            ax.scatter(seg.rb[0], seg.rb[1], marker='o', color=colors[0])

        for seg in self._bottom_segments:
            ax.plot([seg.ra[0], seg.rb[0]], [seg.ra[1], seg.rb[1]], color=colors[seg.color])
            ax.scatter(seg.ra[0], seg.ra[1], marker='o', color=colors[0])
            ax.scatter(seg.rb[0], seg.rb[1], marker='o', color=colors[0])

        if show:
            ax.grid(True)
            ax.set_xlim([-0.25, 0.25])
            ax.set_ylim([-0.25, 0.25])
            ax.set_aspect(1)
            ax.set_xlabel('x (m)')
            ax.set_ylabel('y (m)')
            plt.show()

        return 0


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

        _ra, _rb = np.array([self.ra[0], self.ra[1], 1.0]), np.array([self.rb[0], self.rb[1], 1.0])

        new_ra = np.matmul(rot, _ra[np.newaxis].T)[:2, 0]
        new_rb = np.matmul(rot, _rb[np.newaxis].T)[:2, 0]

        self.ra = new_ra
        self.rb = new_rb
        self.mid = 0.5 * (new_ra + new_rb)

        return 0

    def translate(self, dr):

        self.ra = self.ra + dr
        self.rb = self.rb + dr

        return 0

    def check_intersection(self, rc, rd):

        return check_intersection(self.ra, self.rb, rc, rd)


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
