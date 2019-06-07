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


class Dee(object):
    def __init__(self, dummy=False, voltage=7.0E4):
        self.dummy = dummy
        self.voltage = voltage


class CentralRegion(object):
    def __init__(self, spiral_inflector=None):

        self._initialized = False

        self._spiral_inflector = None

        self._params_analytic = {}
        self._variables_analytic = {}
        self._variables_numerical = {}
        self._variables_track = {}
        self._segments = []

        self._tracked_trjs = []
        self._track_trjs_vel = []

        self.rf_phase = 0.0
        self.rf_freq = 32.8E6
        self.v_peak = 70.0E3

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

    def dee_crossing(self, rc, rd, v, t):
        ion = self._params_analytic["ion"]
        intersected = False
        initial_v = ion.v_m_per_s()
        for segment in self._segments:
            intersected = segment.check_intersection(rc[:2], rd[:2]) or intersected
        if not intersected:
            return v
        ttf = 1.0
        dE = np.abs(ion.q() * np.sin( 2 * np.pi * self.rf_freq * t + self.rf_phase) * self.v_peak * ttf)
        print("Gained energy {:.3f} keV".format(dE / 1E3))
        ion.calculate_from_energy_mev(ion._energy_mev + dE / 1E6)
        return ion.v_m_per_s() * v / initial_v

    def initialize(self):
        trj_design = self._spiral_inflector.analytic_variables["trj_design"]
        v_design = self._spiral_inflector.analytic_variables["v_design"]
        self._xi, self._vi = trj_design[-1, :], v_design[-1, :]

    def create_initial_dees(self):

        thetas = np.linspace(0, 7*np.pi / 4, 8) + np.pi / 8.0
        r_init, r_final = 0.05, 1

        for theta in thetas:
            ra = np.array([r_init * np.cos(theta), r_init * np.sin(theta)])
            rb = np.array([r_final * np.cos(theta), r_final * np.sin(theta)])
            dee_seg = CRSegment(ra=ra, rb=rb)
            self._segments.append(dee_seg)

    def set_inflector(self, spiral_inflector):
        self._spiral_inflector = spiral_inflector
        self._params_analytic["ion"] = self._spiral_inflector.analytic_parameters["ion"]
        self._params_analytic["bf_itp"] = self._spiral_inflector.analytic_parameters["bf_itp"]

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

    def plot_segments(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for seg in self._segments:
            ax.plot([seg.ra[0], seg.rb[0]], [seg.ra[1], seg.rb[1]], color=colors[0])
            ax.scatter(seg.ra[0], seg.ra[1], marker='o', color=colors[0])
            ax.scatter(seg.rb[0], seg.rb[1], marker='o', color=colors[0])

        ax.grid(True)
        ax.set_xlim([-0.25, 0.25])
        ax.set_ylim([-0.25, 0.25])
        ax.set_aspect(1)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')

        si = self._spiral_inflector

        trj_design = si.analytic_variables["trj_design"]
        ax.plot(trj_design[:, 0], trj_design[:, 1], color=colors[1])

        xc, yc = si.calculate_orbit_center()
        ax.scatter(xc, yc, marker='X', color=colors[2])

        for tracked_trj in self._tracked_trjs:
            ax.plot(tracked_trj[:, 0], tracked_trj[:, 1], color=colors[3])

        plt.show()

    def plot_bfield(self, nx=100, ny=100):
        x, y = np.linspace(-0.1, 0.1, nx), np.linspace(-0.1, 0.1, ny)
        B = np.zeros([nx, ny])
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                B[i, j] = self._params_analytic["bf_itp"](np.array([xi, yi, 0.0]))

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x, y, B, color='b')
        plt.show()

    def track(self, **kwargs):
        r, v = central_region_track(self, r_start=self._xi, v_start=self._vi, nsteps=5000, dt=1e-11, **kwargs)
        self._tracked_trjs.append(r)
        self._track_trjs_vel.append(v)
        return r, v


class CRSegment(object):
    def __init__(self, ra, rb):
        self.ra = ra  # Inner coordinate (x,y)
        self.rb = rb  # Outer coordinate (x,y)

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
