from .field_solving import *
from .geometry import *
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

Macro MakeSomething

    // Macro for interfacing the end of the spiral inflector to the beginning of the dee system
    // Suppose there are a start index and end index for the points in the trajectory
    

Return

"""

two_dim_dee_gmsh_str = """

// Top dee electrode
//Point(1) = {-(h_gap / 2 + dee_len), v_gap / 2 + dee_thk, 0.0};
Point(2) = {-(h_gap / 2 + dee_len), v_gap / 2, 0.0};
Point(3) = {-(h_gap / 2), v_gap / 2, 0.0};
Point(4) = {-(h_gap / 2), v_gap / 2 + dee_thk, 0.0};

//Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
//Line(4) = {4, 1};
//Line Loop(1) = {1:4};

// Bottom dee electrode
//Point(5) = {-(h_gap / 2 + dee_len), -(v_gap / 2 + dee_thk), 0.0};
Point(6) = {-(h_gap / 2 + dee_len), -(v_gap / 2), 0.0};
Point(7) = {-(h_gap / 2), -(v_gap / 2), 0.0};
Point(8) = {-(h_gap / 2), -(v_gap / 2 + dee_thk), 0.0};

//Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
//Line(8) = {8, 5};
//Line Loop(2) = {5:8};

// Top dummy dee electrode
//Point(9) = {h_gap / 2 + dee_len, v_gap / 2 + dee_thk, 0.0};
Point(10) = {h_gap / 2 + dee_len, v_gap / 2, 0.0};
Point(11) = {h_gap / 2, v_gap / 2, 0.0};
Point(12) = {h_gap / 2, v_gap / 2 + dee_thk, 0.0};

//Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
//Line(12) = {12, 9};
//Line Loop(3) = {9:12};

// Bottom dummy dee electrode
//Point(13) = {h_gap / 2 + dee_len, -(v_gap / 2 + dee_thk), 0.0};
Point(14) = {h_gap / 2 + dee_len, -(v_gap / 2), 0.0};
Point(15) = {h_gap / 2, -(v_gap / 2), 0.0};
Point(16) = {h_gap / 2, -(v_gap / 2 + dee_thk), 0.0};

//Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
//Line(16) = {16, 13};
//Line Loop(4) = {13:16};

// Neumann boundary lines, clockwise

Line(100) = {6, 2}; // left
Line(101) = {16, 8}; // bottom
Line(102) = {4, 12}; // top
Line(103) = {10, 14}; // right

Line Loop(1) = {2, 3, 102, -11, -10, 103, 14, 15, 101, -7, -6, 100};
Surface(1) = {1};

Physical Curve(1) = {2, 3, 6, 7, 100}; // Dee
Physical Curve(2) = {10, 11, 14, 15, 103}; // Dummy Dee
Physical Curve(3) = {102, 101}; // Top/Bottom boundary
//Physical Curve(4) = {100, 103}; // Left/Right boundary

Physical Surface(1) = {1}; // Domain

"""


class CentralRegion(PyElectrodeAssembly):
    def __init__(self, spiral_inflector=None,
                 r_cr=(0.1, 0.3),
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
        self.dee_gap = 0.056
        self.r_cr = r_cr
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

    def initialize(self, xi=None, vi=None, dee_z_func=None):
        if xi is None and vi is None:
            trj_design = self._spiral_inflector.analytic_variables["trj_design"]
            v_design = self._spiral_inflector.analytic_variables["trj_vel"]
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

    def plot_bfield(self, nx=100, ny=100, xlims=(-0.3, 0.3), ylims=(-0.3, 0.3), z=0.0, fig=None, ax=None):
        from matplotlib import cm

        x, y = np.linspace(xlims[0], xlims[1], nx), \
               np.linspace(ylims[0], ylims[1], ny)
        B = np.zeros([ny, nx])
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                B[j, i] = self._params_analytic["bf_itp"](np.array([xi, yi, z]))[2]

        xx, yy = np.meshgrid(x, y)
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sc = ax.contour(xx, yy, B, levels=15, cmap=cm.coolwarm, vmin=np.min(B), vmax=np.max(B))
            ax.set_aspect(1)
            fig.colorbar(sc)
            plt.show()
        else:
            sc = ax.contour(xx, yy, B, levels=15, cmap=cm.coolwarm, vmin=np.min(B), vmax=np.max(B))
            ax.set_aspect(1)
            cbar = fig.colorbar(sc)
            cbar.set_label("B (T)")

    def track(self, **kwargs):
        r, v = modulated_track(self, r_start=self._xi, v_start=self._vi, dt=1e-11, **kwargs)
        self._tracked_trjs.append(r)
        self._track_trjs_vel.append(v)

        return r, v

    def track_segment(self, **kwargs):

        r, v = track_segment(self, r_start=self._xi, v_start=self._vi, dt=1e-11, **kwargs)
        self._tracked_trjs.append(r)
        self._track_trjs_vel.append(v)

        return r, v

    def make_dees(self, dees, n, gap, thickness, cl=0.1, **kwargs):

        # cl = 0.1

        self._abstract_dees = dees

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

    def track_from_si(self, nsteps=1000, dt=1e-11):
        si = self._spiral_inflector

        r_start, v_start = si.analytic_variables["trj_design"][-1, :], \
                           si.analytic_variables["trj_vel"][-1, :]

        r, v = cr_track(cr=self, r_init=r_start, v_init=v_start,
                        end_type="steps", maxsteps=nsteps, dt=dt,
                        input_bfield=self.analytic_parameters["bf_itp"])

        return r, v

    def simple_track(self, r_start, v_start):
        return central_region_simple_track(self, r_start, v_start)

    def optimize(self):
        pass


class TwoDeeField(object):
    def __init__(self, gap_center_angle=0.0, left_voltage=70e3, right_voltage=0.0, h_gap=0.01, v_gap=0.05):

        assert HAVE_FENICS, "This object needs fenics for the field calculations!"
        # TODO: Assert meshio

        # Calculate a two dimensional field with the dee parameters.
        # Create a copy and multiply the field by -1 to represent DD-->D
        # Dee field is always a perpendicular slice w.r.t dee edge
        # --> Need a geometrical factor for off-angle trajectories

        self._id = uuid.uuid1()

        self._gap_center_angle = gap_center_angle
        self._left_voltage = left_voltage
        self._right_voltage = right_voltage
        self._h_gap = h_gap
        self._v_gap = v_gap

        gmsh_str_prelude = """
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMax = 0.0025;
dee_len = {};
h_gap = {};
v_gap = {};
dee_thk = {};

""".format(0.15, self._h_gap, self._v_gap, 0.02)

        self._gmsh_str = gmsh_str_prelude + two_dim_dee_gmsh_str

        with open(TEMP_DIR + '/{}.geo'.format(self._id), 'w') as f:
            f.write(self._gmsh_str)

        os.system('gmsh -2 ' + TEMP_DIR + '/{}.geo -format msh2 -v 0 -o '.format(self._id) + TEMP_DIR + '/{}.msh'.format(self._id))
        msh = meshio.read(TEMP_DIR + "/{}.msh".format(self._id))
        meshio.write_points_cells(TEMP_DIR + "/{}_markers.xdmf".format(self._id),
                                  msh.points,
                                  {"triangle": msh.cells["triangle"]},
                                  cell_data={"triangle": {"gmsh:physical": msh.cell_data["triangle"]["gmsh:physical"]}})

        meshio.write_points_cells(TEMP_DIR + "/{}_boundaries.xdmf".format(self._id),
                                  msh.points,
                                  {"line": msh.cells["line"]},
                                  cell_data={"line": {"gmsh:physical": msh.cell_data["line"]["gmsh:physical"]}})

        mesh = fn.Mesh()
        fn.XDMFFile(TEMP_DIR + "/{}_markers.xdmf".format(self._id)).read(mesh)

        markers = fn.MeshFunction("size_t", mesh, mesh.topology().dim())
        fn.XDMFFile(TEMP_DIR + "/{}_markers.xdmf".format(self._id)).read(markers, "gmsh:physical")

        _boundaries = fn.MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
        fn.XDMFFile(TEMP_DIR + "/{}_boundaries.xdmf".format(self._id)).read(_boundaries, "gmsh:physical")
        boundaries = fn.MeshFunction("size_t", mesh, _boundaries)

        dx = fn.Measure('dx', domain=mesh, subdomain_data=markers)
        V = fn.FunctionSpace(mesh, 'P', 1)

        bcs = [fn.DirichletBC(V, fn.Constant(self._left_voltage), boundaries, 1),
               fn.DirichletBC(V, fn.Constant(self._right_voltage), boundaries, 2)]

        # Test and trial functions
        u = fn.TrialFunction(V)
        v = fn.TestFunction(V)

        a = fn.dot(fn.grad(u), fn.grad(v)) * dx
        L = fn.Constant("0.0") * v * dx

        u = fn.Function(V)

        fn.solve(a == L, u, bcs, solver_parameters={"linear_solver": "cg", "preconditioner": "ilu"})
        # print("Done!", flush=True)

        # potentialFile = fn.File(TEMP_DIR + '/{}_potential.pvd'.format(self._id))
        # potentialFile << u
        #
        # meshfile = fn.File(TEMP_DIR + '/{}_mesh.pvd'.format(self._id))
        # meshfile << mesh

        fenics_field = fn.project(-fn.grad(u), solver_type='cg', preconditioner_type='ilu')
        electric_field = FenicsField(fenics_field)

        self._potential = u
        self._efield = electric_field

        # xr = np.linspace(-0.1, 0.1, 100)
        # ef = []
        # for x in xr:
        #     ef.append(electric_field(np.array([x, 0.0, 0.0]))[0])
        #
        # plt.plot(xr, ef)
        # plt.plot(xr, np.max(ef) * np.exp(-650 * xr ** 2))
        # plt.show()


class SimpleGap(object):  # TODO: Might not need this class for CR, higher N turns might be better
    def __init__(self, segments, phase=0.0):
        self.segments = segments
        self.phase = phase

    def check_intersection(self, start, end):
        for seg in self.segments:
            if check_intersection(seg.ra, seg.rb, start, end):
                return True

        return False


class Sectors(object):
    def __init__(self, abstract_dees):
        self.abstract_dees = abstract_dees
        self.lookup = None

        self.initialize()

    def initialize(self):

        first_dee = self.abstract_dees[0]
        last_dee = self.abstract_dees[-1]

        lookup = {}

        # [(0.0, 0.05): <dee1 ... deen>|segments in r_range, (0.05, 0.1): ...]

        char_len = 0.05

        n_segments = 4
        for k in range(n_segments):
            tup = (char_len*(k + 1), char_len*(k + 2))  # Skip over [0, char_len] by adding a factor of char_len
            ordered_dee_segments = []
            for dee in self.abstract_dees:

                top_seg = dee._top_segments[k]
                bottom_seg = dee._bottom_segments[k]

                if dee is first_dee:
                    ordered_dee_segments.append(top_seg)
                elif dee is last_dee:
                    ordered_dee_segments.append(bottom_seg)
                    ordered_dee_segments.append(top_seg)
                    ordered_dee_segments.append(first_dee._bottom_segments[k])
                else:
                    ordered_dee_segments.append(bottom_seg)
                    ordered_dee_segments.append(top_seg)

            lookup[tup] = ordered_dee_segments

        self.lookup = lookup

    def get_sector(self, pos, test=False):

        pos_r = np.sqrt(pos[0]**2 + pos[1]**2)  # Particle radial position

        rad_idx, next_gap, prev_gap = None, None, None

        for k in self.lookup.keys():  # Loop through each set of rmin, rmax
            rmin, rmax = k
            if rmin <= pos_r < rmax:  # If pos is in the r range
                for i, v in enumerate(self.lookup[k]):  # Enumerate each set of segments for this r range
                    ra, rb = v.ra, v.rb

                    # For checking if pos is 'before' or 'after' the line
                    m = (rb[1] - ra[1]) / (rb[0] - ra[0])
                    b = ra[1] - m * ra[0]

                    if rb[0] < ra[0]:  # Account for the direction of rb w.r.t. ra to get the direction right
                        dir = -1
                    else:
                        dir = 1

                    line_at_posx = m*pos[0] + b  # Find the y position of the line at the x position of the particle
                    ydiff = dir * (pos[1] - line_at_posx)  # negative = before
                    # print("-> ", ydiff)

                    if i != 0:  # Skip the first gap, since there is no 'previous' gap with these indices
                        if prev_ydiff >= 0 and ydiff < 0.0:
                            # print("^ I'm here!")
                            rad_idx = k
                            next_gap = i  # Save the gap indices
                            prev_gap = i - 1
                            print(rmin, rmax)
                            break  # Break out of the first for loop

                    prev_ydiff = ydiff  # Save the y difference from this iteration for the next one
                else:  # If there is no break, then the particle must be between the last and first gaps
                    # print("I'm in between the first and last!")
                    rad_idx = k
                    next_gap = 0  # In this case, we know what the indices must be
                    prev_gap = len(self.lookup[k]) - 1

            if rad_idx is not None:  # If we found where the particle is, then break
                break

        r_range = rad_idx  # TODO: Clean up these variable names... -PW
        segs_at_r = self.lookup[rad_idx]

        prev_seg = segs_at_r[prev_gap]
        next_seg = segs_at_r[next_gap]

        if test:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlim([-0.2, 0.2])
            ax.set_ylim([-0.2, 0.2])
            ax.set_aspect(1)
            ax.grid(True)

            t = np.deg2rad(np.linspace(0, 359, 360))
            c1 = r_range[0] * np.array([np.cos(t), np.sin(t)]).T
            c2 = r_range[1] * np.array([np.cos(t), np.sin(t)]).T

            ax.plot(c1[:, 0], c1[:, 1], '--', color='g')
            ax.plot(c2[:, 0], c2[:, 1], '--', color='g')

            for dee in self.abstract_dees:
                dee.plot_segments(ax=ax, show=False)

            ax.plot([prev_seg.ra[0], prev_seg.rb[0]], [prev_seg.ra[1], prev_seg.rb[1]], color='orange')
            ax.plot([next_seg.ra[0], next_seg.rb[0]], [next_seg.ra[1], next_seg.rb[1]], color='m')

            ax.scatter(pos[0], pos[1], marker='X', color='k')

            plt.show()


class AbstractDee(PyElectrode):
    def __init__(self,
                 r_init=0.1,
                 char_len=0.03,
                 angle=0.0,
                 opening_angle=30.0,
                 gap=0.05,
                 thickness=0.025,
                 **kwargs):
        super().__init__(name="Dee", **kwargs)

        self._r_init = r_init

        self.opening_angle = opening_angle
        self._opening_angle = np.deg2rad(self.opening_angle)
        self.angle = angle
        self._angle = np.deg2rad(angle)

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

        # For the 2D field tracking method
        # self.angle is the angle of the center of the dee
        # self.opening_angle is the angular span of the dee
        self._top_st = None
        self._bottom_st = None

    def make_gaps(self, phase=0.0):

        bottom_gap = SimpleGap(segments=self._bottom_segments, phase=phase)
        top_gap = SimpleGap(segments=self._top_segments, phase=phase + self._opening_angle)

        return bottom_gap, top_gap

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

        ra = self._r_init * np.array([np.cos(self._opening_angle / 2.0),
                                      np.sin(self._opening_angle / 2.0),
                                      0.0])
        rb = (self._r_init + self._char_len) * np.array([np.cos(self._opening_angle / 2.0),
                                                         np.sin(self._opening_angle / 2.0),
                                                         0.0])
        top_seg = CRSegment(ra, rb)

        # Initial bottom segment
        ra = self._r_init * np.array([np.cos(-self._opening_angle / 2.0),
                                      np.sin(-self._opening_angle / 2.0),
                                      0.0])
        rb = (self._r_init + self._char_len) * np.array([np.cos(-self._opening_angle / 2.0),
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

    def make_transforms(self):

        top_transforms = []
        for segment in self._top_segments:
            ra = segment.ra  # Initial point of the segment
            rb = segment.rb  # Final point of the segment
            dr = rb - ra  # Vector connecting the two points
            s = np.linalg.norm(dr)  # Length of the segment
            dth = -np.arccos(dr[0] / s) * np.sign(dr[1])  # Angle between x axis and segment vector
            dx, dy, _ = -ra  # Offset from origin

            # T = translation matrix
            # S = scaling matrix
            # R = rotation matrix
            T = np.array([[1.0, 0.0, dx], [0.0, 1.0, dy], [0.0, 0.0, 1.0]])
            S = np.array([[1.0 / s, 0.0, 0.0], [0.0, 1.0 / s, 0.0], [0.0, 0.0, 1.0]])
            R = np.array([[np.cos(dth), -np.sin(dth), 0.0], [np.sin(dth), np.cos(dth), 0.0], [0.0, 0.0, 1.0]])

            # M = overall transformation in order: T, R, S
            M = np.matmul(S, np.matmul(R, T))
            top_transforms.append(M)

        bottom_transforms = []
        for segment in self._bottom_segments:
            ra = segment.ra  # Initial point of the segment
            rb = segment.rb  # Final point of the segment
            dr = rb - ra  # Vector connecting the two points
            s = np.linalg.norm(dr)  # Length of the segment
            dth = -np.arccos(dr[0] / s) * np.sign(dr[1])  # Angle between x axis and segment vector
            dx, dy, _ = -ra  # Offset from origin

            # T = translation matrix
            # S = scaling matrix
            # R = rotation matrix
            T = np.array([[1.0, 0.0, dx], [0.0, 1.0, dy], [0.0, 0.0, 1.0]])
            S = np.array([[1.0 / s, 0.0, 0.0], [0.0, 1.0 / s, 0.0], [0.0, 0.0, 1.0]])
            R = np.array([[np.cos(dth), -np.sin(dth), 0.0], [np.sin(dth), np.cos(dth), 0.0], [0.0, 0.0, 1.0]])

            # M = overall transformation in order: T, R, S
            M = np.matmul(S, np.matmul(R, T))
            bottom_transforms.append(M)

        self._top_st = (self._top_segments, top_transforms)
        self._bottom_st = (self._bottom_segments, bottom_transforms)

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
