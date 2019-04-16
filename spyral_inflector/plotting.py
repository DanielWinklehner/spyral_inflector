import bempp.api
from dans_pymodules import *
from mpl_toolkits.mplot3d import proj3d
from .geometry import SITrajectory

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


def draw_geometry(si, freq=10, show=False, filename=None, aux_trajectories=None):

    if si.analytic_variables["geo"] is None or si.analytic_variables["trj_design"] is None:
        print("No geometry yet ... generating")
        si.generate_geometry()

    analytic_pars = si.analytic_parameters
    analytic_vars = si.analytic_variables
    bempp_pars = si.bempp_parameters
    bempp_vars = si.bempp_variables
    track_pars = si.track_parameters
    track_vars = si.track_variables

    trj_design = analytic_vars["trj_design"]  # type: np.ndarray
    geo = analytic_vars["geo"]  # type: np.ndarray
    shift = track_vars["shift"]  # type: np.ndarray

    if shift is None:
        shift = np.zeros(3)

    # --- Plot with pythonocc-core Qt5 Window --- #
    display, start_display = bempp_vars["objects"].show()
    occ_trj = SITrajectory(name="Design Trajectory", voltage=0)
    occ_trj.create_geo_str(trj_design, max_points=25, load=True)
    occ_trj.color = "BLACK"
    occ_trj.show(display=display)

    display.FitAll()
    display.Repaint()

    start_display()

    input()
    exit()

    # Plot with matplotlib
    fig = plt.figure()
    ax = Axes3D(fig)

    proj3d.persp_transformation = orthogonal_proj

    ns = si._params_analytic["ns"]  # type: int

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

    if si._variables_track["trj_tracker"] is not None:
        trj_tracked = si._variables_track["trj_tracker"]

        # Plot the tracked beam trajectory
        ax.plot(trj_tracked[:, 0] + shift[0],
                trj_tracked[:, 1] + shift[1],
                trj_tracked[:, 2] + shift[2],
                color=colors[1], linewidth=3.0)

    if aux_trajectories is not None:
        for trj in aux_trajectories:
            ax.plot(trj[:, 0], trj[:, 1], trj[:, 2], linewidth=1.0)

    plt.axis("equal")
    plt.xlabel("x")
    plt.ylabel("y")
    ax.set_xlim(-0.1, 0.1)
    ax.set_ylim(-0.1, 0.1)
    ax.set_zlabel("z")
    ax.set_aspect("equal")

    ax.view_init(elev=-180.0, azim=-90.0)

    if filename == 'dialog':
        from dans_pymodules import FileDialog
        print("Please select figure filename in dialog...")
        fd = FileDialog()
        filename = fd.get_filename('open')
    elif filename == 'auto':
        filename = os.path.join(si._outp_folder, "geometry.png")

    if filename is not None:
        try:
            fig.savefig(fname=filename)
        except Exception as ex:
            print("Something went wrong when trying to save the file: {}".format(ex))
            return 1

    if show:
        plt.show()

    return 0


def plot_potential(si, limits=(-0.1, 0.1, -0.1, 0.1), orientation="xy", **kwargs):
    if si._variables_bempp["solution"] is None:
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

    sol = si._variables_bempp["solution"]
    f_space = si._variables_bempp["f_space"]

    plot_grid = np.mgrid[limits[0]:limits[1]:n_grid_points * 1j,
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
    else:
        print("Orientation not recognized. Allowed orientations are: xy, yz, xz.")
        return 1

    slp_pot = bempp.api.operators.potential.laplace.single_layer(f_space, points)
    u_evaluated = slp_pot * sol
    u_evaluated = u_evaluated.reshape((n_grid_points, n_grid_points))

    if "vlims" in kwargs.keys():
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
            plt.savefig(kwargs["filename"])

    return 0


def plot_bfield(si, lims=(-0.2, 0.2), num=5000):
    zpts = np.linspace(lims[0], lims[1], num)
    b = np.zeros(num)
    plt.rc('text', usetex=True)
    plt.rc('font', family="serif")
    matplotlib.rcParams.update({'font.size': 16})
    fig, ax = plt.subplots()
    for i in range(num):
        zpt = zpts[i]
        b[i] = si._params_analytic["bf_itp"](np.array([0.0, 0.0, zpt]))[2]
    ax.plot(np.array(zpts), -np.array(b), colors[0], linewidth=3.0)
    plt.xlim([-0.2, 0.2])
    plt.xlabel('z (m)')
    plt.ylabel('B (T)')
    plt.grid(True)
    plt.savefig("bfield.png")


def plot_trajectories(si, trajectories=None):
    # Quick method for plotting trajectories
    fig = plt.figure()
    ax = Axes3D(fig)

    proj3d.persp_transformation = orthogonal_proj

    # Plot the beam trajectory
    for trj in trajectories:
        ax.plot(trj[:, 0], trj[:, 1], trj[:, 2])

    plt.show()


def draw_mesh(si):
    if si._variables_bempp["full mesh"] is not None:
        si._variables_bempp["full mesh"].plot()
