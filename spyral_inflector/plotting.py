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

    analytic_vars = si.analytic_variables
    numerical_vars = si.numerical_variables
    track_vars = si.track_variables
    shift = track_vars["shift"]

    # --- Plot with pythonocc-core Qt5 Window --- #
    # Electrodes (PyElectrodeAssembly.show() returns an instance of the display)
    # for _, elec in numerical_vars["objects"].electrodes.items():
    #     print(elec.name, elec.occ_obj)

    numerical_vars["objects"].set_translation(shift, absolute=True)  # TODO: This should be applied as part of calculation
    display, start_display = numerical_vars["objects"].show()

    # Trajectories
    occ_trj = SITrajectory(name="Tracked Design Trajectory", voltage=0)
    occ_trj.set_translation(shift, absolute=True)
    occ_trj.create_geo_str(analytic_vars["trj_design"], max_points=freq, load=True)
    occ_trj.color = "RED"
    occ_trj.show(display=display)

    if track_vars["trj_tracker"] is not None:
        occ_trj = SITrajectory(name="Tracked Design Trajectory", voltage=0)
        occ_trj.set_translation(shift, absolute=True)
        occ_trj.create_geo_str(track_vars["trj_tracker"], max_points=freq, load=True)
        occ_trj.color = "BLACK"
        occ_trj.show(display=display)

    if aux_trajectories is not None:
        for i, _trj in enumerate(aux_trajectories):
            occ_trj = SITrajectory(name="Aux Trajectory {}".format(i), voltage=0)
            occ_trj.set_translation(shift, absolute=True)
            occ_trj.create_geo_str(_trj, max_points=freq, load=True)
            occ_trj.color = "BLUE"
            occ_trj.show(display=display)

    # Repaint
    display.FitAll()
    display.Repaint()

    # Show interactive display. Currently, this locks the execution of the script.
    # without start_display, it is still created and pops up briefly, then is destroyed immediately...
    if show:
        start_display()

    if filename is not None:
        if filename == 'auto':
            _fd = FileDialog()
            filename = _fd.get_filename(action="save")
        try:
            display.ExportToImage(filename)
        except Exception as ex:
            print("Something went wrong when trying to save the file: {}".format(ex))
            return 1

    return 0


# TODO: Revamp potential plotting
# def plot_potential(si, limits=(-0.1, 0.1, -0.1, 0.1), orientation="xy", **kwargs):
#     if si._variables_numerical["solution"] is None:
#         print("No BEM++ solution in memory to plot from.")
#         return 1
#
#     if "n_grid_points" in kwargs.keys():
#         n_grid_points = kwargs["n_grid_points"]
#     else:
#         n_grid_points = 200
#
#     if "offset" in kwargs.keys():
#         offset = kwargs["offset"]
#     else:
#         offset = 0.0
#
#     sol = si._variables_numerical["solution"]
#     f_space = si._variables_numerical["f_space"]
#
#     plot_grid = np.mgrid[limits[0]:limits[1]:n_grid_points * 1j,
#                 limits[2]:limits[3]:n_grid_points * 1j]
#
#     e1 = plot_grid[0].ravel()
#     e2 = plot_grid[1].ravel()
#     e3 = offset * np.ones(plot_grid[0].size)
#     if orientation == "xy":
#         points = np.vstack((e1, e2, e3))
#     elif orientation == "yz":
#         points = np.vstack((e3, e1, e2))
#     elif orientation == "xz":
#         points = np.vstack((e1, e3, e2))
#     else:
#         print("Orientation not recognized. Allowed orientations are: xy, yz, xz.")
#         return 1
#
#     slp_pot = bempp.api.operators.potential.laplace.single_layer(f_space, points)
#     u_evaluated = slp_pot * sol
#     u_evaluated = u_evaluated.reshape((n_grid_points, n_grid_points))
#
#     if "vlims" in kwargs.keys():
#         vlims = kwargs["vlims"]
#     else:
#         vlims = [np.min(u_evaluated), np.max(u_evaluated)]
#
#     plt.imshow(u_evaluated.T,
#                extent=(limits[0], limits[1], limits[2], limits[3]),
#                vmin=vlims[0],
#                vmax=vlims[1])
#
#     if "colorbar" in kwargs.keys():
#         if kwargs["colorbar"]:
#             plt.colorbar()
#
#     plt.show()
#
#     if "save_fig" in kwargs.keys() and "filename" in kwargs.keys():
#         if kwargs["save_fig"]:
#             plt.savefig(kwargs["filename"])
#
#     return 0


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
    if si._variables_numerical["full mesh"] is not None:
        si._variables_numerical["full mesh"].plot()
