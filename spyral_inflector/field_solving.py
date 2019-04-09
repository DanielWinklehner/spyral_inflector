from dans_pymodules import *
import bempp.api
# noinspection PyUnresolvedReferences
from bempp.api.shapes.shapes import __generate_grid_from_geo_string as generate_from_string

# Define the directions:
X = 0
Y = 1
Z = 2


def calculate_efield(si,
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

    if si._variables_bempp["solution"] is None:
        print("Please solve with BEM++ before calculating the E-Field")
        return 1

    _ts = time.time()

    sol = si._variables_bempp["solution"]
    fsp = si._variables_bempp["f_space"]

    # noinspection PyUnresolvedReferences
    all_vert = si._variables_bempp["full mesh"].leaf_view.vertices

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
    if si._debug:
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

                if si._debug:
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

    si._variables_track["ef_itp"] = Field("Spiral Inflector E-Field",
                                          dim=3,
                                          field={"x": RegularGridInterpolator(points=_r, values=-ex,
                                                                              bounds_error=False, fill_value=0.0),
                                                 "y": RegularGridInterpolator(points=_r, values=-ey,
                                                                              bounds_error=False, fill_value=0.0),
                                                 "z": RegularGridInterpolator(points=_r, values=-ez,
                                                                              bounds_error=False, fill_value=0.0)
                                                 })

    return 0


def solve_bempp(si):
    if si._variables_bempp["full mesh"] is None:
        print("Please generate a mesh before solving with BEM++!")
        return 1

    print("Generating necessary BEM++ operators, function spaces. Solving... ", end="")

    dp0_space = bempp.api.function_space(si._variables_bempp["full mesh"], "DP", 0)
    slp = bempp.api.operators.boundary.laplace.single_layer(dp0_space, dp0_space, dp0_space)

    domain_mapping = {}
    for name, electrode in si._variables_bempp["objects"].items():
        domain_mapping[electrode["domain"]] = electrode["voltage"]

    def f(*args):

        domain_index = args[2]
        result = args[3]

        result[0] = domain_mapping[domain_index]

    dirichlet_fun = bempp.api.GridFunction(dp0_space, fun=f)

    if si._debug:
        dirichlet_fun.plot()

    # Solve
    sol, info = bempp.api.linalg.gmres(slp, dirichlet_fun, tol=1e-5, use_strong_form=True)

    print("Done!")

    # Save results
    si._variables_bempp["solution"] = sol
    si._variables_bempp["f_space"] = dp0_space
    si._variables_bempp["operator"] = slp
    si._variables_bempp["grid_fun"] = dirichlet_fun

    return 0
