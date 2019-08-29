from dans_pymodules import *
from py_electrodes.py_electrodes import *

# Define the directions:
X, Y, Z = 0, 1, 2

XYZ = range(3)


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
    # fn.set_log_level(60)
    HAVE_FENICS = True
except ImportError as e:
    print(e)
    fn = None


HAVE_MESHIO = False
try:
    import meshio
    HAVE_MESHIO = True
except ImportError as e:
    print(e)
    meshio = None


class FenicsField(object):  # TODO: Implement this into the dans_pymodules Field class -PW
    def __init__(self, field):
        self._field = field
        self._field.set_allow_extrapolation(True)

    def __call__(self, position):
        return self._field(position[0], position[1], position[2])


def calculate_efield_bempp(si):  # TODO: Change 'si' to a more generic name (works for cr too) -PW

    numerical_vars = si.numerical_variables
    phi = numerical_vars["ef_phi"]

    assert phi is not None, "Please calculate the potential first!"

    print("Calculating the electric field... ", end="", flush=True)

    _d = numerical_vars["d"]
    _n = numerical_vars["n"]

    limits = numerical_vars["limits"]

    _r = np.array([np.linspace(limits[i, 0], limits[i, 1], _n[i]) for i in XYZ])

    ex, ey, ez = np.gradient(phi, _d[X], _d[Y], _d[Z])

    _field = Field("E-Field",
                   dim=3,
                   field={"x": RegularGridInterpolator(points=_r, values=-ex,
                                                       bounds_error=False, fill_value=0.0),
                          "y": RegularGridInterpolator(points=_r, values=-ey,
                                                       bounds_error=False, fill_value=0.0),
                          "z": RegularGridInterpolator(points=_r, values=-ez,
                                                       bounds_error=False, fill_value=0.0)
                          })

    numerical_vars["ef_itp"] = _field
    si.numerical_variables = numerical_vars

    print("Done!", flush=True)

    return 0


def calculate_efield_fenics(si):

    numerical_vars = si.numerical_variables
    u = numerical_vars["ef_phi"]

    assert u is not None, "Please calculate the potential first!"

    print("Calculating the electric field... ", end="", flush=True)

    fenics_field = fn.project(-fn.grad(u), solver_type='cg', preconditioner_type='ilu')
    electric_field = FenicsField(fenics_field)

    numerical_vars["ef_itp"] = electric_field
    si.numerical_variables = numerical_vars

    if si.debug:
        fieldfile = fn.File(TEMP_DIR + '/e_field.pvd')
        fieldfile << fenics_field

    print("Done!", flush=True)

    return 0


def calculate_potential(si,  # TODO: Change 'si' to a more generic name (works for cr too) -PW
                        limits=((None, None), (None, None), (None, None)),
                        res=0.002,
                        domain_decomp=(4, 4, 4),
                        overlap=0):
    # TODO: Add some of the debug stuff back in

    assert HAVE_BEMPP, "BEMPP not found. Aborting!"

    limits = np.array(limits)

    numerical_vars = si.numerical_variables

    if None in limits and numerical_vars["limits"] is not None:
        limits = numerical_vars["limits"]

    if limits.shape != (3, 2):
        print("Wrong shape of limits: {}. "
              "Must be ((xmin, xmax), (ymin, ymax), (zmin, zmax)) = (3, 2).".format(limits.shape))
        return 1

    print("Now calculating the electrostatic potential... ", end="", flush=True)
    _ts = time.time()

    _mesh_data = numerical_vars["full mesh"]
    _n_fun_coeff = numerical_vars["n_fun_coeff"]

    _mesh = bempp.api.grid.grid_from_element_data(_mesh_data["verts"],
                                                  _mesh_data["elems"],
                                                  _mesh_data["domns"])

    dp0_space = bempp.api.function_space(_mesh, "DP", 0)

    if numerical_vars["solution"] is None:
        n_fun = bempp.api.GridFunction(dp0_space, coefficients=_n_fun_coeff)
    else:
        n_fun = numerical_vars["solution"]

    # noinspection PyUnresolvedReferences
    all_vert = _mesh_data["verts"]

    # get limits from electrodes
    limits_elec = np.array([[np.min(all_vert[i, :]), np.max(all_vert[i, :])] for i in XYZ])

    # replace None limits with electrode limits
    limits[np.where(limits is None)] = limits_elec[np.where(limits is None)]

    res = np.array([res]).ravel()
    _n = np.array(np.round((limits[:, 1] - limits[:, 0]) / res, 10), int) + 1

    # Recalculate resolution to match integer n's
    _d = (limits[:, 1] - limits[:, 0]) / (_n - 1)

    # Generate a full spatial mesh to be indexed later
    _r = np.array([np.linspace(limits[i, 0], limits[i, 1], _n[i]) for i in XYZ])
    mesh = np.meshgrid(_r[X], _r[Y], _r[Z], indexing='ij')  # type: np.ndarray

    # Initialize potential array
    pot = np.zeros(mesh[0].shape)

    # Index borders (can be float)
    borders = np.array([np.linspace(0, _n[i], domain_decomp[i] + 1) for i in XYZ])

    # Indices (must be int)
    # note: rounding will likely lead to domains that are off in size by one index, but that's fine
    start_idxs = np.array([np.array(borders[i][:-1], int) - overlap for i in XYZ])
    end_idxs = np.array([np.array(borders[i][1:], int) + overlap for i in XYZ])

    for i in XYZ:
        start_idxs[i][0] = 0
        end_idxs[i][-1] = int(borders[i][-1])

    _ts = time.time()

    # Iterate over all the dimensions, calculate the subset of potential

    domain_idx = 1
    for x1, x2 in zip(start_idxs[X], end_idxs[X]):
        for y1, y2 in zip(start_idxs[Y], end_idxs[Y]):
            for z1, z2 in zip(start_idxs[Z], end_idxs[Z]):
                grid_pts = np.vstack([_mesh[x1:x2, y1:y2, z1:z2].ravel() for _mesh in mesh])
                grid_pts_len = grid_pts.shape[1]  # save shape for later

                # temp_pot = bempp.api.operators.potential.laplace.single_layer(dp0_space, grid_pts) * n_fun
                slp_pot = bempp.api.operators.potential.laplace.single_layer(dp0_space, grid_pts)
                temp_pot = slp_pot * n_fun

                # Create array of original shape and fill with result at right place,
                # then move into master array
                _pot = np.zeros(grid_pts_len)
                _pot = temp_pot[0]

                pot[x1:x2, y1:y2, z1:z2] = _pot.reshape([x2 - x1, y2 - y1, z2 - z1])

                domain_idx += 1

    try:

        del grid_pts
        del _pot
        del temp_pot

    except Exception as _e:

        print("Exception {} happened, but trying to carry on...".format(_e))

    # TODO: Distribute results to other nodes -DW

    numerical_vars["ef_phi"] = pot
    numerical_vars["d"] = _d
    numerical_vars["n"] = _n
    numerical_vars["limits"] = limits

    si.numerical_variables = numerical_vars

    print("Done!", flush=True)

    return 0


def solve_bempp(electrode_assm):

    assert HAVE_BEMPP, "BEMPP not found. Aborting!"

    numerical_vars = electrode_assm.numerical_variables
    bempp_params = electrode_assm.numerical_parameters

    try:
        electrodes = numerical_vars["objects"].electrodes
    except KeyError:
        electrodes = electrode_assm.electrodes

    gmres_tol = bempp_params["gmres_tol"]

    if numerical_vars["full mesh"] is None:
        print("Please generate a mesh before solving with BEM++!")
        return 1

    print("Generating necessary BEM++ operators, function spaces. Solving... ", end="", flush=True)

    _mesh_data = numerical_vars["full mesh"]

    _mesh = bempp.api.grid.grid_from_element_data(_mesh_data["verts"],
                                                  _mesh_data["elems"],
                                                  _mesh_data["domns"])

    dp0_space = bempp.api.function_space(_mesh, "DP", 0)

    slp = bempp.api.operators.boundary.laplace.single_layer(dp0_space, dp0_space, dp0_space)

    domain_mapping = {}
    for name, electrode in electrodes.items():
        domain_mapping[electrode.bempp_domain] = electrode.voltage

    def f(*args):

        domain_index = args[2]
        result = args[3]

        result[0] = domain_mapping[domain_index]

    dirichlet_fun = bempp.api.GridFunction(dp0_space, fun=f)
    numerical_vars["grid_fun"] = dirichlet_fun
    numerical_vars["d_fun_coeff"] = dirichlet_fun.coefficients

    if electrode_assm.debug:
        dirichlet_fun.plot()

    sol, info, res = bempp.api.linalg.gmres(slp, dirichlet_fun, tol=gmres_tol, return_residuals=True)
    print("Done!", flush=True)

    # Save results
    numerical_vars["solution"] = sol
    numerical_vars["d_fun"] = dirichlet_fun
    numerical_vars["n_fun_coeff"] = sol.coefficients
    numerical_vars["f_space"] = dp0_space
    numerical_vars["operator"] = slp

    electrode_assm.numerical_variables = numerical_vars

    return 0


def solve_fenics(si):

    assert HAVE_FENICS, "Fenics not found. Aborting!"

    print("Solving with FEniCS... ", end="", flush=True)

    numerical_vars = si.numerical_variables
    full_mesh = numerical_vars["full_mesh"]

    mesh = full_mesh[0]
    markers = full_mesh[1]
    boundaries = full_mesh[2]

    dx = fn.Measure('dx', domain=mesh, subdomain_data=markers)
    V = fn.FunctionSpace(mesh, 'P', 1)

    volt = si.get_parameter("volt")

    # Build boundary conditions
    # TODO: Get the #'s from the meshing stuff -PW
    bcs = []
    for i in range(1, 1000):  # Anode
        bcs.append(fn.DirichletBC(V, fn.Constant(volt), boundaries, i))
    for i in range(1000, 2000):  # Cathode
        bcs.append(fn.DirichletBC(V, fn.Constant(-volt), boundaries, i))

    # TODO: Assuming voltages on apertures, housing, etc. are all zero -PW
    # for i in range(2000, 3000):  # Exit Aperture/Housing
    #     bcs.append(fn.DirichletBC(V, fn.Constant(0.0), boundaries, i))
    # for i in range(3000, 4000):  # Entrance Aperture
    #     bcs.append(fn.DirichletBC(V, fn.Constant(0.0), boundaries, i))
    # for i in range(4000, 5000):  # Cylinder/Boundary
    #     bcs.append(fn.DirichletBC(V, fn.Constant(0.0), boundaries, i))

    # Test and trial functions
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    a = fn.dot(fn.grad(u), fn.grad(v)) * dx
    L = fn.Constant('0.0') * v * dx

    u = fn.Function(V)

    _tsi = time.time()
    fn.solve(a == L, u, bcs, solver_parameters={"linear_solver": "cg", "preconditioner": "ilu"})
    _tsf = time.time()

    print("Done!", flush=True)

    print("Potential solving time: {:.4f} s".format(_tsf - _tsi), flush=True)

    numerical_vars["ef_phi"] = u
    numerical_vars["f_space"] = V

    si.numerical_variables = numerical_vars

    if si.debug:
        potentialFile = fn.File(TEMP_DIR + '/potential.pvd')
        potentialFile << u

        meshfile = fn.File(TEMP_DIR + '/mesh.pvd')
        meshfile << mesh

    return 0
