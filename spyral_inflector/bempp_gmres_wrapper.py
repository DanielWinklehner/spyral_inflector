"""Iterative solver interfaces with GPU acceleration via CuPy."""

import numpy as _np

# Try to import CuPy for GPU acceleration
try:
    import cupy as _cp
    import cupyx.scipy.sparse.linalg as _cu_linalg
    CUPY_VERSION = tuple(int(x) for x in _cp.__version__.split('.')[:2])
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_VERSION = None
    CUPY_AVAILABLE = False


# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals


class IterationCounter(object):
    """Iteration Counter class."""

    def __init__(self, store_residuals, iteration_is_cg=False, operator=None, rhs=None):
        self._count = 0
        self._store_residuals = store_residuals
        self._residuals = []
        self._iteration_is_cg = iteration_is_cg
        self._operator = operator
        self._rhs = rhs

    def __call__(self, x):
        """Call."""
        from bempp_cl.api import log

        self._count += 1
        if self._store_residuals:
            if self._iteration_is_cg:
                res = self._rhs - self._operator * x
            else:
                res = x
            self._residuals.append(_np.linalg.norm(res))
            log(f"GMRES Iteration {self._count} with residual {self._residuals[-1]}")
        else:
            log(f"GMRES Iteration {self._count}")

    @property
    def count(self):
        """Return the number of iterations."""
        return self._count

    @property
    def residuals(self):
        """Return the vector of residuals."""
        return self._residuals


def _dense_array(A_op):
    """Return the dense matrix behind a bempp discrete operator, or None.

    Only dense assembly exposes one. FMM and H-matrix operators do not, and both
    the GPU path and the Jacobi preconditioner need it.
    """
    for attr in ("_A", "A"):
        if hasattr(A_op, attr):
            try:
                dense = _np.asarray(getattr(A_op, attr))
            except Exception:
                return None
            return dense if dense.ndim == 2 else None

    return None


def _jacobi_preconditioner(A_op, on_gpu):
    """Diagonal (Jacobi) preconditioner for a dense operator, or None.

    The Laplace single-layer matrix is badly scaled -- on a 15.5k element spiral
    inflector mesh its diagonal spans 8.0e-11 to 1.7e-07 -- so dividing the
    diagonal out cuts the iteration count sharply on both backends.
    """
    dense = _dense_array(A_op)

    if dense is None:
        return None

    diag = _np.diag(dense).copy()

    if not _np.all(_np.isfinite(diag)) or _np.any(diag == 0.0):
        return None

    if on_gpu:
        diag_gpu = _cp.asarray(diag)

        return _cu_linalg.LinearOperator(dense.shape,
                                         matvec=lambda v: v / diag_gpu,
                                         dtype=_cp.float64)

    import scipy.sparse.linalg

    return scipy.sparse.linalg.LinearOperator(dense.shape, matvec=lambda v: v / diag)


def _resolve_preconditioner(name, A_op, on_gpu):
    """Turn the preconditioner setting into a LinearOperator, or None.

    "auto" means jacobi on both backends. Measured on a 15,508 element spiral
    inflector mesh at restart=200, seconds for the GMRES solve:

        tol     GPU none   GPU jacobi   CPU none   CPU jacobi
        1e-05      2.71        0.90        9.84        2.80
        1e-06      3.10        0.92       15.70        3.24
        1e-08      3.80        0.90       25.60        4.35
        1e-10      5.73        2.20           -        5.09

    It is also far more accurate, which is the stronger reason to keep it on.
    Against a dense direct solve, relative error in the surface-charge
    coefficients and in the potential sampled along the beam axis:

        tol      none: coeff / potential      jacobi: coeff / potential
        1e-05      1.97e-03 / 8.19e-07          6.14e-09 / 8.04e-12
        1e-06      1.01e-05 / 4.39e-09          6.14e-09 / 8.04e-12
        1e-08      8.63e-08 / 3.19e-11          6.14e-09 / 8.04e-12

    Note the jacobi column is flat: with restart=200 it converges to ~1e-9 in
    the first restart cycle whatever tol says, so gmres_tol stops being the
    thing that controls accuracy. Unpreconditioned at the default tol=1e-5 the
    coefficients are only good to 2e-3.
    """
    import bempp_cl.api

    name = "none" if name is None else str(name).lower()

    if name == "auto":
        name = "jacobi"

    if name == "none":
        return None

    if name != "jacobi":
        raise ValueError("Unknown preconditioner %r; use 'auto', 'jacobi' or 'none'" % name)

    preconditioner = _jacobi_preconditioner(A_op, on_gpu)

    if preconditioner is None:
        bempp_cl.api.log("Jacobi preconditioner unavailable (no usable dense diagonal); "
                         "continuing without one")

    return preconditioner


def gmres(
        A,
        b,
        tol=1e-5,
        restart=None,
        maxiter=None,
        use_strong_form=False,
        return_residuals=False,
        return_iteration_count=False,
        use_gpu=True,
        preconditioner="auto",
):
    """Perform GMRES solve via CuPy (GPU) or scipy (CPU).

    This function behaves like the scipy.sparse.linalg.gmres function. But
    instead of a linear operator and a vector b it takes a boundary operator
    and a grid function or a blocked operator and a list of grid functions.
    The result is returned as a grid function or as a list of grid functions
    in the correct spaces.

    Parameters
    ----------
    use_gpu : bool
        If True and CuPy is available, use GPU acceleration. Default: True.
    restart : int or None
        GMRES restart depth. None uses the backend default of 20, which is far
        from optimal for these operators: on a 15,508 element spiral inflector
        mesh, restart=200 took the GPU solve from 7.6 s to 2.1 s and the CPU
        solve from 60.6 s to 9.8 s.
    preconditioner : str
        'auto' (jacobi), 'jacobi' or 'none'. Jacobi is roughly 3-4x faster on
        both backends here and needs a dense operator; it is skipped with a log
        message otherwise.

    """
    from bempp_cl.api.assembly.boundary_operator import BoundaryOperator
    from bempp_cl.api.assembly.blocked_operator import BlockedOperatorBase

    if isinstance(A, BoundaryOperator):
        return _gmres_single_op_imp(
            A,
            b,
            tol,
            restart,
            maxiter,
            use_strong_form,
            return_residuals,
            return_iteration_count,
            use_gpu,
            preconditioner,
        )

    if isinstance(A, BlockedOperatorBase):
        return _gmres_block_op_imp(
            A,
            b,
            tol,
            restart,
            maxiter,
            use_strong_form,
            return_residuals,
            return_iteration_count,
            use_gpu,
            preconditioner,
        )

    raise ValueError("A must be a BoundaryOperator or BlockedBoundaryOperator")


def _gmres_single_op_imp(
        A,
        b,
        tol=1e-5,
        restart=None,
        maxiter=None,
        use_strong_form=False,
        return_residuals=False,
        return_iteration_count=False,
        use_gpu=True,
        preconditioner="auto",
):
    """Run implementation of GMRES for single operators (CPU or GPU).

    restart is the GMRES restart depth. Leaving it None uses the backend default
    of 20, which is far from optimal here: on a 15.5k element mesh, restart=200
    took the GPU solve from 7.6 s to 2.1 s and the CPU solve from 60.6 s to 9.8 s.

    preconditioner is 'auto' (jacobi), 'jacobi' or 'none'.
    """
    from bempp_cl.api.assembly.grid_function import GridFunction
    import scipy.sparse.linalg
    import bempp_cl.api
    import time

    if not isinstance(b, GridFunction):
        raise ValueError("b must be of type GridFunction")

    # Assemble weak form
    if use_strong_form:
        if not A.range.is_compatible(b.space):
            raise ValueError(
                "The range of A and the domain of A must have"
                + "the same number of unknowns if the strong form is used."
            )
        A_op = A.strong_form()
        b_vec = b.coefficients
    else:
        A_op = A.weak_form()
        b_vec = b.projections(A.dual_to_range)

    # Decide whether to use GPU
    use_gpu_actual = use_gpu and CUPY_AVAILABLE

    M = _resolve_preconditioner(preconditioner, A_op, use_gpu_actual)

    if use_gpu_actual:
        bempp_cl.api.log("Starting GMRES iteration on GPU (CuPy)")
        x, info, res = _gmres_gpu(A_op, b_vec, tol, restart, maxiter, return_residuals, M=M)
    else:
        if use_gpu and not CUPY_AVAILABLE:
            bempp_cl.api.log("CuPy not available, falling back to CPU scipy")
        bempp_cl.api.log("Starting GMRES iteration on CPU (scipy)")

        callback = IterationCounter(return_residuals)
        start_time = time.time()
        x, info = scipy.sparse.linalg.gmres(
            A_op, b_vec, rtol=tol, restart=restart, maxiter=maxiter, M=M, callback=callback
        )
        end_time = time.time()
        bempp_cl.api.log("GMRES finished in %i iterations and took %.2E sec." % (callback.count, end_time - start_time))

        res = callback.residuals if return_residuals else None

    res_fun = GridFunction(A.domain, coefficients=x.ravel())

    if return_residuals and return_iteration_count:
        iteration_count = len(res) if res else 0
        return res_fun, info, res, iteration_count

    if return_residuals:
        return res_fun, info, res

    if return_iteration_count:
        iteration_count = len(res) if res else 0
        return res_fun, info, iteration_count

    return res_fun, info


def _gmres_block_op_imp(
        A,
        b,
        tol=1e-5,
        restart=None,
        maxiter=None,
        use_strong_form=False,
        return_residuals=False,
        return_iteration_count=False,
        use_gpu=True,
        preconditioner="auto",
):
    """Run implementation of GMRES for blocked operators (CPU or GPU)."""
    import scipy.sparse.linalg
    import bempp_cl.api
    import time
    from bempp_cl.api.assembly.blocked_operator import (
        coefficients_from_grid_functions_list,
        projections_from_grid_functions_list,
        grid_function_list_from_coefficients,
    )

    # Assemble weak form
    if use_strong_form:
        b_vec = coefficients_from_grid_functions_list(b)
        A_op = A.strong_form()
    else:
        A_op = A.weak_form()
        b_vec = projections_from_grid_functions_list(b, A.dual_to_range_spaces)

    # Decide whether to use GPU
    use_gpu_actual = use_gpu and CUPY_AVAILABLE

    M = _resolve_preconditioner(preconditioner, A_op, use_gpu_actual)

    if use_gpu_actual:
        bempp_cl.api.log("Starting GMRES iteration on GPU (CuPy)")
        x, info, res = _gmres_gpu(A_op, b_vec, tol, restart, maxiter, return_residuals, M=M)
    else:
        if use_gpu and not CUPY_AVAILABLE:
            bempp_cl.api.log("CuPy not available, falling back to CPU scipy")
        bempp_cl.api.log("Starting GMRES iteration on CPU (scipy)")

        callback = IterationCounter(return_residuals)
        start_time = time.time()
        x, info = scipy.sparse.linalg.gmres(
            A_op, b_vec, rtol=tol, restart=restart, maxiter=maxiter, M=M, callback=callback
        )
        end_time = time.time()
        bempp_cl.api.log("GMRES finished in %i iterations and took %.2E sec." % (callback.count, end_time - start_time))

        res = callback.residuals if return_residuals else None

    res_fun = grid_function_list_from_coefficients(x.ravel(), A.domain_spaces)

    if return_residuals and return_iteration_count:
        iteration_count = len(res) if res else 0
        return res_fun, info, res, iteration_count

    if return_residuals:
        return res_fun, info, res

    if return_iteration_count:
        iteration_count = len(res) if res else 0
        return res_fun, info, iteration_count

    return res_fun, info


def _gmres_gpu(A_op, b_vec, tol, restart, maxiter, return_residuals, M=None):
    """Perform GMRES on GPU using CuPy.

    Requires a dense operator: the whole matrix is copied to the device, which is
    O(N^2) memory (1.92 GB at N = 15,508). FMM or H-matrix operators cannot use
    this path.
    """
    import bempp_cl.api
    import time

    if not CUPY_AVAILABLE:
        raise RuntimeError("CuPy is not available")

    # Transfer to GPU
    start_time = time.time()

    # Extract dense matrix from bempp operator
    A_dense = _dense_array(A_op)

    if A_dense is None:
        raise RuntimeError(f"Could not extract a dense matrix from {type(A_op)}")

    # Convert to CuPy
    A_gpu = _cp.asarray(A_dense, dtype=_cp.float64)
    b_gpu = _cp.asarray(b_vec, dtype=_cp.float64)
    transfer_time_to = time.time() - start_time
    bempp_cl.api.log(f"Transferred data to GPU in {transfer_time_to:.4f}s")

    # Track residuals (CuPy callback receives residual norm, not solution vector)
    residuals = []
    iteration_count = [0]

    def gpu_callback(residual_norm):
        """CuPy GMRES callback - receives residual norm (scalar)"""
        iteration_count[0] += 1
        if return_residuals:
            residuals.append(float(residual_norm))
            bempp_cl.api.log(f"GMRES Iteration {iteration_count[0]} with residual {residual_norm:.6e}")
        else:
            bempp_cl.api.log(f"GMRES Iteration {iteration_count[0]}")

    # Solve on GPU
    start_time = time.time()

    if CUPY_VERSION >= (14, 0):
        x_gpu, info = _cu_linalg.gmres(
            A_gpu,
            b_gpu,
            rtol=tol,
            atol=0.0,
            restart=restart,
            maxiter=maxiter,
            M=M,
            callback=gpu_callback,
        )
    else:
        x_gpu, info = _cu_linalg.gmres(
            A_gpu,
            b_gpu,
            tol=tol,
            restart=restart,
            maxiter=maxiter,
            M=M,
            callback=gpu_callback,
        )

    solve_time = time.time() - start_time
    bempp_cl.api.log(f"GPU GMRES finished in {iteration_count[0]} iterations and took {solve_time:.4E}s")

    # Transfer back to CPU
    start_time = time.time()
    x = _cp.asnumpy(x_gpu)
    transfer_time_from = time.time() - start_time
    bempp_cl.api.log(f"Transferred solution to CPU in {transfer_time_from:.4f}s")

    return x, info, residuals if return_residuals else None