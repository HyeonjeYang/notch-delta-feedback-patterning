"""Mode-resolved spatial stability screening.

For a homogeneous equilibrium, scan every graph eigenvalue mu of the contact
adjacency operator W (W v_mu = mu v_mu) and record the leading eigenvalue of
the 4x4 block J_mu:

    lambda_pattern = max over nonuniform mu of max Re eig(J_mu)

lambda_pattern < 0: the homogeneous state is linearly stable to all
spatial perturbations. lambda_pattern > 0: some spatial mode grows.

The 576 eigenvalues of the 24x24 lattice take only 54 distinct values, and
the scan runs over those.
"""

from __future__ import annotations

import numpy as np

from geometry.graph_spectrum import graph_eigenvalues
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams

UNIFORM_MU = 6.0
DOMINANT_PATTERN_MU = -3.0  # minimum graph eigenvalue: the three-sublattice (K-point) mode
UNIQUE_MU = np.unique(np.round(graph_eigenvalues(), 12))
NONUNIFORM_MU = UNIQUE_MU[UNIQUE_MU < UNIFORM_MU - 1e-9]


def leading_eigenvalue(M: np.ndarray) -> complex:
    eigvals = np.linalg.eigvals(M)
    return complex(eigvals[int(np.argmax(eigvals.real))])


def rhp_count(M: np.ndarray, tol: float = 0.0) -> int:
    """Number of eigenvalues with real part > tol."""
    return int(np.sum(np.linalg.eigvals(M).real > tol))


def lambda_pattern(fixed_point: np.ndarray, p: NotchDeltaParams, eigs: np.ndarray | None = None,
                   reduce=None) -> dict:
    """Mode-resolved spatial stability at one homogeneous equilibrium.

    ``reduce`` optionally maps J_mu to a reduced Jacobian (model/reduced.py)
    so that the same scan can be run on a matched reduced model."""
    if eigs is None:
        eigs = NONUNIFORM_MU
    transform = reduce if reduce is not None else (lambda J: J)
    best, best_mu = complex(-np.inf), None
    for mu in eigs:
        lead = leading_eigenvalue(transform(jacobian_block(mu, fixed_point, p)))
        if lead.real > best.real:
            best, best_mu = lead, float(mu)
    uniform = leading_eigenvalue(transform(jacobian_block(UNIFORM_MU, fixed_point, p)))
    return {
        "nonuniform_max_re_lambda": best.real,
        "nonuniform_best_mu": best_mu,
        "nonuniform_im_lambda": abs(best.imag),
        "nonuniform_oscillatory": bool(abs(best.imag) > 1e-9),
        "uniform_re_lambda": uniform.real,
        "uniform_im_lambda": abs(uniform.imag),
    }
