"""Mode-resolved spatial stability screening: for a given homogeneous fixed
point, scan every nonuniform graph eigenvalue and find the fastest-growing
mode. Logic copied unmodified from the source repository's
``patternability_map/scripts/pm_core.py`` (``lambda_pattern``).

This is the primary screening quantity used throughout the final analysis:

    lambda_pattern(p) = max over nonuniform graph modes mu
                         of [ max eigenvalue real part of J_mu ]

lambda_pattern < 0: homogeneous state linearly stable.
lambda_pattern = 0: candidate bifurcation boundary.
lambda_pattern > 0: linearly unstable to a spatial perturbation.
"""

from __future__ import annotations

import numpy as np

from geometry.graph_spectrum import graph_eigenvalues
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams

_ALL_EIGS = graph_eigenvalues()
_NONUNIFORM_EIGS = _ALL_EIGS[~np.isclose(_ALL_EIGS, 6.0, atol=1e-9)]
UNIFORM_MU = 6.0
DOMINANT_PATTERN_MU = -3.0  # the three-sublattice mode, found dominant at every tested point


def lambda_pattern(fixed_point: np.ndarray, p: NotchDeltaParams, eigs: np.ndarray | None = None) -> dict:
    """Mode-resolved spatial stability at one homogeneous fixed point."""
    if eigs is None:
        eigs = _NONUNIFORM_EIGS
    best_re, best_mu, best_im = -np.inf, None, 0.0
    for mu in eigs:
        eigvals = np.linalg.eigvals(jacobian_block(mu, fixed_point, p))
        idx = int(np.argmax(eigvals.real))
        if eigvals[idx].real > best_re:
            best_re, best_mu, best_im = float(eigvals[idx].real), float(mu), float(eigvals[idx].imag)
    uniform_re = float(np.max(np.linalg.eigvals(jacobian_block(UNIFORM_MU, fixed_point, p)).real))
    return {
        "nonuniform_max_re_lambda": best_re,
        "nonuniform_best_mu": best_mu,
        "nonuniform_oscillatory": bool(abs(best_im) > 1e-9),
        "uniform_re_lambda": uniform_re,
    }
