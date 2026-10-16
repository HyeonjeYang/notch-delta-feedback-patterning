"""Homogeneous equilibrium of the D-N-A-S model on a z-regular contact graph.

At a homogeneous state every neighbour sum equals z times the local value
(s_D = z D, s_N = z N; z = 6 on the triangular contact lattice). The
remaining four algebraic equations collapse to one scalar equation in D:

    N(D) = lam_N / (d_N + (z f_N + k_cis) D)              (decreasing in D)
    A(D) = z f_N N(D) D / d_A                              (increasing in D)
    S(D) = (lam_S / d_S) h_S(A(D))                         (decreasing in D)
    F(D) = lam_D h_D(S(D)) - D [d_D + (z f_D + k_cis) N(D)] = 0

The first term of F is non-increasing in D (h_D increases with S, S
decreases with D) and the second strictly increasing, so F is strictly
decreasing, with F(0) = lam_D h_D(lam_S/d_S) > 0 and
F(lam_D/d_D) < 0 (because h_D < 1). Hence, for every positive parameter set, there is
exactly one homogeneous equilibrium, it is positive, and it lies in
D in (0, lam_D/d_D). It is found here by Brent bracketing on that interval.

This replaces the multistart ``scipy.optimize.root`` search used in the v1
analysis (archive/v1_physbio_submission/model/fixed_point.py). That search
returned no root for 733 of the 4000 global-sample points. Those were
solver failures, not missing equilibria.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from model.equations import drift, hill_ASC, hill_D
from model.params import NotchDeltaParams

Z_CONTACT = 6.0


def _chain(D: float, p: NotchDeltaParams, z: float):
    N = p.lam_N / (p.d_N + (z * p.f_N + p.k_cis) * D)
    A = z * p.f_N * N * D / p.d_A
    S = p.lam_ASC / p.d_ASC * float(hill_ASC(np.asarray(A), p))
    return N, A, S


def _scalar_residual(D: float, p: NotchDeltaParams, z: float) -> float:
    N, _, S = _chain(D, p, z)
    return p.lam_D * float(hill_D(np.asarray(S), p)) - D * (p.d_D + (z * p.f_D + p.k_cis) * N)


def rhs_single_cell(y: np.ndarray, p: NotchDeltaParams, z: float = Z_CONTACT) -> np.ndarray:
    """Homogeneous single-cell reduction: neighbour sum = z * local value."""
    D, N, A, S = (np.array([[v]]) for v in y)
    dD, dN, dA, dS = drift(D, N, A, S, p, z * D, z * N)
    return np.array([dD[0, 0], dN[0, 0], dA[0, 0], dS[0, 0]])


def homogeneous_fixed_point(p: NotchDeltaParams, z: float = Z_CONTACT) -> np.ndarray:
    """The unique homogeneous equilibrium (D*, N*, A*, S*)."""
    D = brentq(_scalar_residual, 0.0, p.lam_D / p.d_D, args=(p, z),
               xtol=np.finfo(float).tiny, rtol=4 * np.finfo(float).eps, maxiter=500)
    N, A, S = _chain(D, p, z)
    return np.array([D, N, A, S])
