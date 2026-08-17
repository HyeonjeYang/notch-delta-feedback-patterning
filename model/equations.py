"""Exact deterministic D-N-A-S drift equations, copied unmodified from the
source repository's
``src/model1_notch_delta_stochastic/literature_ver2.py::_drift_ver2``
(``model_variant="modern"`` -- the only variant used anywhere in the final
analysis). No equation, sign, or coefficient was changed.

    D_dot = lam_D * h_D(S) - d_D*D - f_D*D*s_N - k_cis*N*D
    N_dot = lam_N          - d_N*N - f_N*N*s_D - k_cis*N*D
    A_dot = -d_A*A + f_N*N*s_D
    S_dot = -d_S*S + lam_S*h_S(A)

    h_D(S) = 1 / (1 + (Kd/S)^n_d)         [proneural -> Delta, increasing in S]
    h_S(A) = 1 / (1 + (A/Kd)^n_ASC)       [NICD -| proneural, decreasing in A]

    s_D = sum_{j~i} D_j,   s_N = sum_{j~i} N_j   (raw six-neighbour contact sums)
"""

from __future__ import annotations

import numpy as np

from model.params import NotchDeltaParams


def hill_D(S: np.ndarray, p: NotchDeltaParams) -> np.ndarray:
    S_safe = np.maximum(S, 1e-8)
    return 1.0 / (1.0 + (p.Kd / S_safe) ** p.n_d)


def hill_ASC(A: np.ndarray, p: NotchDeltaParams) -> np.ndarray:
    a_ratio = np.minimum(np.asarray(A) / p.Kd, 1e5)
    return 1.0 / (1.0 + a_ratio ** p.n_ASC)


def dhillD_dS(S: float, p: NotchDeltaParams) -> float:
    """Analytic derivative of hill_D wrt S (S > 0)."""
    u = (p.Kd / S) ** p.n_d
    h = 1.0 / (1.0 + u)
    return p.n_d * u / S * h ** 2


def dhillS_dA(A: float, p: NotchDeltaParams) -> float:
    """Analytic derivative of hill_ASC wrt A (A >= 0)."""
    if A <= 0.0:
        return 0.0
    v = (A / p.Kd) ** p.n_ASC
    h = 1.0 / (1.0 + v)
    return -p.n_ASC * v / A * h ** 2


def drift(D: np.ndarray, N: np.ndarray, A: np.ndarray, S: np.ndarray,
         p: NotchDeltaParams, s_D: np.ndarray, s_N: np.ndarray):
    """Return (D_dot, N_dot, A_dot, S_dot) given local fields and neighbour
    sums s_D, s_N (see geometry.lattice.hex_adj)."""
    h_D = hill_D(S, p)
    h_S = hill_ASC(A, p)
    D_dot = p.lam_D * h_D - p.d_D * D - p.f_D * D * s_N - p.k_cis * N * D
    N_dot = p.lam_N - p.d_N * N - p.f_N * N * s_D - p.k_cis * N * D
    A_dot = -p.d_A * A + p.f_N * N * s_D
    S_dot = -p.d_ASC * S + p.lam_ASC * h_S
    return D_dot, N_dot, A_dot, S_dot
