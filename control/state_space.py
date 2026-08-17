"""Exact mode-resolved state-space / feedback representation, copied
(logic unchanged) from the source repository's
``preprint_final_analysis/scripts/control_core.py``.

Partition x = (D, N, A, S) as z = (D, N, A)^T, S = ASC (scalar), read
directly off the existing 4x4 Jacobian block J_mu (model.jacobian) --
no re-derivation of the biology, just a different grouping of the same
matrix entries:

    z_dot = A_mu z + b S
    S_dot = c^T z + a_S S

    A_mu = J_mu[0:3, 0:3]      b = J_mu[0:3, 3]
    c^T  = J_mu[3, 0:3]        a_S = J_mu[3, 3]

Laplace-domain plant and loop transfer function:

    P_mu(s) = c^T (sI - A_mu)^{-1} b
    L_mu(s) = P_mu(s) / (s - a_S)

Exact closed-loop characteristic equation (Schur-complement identity,
verified in control/schur_verification.py to 1.4e-13):

    det(sI - J_mu) = det(sI - A_mu) * (s - a_S) * [ 1 - L_mu(s) ]

Sign convention (1 - L, not 1 + L) was derived from the block-determinant
expansion, not assumed -- see docs/CONTROL_THEORY_DERIVATION.md section 4.
"""

from __future__ import annotations

import numpy as np


def partition(J_mu: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Split the 4x4 (D,N,A,S) Jacobian block into (A_mu, b, c, a_S)."""
    A_mu = J_mu[0:3, 0:3]
    b = J_mu[0:3, 3]
    c = J_mu[3, 0:3]
    a_S = J_mu[3, 3]
    return A_mu, b, c, a_S


def P_mu_of_s(s: complex, A_mu: np.ndarray, b: np.ndarray, c: np.ndarray) -> complex:
    M = np.linalg.solve(s * np.eye(3) - A_mu, b)
    return complex(c @ M)


def L_mu_of_s(s: complex, A_mu: np.ndarray, b: np.ndarray, c: np.ndarray, a_S: float) -> complex:
    return P_mu_of_s(s, A_mu, b, c) / (s - a_S)


def G_exact(J_mu: np.ndarray) -> float:
    """Exact DC (zero-frequency) loop gain L_mu(0) for one mode's Jacobian block."""
    A_mu, b, c, a_S = partition(J_mu)
    return float(L_mu_of_s(0.0 + 0.0j, A_mu, b, c, a_S).real)
