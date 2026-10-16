"""Matched reduced models obtained by quasi-steady-state (QSS) elimination
of the explicit downstream layer.

* ``DNA`` (3 variables): proneural S is slaved to activated Notch,
  S = S_q(A) = (lam_S/d_S) h_S(A), so Delta production becomes
  lam_D h_D(S_q(A)). This is a Sprinzak-type Notch/Delta/reporter model
  with a composite repression function.
* ``DN`` (2 variables): A and S are both slaved,
  A = f_N N s_D / d_A and S = S_q(A).

Both reductions keep every equilibrium of the full model, homogeneous or
patterned, because they only impose dS/dt = 0 (and dA/dt = 0). They also
keep every static (zero-frequency) response. Their mode-resolved Jacobians
are the Schur complements of J_mu with respect to the eliminated block.
Because

    det J_mu = det(J_ff) * det(J_red),   J_ff = eliminated block (Hurwitz),

the full and reduced models share their zero-eigenvalue (stationary
bifurcation) set exactly. They can differ in everything that depends on
dynamics: growth rates and instabilities that do not pass through a zero
eigenvalue.

``scale_proneural_timescale`` gives a one-parameter family that interpolates
between the full model (alpha = 1) and the DNA reduction (alpha -> inf) while
keeping every equilibrium and every static gain fixed.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from model.equations import hill_ASC, hill_D
from model.params import NotchDeltaParams


def reduce_S(J_mu: np.ndarray) -> np.ndarray:
    """3x3 Jacobian of the DNA reduction (Schur complement w.r.t. S)."""
    A_mu, b, c, a_S = J_mu[:3, :3], J_mu[:3, 3], J_mu[3, :3], J_mu[3, 3]
    return A_mu - np.outer(b, c) / a_S


def reduce_AS(J_mu: np.ndarray) -> np.ndarray:
    """2x2 Jacobian of the DN reduction (Schur complement w.r.t. (A, S))."""
    J_ss, J_sf, J_fs, J_ff = J_mu[:2, :2], J_mu[:2, 2:], J_mu[2:, :2], J_mu[2:, 2:]
    return J_ss - J_sf @ np.linalg.solve(J_ff, J_fs)


def scale_proneural_timescale(p: NotchDeltaParams, alpha: float) -> NotchDeltaParams:
    """(lam_S, d_S) -> alpha * (lam_S, d_S): S* and the static gain
    lam_S |h_S'(A*)| / d_S are unchanged; only the proneural relaxation
    rate d_S changes."""
    return dataclasses.replace(p, lam_ASC=alpha * p.lam_ASC, d_ASC=alpha * p.d_ASC)


def drift_DNA(D: np.ndarray, N: np.ndarray, A: np.ndarray, p: NotchDeltaParams,
              s_D: np.ndarray, s_N: np.ndarray):
    """Drift of the 3-variable DNA reduction (S at quasi-steady state)."""
    S_q = p.lam_ASC / p.d_ASC * hill_ASC(A, p)
    D_dot = p.lam_D * hill_D(S_q, p) - p.d_D * D - p.f_D * D * s_N - p.k_cis * N * D
    N_dot = p.lam_N - p.d_N * N - p.f_N * N * s_D - p.k_cis * N * D
    A_dot = -p.d_A * A + p.f_N * N * s_D
    return D_dot, N_dot, A_dot
