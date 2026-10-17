"""Mode-resolved feedback partition of the 4x4 Jacobian block J_mu.

Partition x = (D, N, A, S) as z = (D, N, A) and the scalar S:

    z_dot = A_mu z + b S,      S_dot = c^T z + a_S S
    A_mu = J_mu[0:3, 0:3]   b = J_mu[0:3, 3]   c^T = J_mu[3, 0:3]   a_S = J_mu[3, 3]

    L_mu(s) = c^T (sI - A_mu)^{-1} b / (s - a_S)
    det(sI - J_mu) = det(sI - A_mu) (s - a_S) [1 - L_mu(s)]

This is the standard Schur-complement factorisation and is not new. Its value
at s = 0 gives the static (DC) loop gain G(mu) = L_mu(0) and, whenever A_mu is
nonsingular and a_S != 0,

    det J_mu = det(A_mu) * a_S * [1 - G(mu)].

In this model A_mu is Hurwitz for every graph eigenvalue mu in [-z, z], so
det(A_mu) < 0, and a_S = -d_S < 0 (proof: docs/THEORY.md). Hence

    G(mu) = 1  <=>  J_mu has a zero eigenvalue,
    G(mu) > 1  <=>  J_mu has an odd number of positive real eigenvalues.

G(mu) < 1 therefore does NOT imply stability. A complex pair, or two real
eigenvalues, in the right half-plane leaves det J_mu > 0. The complete
per-mode test is the Nyquist criterion: the open loop is stable, so the
number of right-half-plane eigenvalues of J_mu equals the number of
clockwise encirclements of +1 by L_mu(i omega).
"""

from __future__ import annotations

import numpy as np


def partition(J_mu: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Split the 4x4 (D,N,A,S) Jacobian block into (A_mu, b, c, a_S)."""
    return J_mu[0:3, 0:3], J_mu[0:3, 3], J_mu[3, 0:3], J_mu[3, 3]


def L_mu_of_s(s: complex, A_mu: np.ndarray, b: np.ndarray, c: np.ndarray, a_S: float) -> complex:
    return complex(c @ np.linalg.solve(s * np.eye(3) - A_mu, b)) / (s - a_S)


def G_exact(J_mu: np.ndarray) -> float:
    """Static (DC) loop gain G(mu) = L_mu(0)."""
    A_mu, b, c, a_S = partition(J_mu)
    return float(L_mu_of_s(0.0 + 0.0j, A_mu, b, c, a_S).real)


def stage_gains(J_mu: np.ndarray) -> dict[str, float]:
    """Factor G(mu) into the three static stage gains of the S->D->A-|S loop.

    kappa_D = lam_D h_D'(S*)                 S -> Delta production
    T_mu    = [A_mu^{-1}]_{A,D}              Delta production -> neighbour A,
                                             including the N-mediated
                                             cis/trans paths (mode dependent)
    kappa_S = lam_S |h_S'(A*)| / d_S         A -| S, static repression gain

    G(mu) = kappa_D * T_mu * kappa_S.
    """
    A_mu, b, c, a_S = partition(J_mu)
    kappa_D = float(b[0])
    kappa_S = float(-c[2] / -a_S)
    T_mu = float(np.linalg.inv(A_mu)[2, 0])
    return {"kappa_D": kappa_D, "T_mu": T_mu, "kappa_S": kappa_S, "G": kappa_D * T_mu * kappa_S}


def nyquist_curve(J_mu: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """L_mu(i omega) on a frequency grid."""
    A_mu, b, c, a_S = partition(J_mu)
    s = 1j * np.asarray(omega, dtype=float)
    M = s[:, None, None] * np.eye(3) - A_mu
    x = np.linalg.solve(M, np.broadcast_to(b.astype(complex), (len(s), 3))[..., None])[..., 0]
    return (x @ c) / (s - a_S)


def nyquist_rhp_count(J_mu: np.ndarray, n: int = 4000) -> int:
    """Number of right-half-plane closed-loop poles from the winding of
    1 - L_mu(i omega) about 0 (open loop stable, so P = 0). Used only as an
    independent check of the eigenvalue count away from marginal stability.
    As with any sampled winding computation, the frequency grid must resolve
    the curve; the all-sample audit checks it against direct eigenvalues."""
    scale = max(1.0, float(np.max(np.abs(np.linalg.eigvals(J_mu)))))
    w_pos = scale * np.logspace(-8, 4, n)
    omega = np.concatenate([-w_pos[::-1], [0.0], w_pos])
    f = 1.0 - nyquist_curve(J_mu, omega)
    winding = (np.unwrap(np.angle(f))[-1] - np.unwrap(np.angle(f))[0]) / (2 * np.pi)
    return int(round(-winding))
