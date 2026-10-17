"""Graph input-output formulation of the same linearisation, used for the
comparison with existing graph-based lateral-inhibition criteria (e.g.
Arcak 2013, IEEE TAC 58:1250; Rufino Ferreira & Arcak 2013, SIAM J Appl
Dyn Syst 12:2012).

Every cell is a 4-state system driven by its neighbour sums
u = (s_D, s_N) and emitting y = (D, N). Linearised at the homogeneous
equilibrium,

    x_dot = J_loc x + B u,    y = C x,    u = W y  ->  J_mu = J_loc + mu B C,

where J_loc is J_mu with every mu-proportional (neighbour) entry removed.
When J_loc is nonsingular the local static input-output gain is
K = -C J_loc^{-1} B (2x2), and Sylvester's determinant identity gives

    det J_mu = det(J_loc) * det(I - mu K).

With a single coupling channel (f_D = 0, so only s_D matters), K is a scalar
and the stationary condition mu K = 1 is the static-gain graph condition
used in those works. Here trans-depletion f_D D s_N adds a second channel,
so K is 2x2. The zero set {det(I - mu K) = 0} equals {G(mu) = 1}
(control/state_space.py). This determinant tests only the zero-frequency loop;
the cited graph-based theorems have additional assumptions and conclusions.
"""

from __future__ import annotations

import numpy as np

from model.jacobian import jacobian_block
from model.params import NotchDeltaParams


def coupling_decomposition(fixed_point: np.ndarray, p: NotchDeltaParams):
    """(J_loc, B, C) with J_mu = J_loc + mu * B @ C."""
    D, N, _, _ = fixed_point
    J_loc = jacobian_block(0.0, fixed_point, p)
    B = np.zeros((4, 2))
    B[1, 0] = -p.f_N * N   # s_D -> N_dot  (trans depletion of receptor)
    B[2, 0] = p.f_N * N    # s_D -> A_dot  (trans activation)
    B[0, 1] = -p.f_D * D   # s_N -> D_dot  (trans depletion of ligand)
    C = np.zeros((2, 4))
    C[0, 0] = 1.0          # y_1 = D
    C[1, 1] = 1.0          # y_2 = N
    return J_loc, B, C


def static_io_gain(fixed_point: np.ndarray, p: NotchDeltaParams) -> np.ndarray:
    """K = -C J_loc^{-1} B: static response of (D, N) to neighbour sums."""
    J_loc, B, C = coupling_decomposition(fixed_point, p)
    return -C @ np.linalg.solve(J_loc, B)


def stationary_margin(mu: float, fixed_point: np.ndarray, p: NotchDeltaParams) -> float:
    """det(I - mu K); changes sign exactly where G(mu) crosses 1."""
    return float(np.linalg.det(np.eye(2) - mu * static_io_gain(fixed_point, p)))
