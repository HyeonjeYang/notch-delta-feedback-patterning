"""Mode-resolved 4x4 Jacobian block, copied (logic unchanged) from the
source repository's
``continuous_time_validation/scripts/linear_stability.py::jacobian_block``.

Because the contact graph is regular and its adjacency operator W is
diagonalized by its own eigenmodes (W v_mu = mu v_mu), a perturbation
delta_x_i = v_mu,i * xi_mu(t) decouples the full multicellular
linearization into one 4x4 ODE per graph eigenvalue mu:

    d xi_mu / dt = J_mu xi_mu

State order is (D, N, A, S). This construction was independently verified
against a full finite-difference Jacobian of the 2304-dim vectorized
lattice RHS to 6e-10 (see docs/GEOMETRY_DERIVATION.md /
docs/CONTROL_THEORY_DERIVATION.md), and against the exact Schur-complement
characteristic polynomial to 1.4e-13 (control/schur_verification.py).
"""

from __future__ import annotations

import numpy as np

from model.equations import dhillD_dS, dhillS_dA
from model.params import NotchDeltaParams


def jacobian_block(mu: float, fixed_point: np.ndarray, p: NotchDeltaParams) -> np.ndarray:
    """4x4 Jacobian block (D,N,A,S order) for graph eigenvalue mu at the
    given homogeneous fixed point fp=(D*,N*,A*,S*)."""
    D, N, A, S = fixed_point
    J = np.zeros((4, 4))
    # row D: D_dot = lam_D*h_D(S) - d_D*D - f_D*D*s_N - k_cis*N*D
    J[0, 0] = -p.d_D - p.f_D * 6.0 * N - p.k_cis * N              # dD_dot/dD
    J[0, 1] = -p.k_cis * D - p.f_D * D * mu                        # dD_dot/dN
    J[0, 2] = 0.0
    J[0, 3] = p.lam_D * dhillD_dS(S, p)                            # dD_dot/dS
    # row N: N_dot = lam_N - d_N*N - f_N*N*s_D - k_cis*N*D
    J[1, 0] = -p.k_cis * N - p.f_N * N * mu                        # dN_dot/dD
    J[1, 1] = -p.d_N - p.f_N * 6.0 * D - p.k_cis * D               # dN_dot/dN
    J[1, 2] = 0.0
    J[1, 3] = 0.0
    # row A: A_dot = -d_A*A + f_N*N*s_D
    J[2, 0] = p.f_N * N * mu                                       # dA_dot/dD
    J[2, 1] = p.f_N * 6.0 * D                                      # dA_dot/dN
    J[2, 2] = -p.d_A                                               # dA_dot/dA
    J[2, 3] = 0.0
    # row S: S_dot = -d_S*S + lam_S*h_S(A)
    J[3, 0] = 0.0
    J[3, 1] = 0.0
    J[3, 2] = p.lam_ASC * dhillS_dA(A, p)                          # dS_dot/dA
    J[3, 3] = -p.d_ASC                                             # dS_dot/dS
    return J
