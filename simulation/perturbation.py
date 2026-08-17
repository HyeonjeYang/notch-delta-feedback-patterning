"""Controlled small-amplitude perturbation along the dominant unstable
eigenmode, used for the nonlinear validation runs (Figure 5). Copied
(logic unchanged) from the source repository's
``preprint_final_analysis/scripts/rerun_p2_snapshots.py`` and
``patternability_map/scripts/phase6_nonlinear_validation.py``.

Deliberately NOT a random/stochastic initial condition (task requirement
throughout the source project): the perturbation direction is the exact
dominant eigenvector of the mode-resolved Jacobian at mu=-3, scaled by a
small relative amplitude, so the initial linear growth rate can be
compared quantitatively to the predicted pole.
"""

from __future__ import annotations

import numpy as np

from geometry.lattice import COLS, ROWS
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams
from stability.mode_scan import DOMINANT_PATTERN_MU

MODE_M, MODE_N = 8, 16  # the (m,n) indices realizing mu=-3 on the 24x24 grid


def dominant_eigvec(fixed_point: np.ndarray, p: NotchDeltaParams, mu: float = DOMINANT_PATTERN_MU) -> np.ndarray:
    J_mu = jacobian_block(mu, fixed_point, p)
    eigvals, eigvecs = np.linalg.eig(J_mu)
    vec = eigvecs[:, np.argmax(eigvals.real)].real
    return vec / np.linalg.norm(vec)


def perturbed_field(fixed_point: np.ndarray, eigvec: np.ndarray, rel_amp: float,
                    rows: int = ROWS, cols: int = COLS):
    """Homogeneous fixed point + rel_amp * eigvec * cos(mode (m,n)=(8,16) spatial pattern)."""
    D0, N0, A0, S0 = fixed_point
    r = np.arange(rows).reshape(-1, 1)
    c = np.arange(cols).reshape(1, -1)
    phi = 2 * np.pi * (MODE_M * r / rows + MODE_N * c / cols)
    spatial = np.cos(phi)
    scale = rel_amp * np.array([D0, max(N0, 1e-6), max(A0, 1e-6), S0])
    D = np.clip(D0 + scale[0] * eigvec[0] * spatial, 0, None)
    N = np.clip(N0 + scale[1] * eigvec[1] * spatial, 0, None)
    A = np.clip(A0 + scale[2] * eigvec[2] * spatial, 0, None)
    S = np.clip(S0 + scale[3] * eigvec[3] * spatial, 0, None)
    return D, N, A, S
