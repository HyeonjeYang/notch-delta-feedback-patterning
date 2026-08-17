"""Homogeneous fixed-point solver, copied (logic unchanged) from the source
repository's ``patternability_map/scripts/pm_core.py`` (``rhs_single_cell``,
``find_fixed_points``). Exploits the regular six-neighbour graph: at a
homogeneous state every cell's neighbour sum equals 6x its own value.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import root

from model.equations import drift
from model.params import NotchDeltaParams

_GUESS_SCALES = [
    (0.4, 0.5, 0.0, 0.3),
    (0.1, 0.1, 0.1, 0.1),
    (2.0, 2.0, 2.0, 0.9),
    (0.01, 1.0, 0.01, 0.01),
    (5.0, 0.05, 5.0, 0.01),
    (30.0, 0.01, 1.0, 0.4),
]


def rhs_single_cell(y: np.ndarray, p: NotchDeltaParams, z: float = 6.0) -> np.ndarray:
    """Homogeneous single-cell reduction: neighbour sum = z * local value."""
    D, N, A, S = y
    D_a = np.array([[D]]); N_a = np.array([[N]]); A_a = np.array([[A]]); S_a = np.array([[S]])
    s_D = z * D_a
    s_N = z * N_a
    dD, dN, dA, dS = drift(D_a, N_a, A_a, S_a, p, s_D, s_N)
    return np.array([dD[0, 0], dN[0, 0], dA[0, 0], dS[0, 0]])


def find_fixed_points(p: NotchDeltaParams, tol: float = 1e-10,
                      extra_guesses: list | None = None) -> list[np.ndarray]:
    """Return every distinct nonnegative homogeneous fixed point found from
    a battery of initial guesses. Empty list if none found."""
    guesses = list(_GUESS_SCALES) + list(extra_guesses or [])
    found: list[np.ndarray] = []
    for guess in guesses:
        try:
            sol = root(lambda y: rhs_single_cell(y, p), np.array(guess, dtype=float),
                      method="hybr", tol=tol)
        except Exception:
            continue
        if not sol.success:
            continue
        y = sol.x
        if np.any(y < -1e-7) or not np.all(np.isfinite(y)):
            continue
        y = np.clip(y, 0, None)
        resid = np.max(np.abs(rhs_single_cell(y, p)))
        if resid > 1e-6 * max(1.0, np.max(np.abs(y))):
            continue
        if not any(np.allclose(y, u, atol=1e-6, rtol=1e-6) for u in found):
            found.append(y)
    return found
