"""Reciprocal-space growth-rate map: lambda_max(k) over the full graph
spectrum at a given parameter point, interpreting the tissue as a
mode-selective spatial filter. Logic copied unmodified from the source
repository's ``preprint_final_analysis/scripts/reciprocal_space.py``.
"""

from __future__ import annotations

import numpy as np

from geometry.graph_spectrum import ROWS, COLS, circulant_eigenvalue
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams


def growth_rate_map(fixed_point: np.ndarray, p: NotchDeltaParams,
                    rows: int = ROWS, cols: int = COLS) -> np.ndarray:
    """(rows, cols) array of lambda_max(k1,k2) = max Re(eig(J_mu(k))),
    indexed by (m, n)."""
    out = np.zeros((rows, cols))
    for m in range(rows):
        for n in range(cols):
            mu = circulant_eigenvalue(m, n, rows, cols)
            J_mu = jacobian_block(mu, fixed_point, p)
            out[m, n] = float(np.max(np.linalg.eigvals(J_mu).real))
    return out
