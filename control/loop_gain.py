"""Exact DC loop gain and the reduced G_eff heuristic, copied (logic
unchanged) from the source repository's
``preprint_final_analysis/scripts/compute_exact_loop_gain.py`` and
``explain_geff.py``.

G_eff is the dominant-path (S->D->A->S) approximation of the exact
G_exact(mu) = L_mu(0): it drops the N-mediated coupling terms and
approximates [A_mu^{-1}]_{A,D} by J_AD(mu)/(J_DD*J_AA). Over the frozen
patternability-map global sample (n=3267 valid points), G_exact(-3)>1
classifies true spatial instability with 99.79% accuracy vs. 93.88% for
G_eff at its best threshold -- see docs/CONTROL_THEORY_DERIVATION.md
section 8. Reproducing that full comparison requires the large raw global
sample, which is intentionally not included in this repository (see
README.md "Reproducibility note"); the functions below are included so the
two gain quantities can be recomputed for any new fixed point.
"""

from __future__ import annotations

import numpy as np

from control.state_space import G_exact as _G_exact_dc  # re-exported below


def G_exact(J_mu: np.ndarray) -> float:
    """Exact DC loop gain L_mu(0), evaluated at mode mu's Jacobian block."""
    return _G_exact_dc(J_mu)


def G_eff(J_mu: np.ndarray) -> float:
    """Reduced dominant-path approximation of G_exact (drops the
    N-mediated coupling; normalizes by |restoring rates| instead of the
    signed diagonal)."""
    J_DD, J_AD, J_AA = J_mu[0, 0], J_mu[2, 0], J_mu[2, 2]
    J_DS, J_SA, J_SS = J_mu[0, 3], J_mu[3, 2], J_mu[3, 3]
    return float((J_DS * J_AD * J_SA) / (abs(J_DD) * abs(J_AA) * abs(J_SS)))
