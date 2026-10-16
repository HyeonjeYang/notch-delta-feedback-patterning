"""Root-locus continuation vs. lam_N, copied (logic unchanged) from the
source repository's ``preprint_final_analysis/scripts/root_locus_and_gain.py``.

Tracks all four poles of J_mu for a given mode mu across a lam_N sweep. In
the final analysis: the mu=-3 (dominant, three-sublattice) leading pole
crosses Re(s)=0 at lam_N ~= 1.305, while mu=6 (uniform) stays stable
(max Re(s) = -0.1) across the whole tested range -- see
docs/CONTROL_THEORY_DERIVATION.md section 6.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from control.state_space import G_exact
from model.fixed_point import find_fixed_points
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams


def poles_at(mu: float, lam_N: float, p0: NotchDeltaParams) -> np.ndarray:
    """The four poles (sorted, most unstable first) of J_mu at this lam_N."""
    p = dataclasses.replace(p0, lam_N=lam_N)
    fixed_points = find_fixed_points(p)
    assert len(fixed_points) == 1, f"expected 1 fixed point at lam_N={lam_N}, found {len(fixed_points)}"
    J_mu = jacobian_block(mu, fixed_points[0], p)
    eigs = np.linalg.eigvals(J_mu)
    return eigs[np.argsort(-eigs.real)]


def root_locus(mu: float, lam_N_values, p0: NotchDeltaParams) -> list[dict]:
    rows = []
    for lam_N in lam_N_values:
        p = dataclasses.replace(p0, lam_N=lam_N)
        fixed_points = find_fixed_points(p)
        if len(fixed_points) != 1:
            continue
        fp = fixed_points[0]
        J_mu = jacobian_block(mu, fp, p)
        eigs = np.linalg.eigvals(J_mu)
        eigs = eigs[np.argsort(-eigs.real)]
        g_exact = G_exact(J_mu)
        for i, e in enumerate(eigs):
            rows.append({"lam_N": lam_N, "mu": mu, "pole_index": i,
                        "re": float(e.real), "im": float(e.imag), "G_exact_DC": g_exact})
    return rows


def find_critical_lam_N(mu: float, lam_N_values, p0: NotchDeltaParams) -> float | None:
    """Linear-interpolated lam_N where the leading pole of mode mu crosses Re(s)=0."""
    prev_re, prev_lam = None, None
    for lam_N in sorted(lam_N_values):
        re = float(poles_at(mu, lam_N, p0)[0].real)
        if prev_re is not None and prev_re < 0 <= re:
            frac = -prev_re / (re - prev_re)
            return prev_lam + frac * (lam_N - prev_lam)
        prev_re, prev_lam = re, lam_N
    return None
