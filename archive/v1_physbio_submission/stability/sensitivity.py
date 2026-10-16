"""Local log-sensitivity and one-dimensional parameter continuation of
lambda_pattern, condensed from the source repository's
``patternability_map/scripts/phase1_sensitivity.py`` and
``phase2_1d_scans.py``. Logic unchanged; the two scripts are combined here
since the reduced results (figures/data/*.csv) are what the paper actually
cites -- this module documents and can regenerate that data, it is not
required to reproduce Figures 1-5 (which read the already-included CSVs).
"""

from __future__ import annotations

import dataclasses

import numpy as np

from model.fixed_point import find_fixed_points
from model.params import NotchDeltaParams
from stability.mode_scan import lambda_pattern

PARAM_NAMES = ["lam_D", "lam_N", "lam_ASC", "Kd", "n_d", "n_ASC",
              "d_D", "d_N", "d_A", "d_ASC", "f_D", "f_N", "k_cis"]


def _lambda_pattern_at(p: NotchDeltaParams) -> tuple[float, float]:
    branches = find_fixed_points(p)
    if not branches:
        return float("nan"), float("nan")
    best = max((lambda_pattern(fp, p) for fp in branches),
              key=lambda d: d["nonuniform_max_re_lambda"])
    return best["nonuniform_max_re_lambda"], best["uniform_re_lambda"]


def local_sensitivity(p0: NotchDeltaParams, h: float = 0.01) -> dict[str, float]:
    """S_i = d(lambda_pattern) / d(log p_i), central difference in log-space."""
    out = {}
    for name in PARAM_NAMES:
        base = getattr(p0, name)
        p_plus = dataclasses.replace(p0, **{name: base * np.exp(h)})
        p_minus = dataclasses.replace(p0, **{name: base * np.exp(-h)})
        nonuni_plus, _ = _lambda_pattern_at(p_plus)
        nonuni_minus, _ = _lambda_pattern_at(p_minus)
        out[name] = (nonuni_plus - nonuni_minus) / (2 * h)
    return out


def one_dimensional_scan(p0: NotchDeltaParams, name: str, multipliers) -> list[dict]:
    """Continuation of a single parameter (multiplicative for rate
    constants; pass raw values via `multipliers` for Hill exponents)."""
    rows = []
    base = getattr(p0, name)
    for mult in multipliers:
        p = dataclasses.replace(p0, **{name: base * mult})
        branches = find_fixed_points(p)
        if not branches:
            rows.append({"multiplier": mult, "n_branches": 0})
            continue
        for i, fp in enumerate(branches):
            stab = lambda_pattern(fp, p)
            rows.append({"multiplier": mult, "branch": i, "n_branches": len(branches),
                        "D_star": fp[0], "N_star": fp[1], "A_star": fp[2], "S_star": fp[3], **stab})
    return rows
