"""Local log-sensitivities and one-parameter continuations.

Two margins are differentiated:

* lambda_pattern, the leading spatial growth (or decay) rate. It depends on
  static gains and on the kinetic timescales.
* ln G(-3), the static loop gain of the three-sublattice mode. Its sign
  relative to 0 sets the stationary boundary. It splits exactly into
  three stage contributions (control.state_space.stage_gains):

      d ln G / d ln p = d ln kappa_D + d ln T_-3 + d ln kappa_S

  Each term includes the parameter's indirect effect through the shift of
  the homogeneous equilibrium.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from control.state_space import stage_gains
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams
from stability.mode_scan import DOMINANT_PATTERN_MU, lambda_pattern

PARAM_NAMES = ["lam_D", "lam_N", "lam_ASC", "Kd", "n_d", "n_ASC",
               "d_D", "d_N", "d_A", "d_ASC", "f_D", "f_N", "k_cis"]


def _lambda_pattern_at(p: NotchDeltaParams) -> float:
    return lambda_pattern(homogeneous_fixed_point(p), p)["nonuniform_max_re_lambda"]


def _log_stage_gains(p: NotchDeltaParams, mu: float = DOMINANT_PATTERN_MU) -> np.ndarray:
    g = stage_gains(jacobian_block(mu, homogeneous_fixed_point(p), p))
    return np.log([g["kappa_D"], g["T_mu"], g["kappa_S"]])


def _perturbed(p0: NotchDeltaParams, name: str, h: float):
    base = getattr(p0, name)
    return (dataclasses.replace(p0, **{name: base * np.exp(h)}),
            dataclasses.replace(p0, **{name: base * np.exp(-h)}))


def local_sensitivity(p0: NotchDeltaParams, h: float = 0.01) -> dict[str, float]:
    """S_i = d lambda_pattern / d ln p_i (central difference)."""
    out = {}
    for name in PARAM_NAMES:
        p_plus, p_minus = _perturbed(p0, name, h)
        out[name] = (_lambda_pattern_at(p_plus) - _lambda_pattern_at(p_minus)) / (2 * h)
    return out


def stage_gain_sensitivity(p0: NotchDeltaParams, h: float = 1e-4) -> list[dict]:
    """d ln G(-3) / d ln p_i, split into (kappa_D, T_-3, kappa_S) parts."""
    rows = []
    for name in PARAM_NAMES:
        p_plus, p_minus = _perturbed(p0, name, h)
        d = (_log_stage_gains(p_plus) - _log_stage_gains(p_minus)) / (2 * h)
        rows.append({"parameter": name, "dlnG_total": float(d.sum()),
                     "via_kappa_D": float(d[0]), "via_T": float(d[1]), "via_kappa_S": float(d[2])})
    return rows


def one_dimensional_scan(p0: NotchDeltaParams, name: str, multipliers) -> list[dict]:
    """Continuation of a single parameter, multiplicatively about p0."""
    rows = []
    base = getattr(p0, name)
    for mult in multipliers:
        p = dataclasses.replace(p0, **{name: base * mult})
        fp = homogeneous_fixed_point(p)
        g = stage_gains(jacobian_block(DOMINANT_PATTERN_MU, fp, p))["G"]
        rows.append({"parameter": name, "multiplier": mult, "value": base * mult,
                     "D_star": fp[0], "N_star": fp[1], "A_star": fp[2], "S_star": fp[3],
                     "G_minus3": g, **lambda_pattern(fp, p)})
    return rows
