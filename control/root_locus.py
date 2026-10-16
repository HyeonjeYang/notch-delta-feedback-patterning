"""Pole continuation in lam_N and direct location of the critical point.

All four poles of J_mu are tracked for a given mode mu across a lam_N
sweep. The critical lam_N is located by Brent root-finding on the
leading pole's real part. It is located a second time, independently, by
Brent root-finding on G(mu) - 1. Grid interpolation is not used.
"""

from __future__ import annotations

import dataclasses

import numpy as np
from scipy.optimize import brentq

from control.state_space import G_exact
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams


def block_at(mu: float, lam_N: float, p0: NotchDeltaParams) -> np.ndarray:
    p = dataclasses.replace(p0, lam_N=lam_N)
    return jacobian_block(mu, homogeneous_fixed_point(p), p)


def poles_at(mu: float, lam_N: float, p0: NotchDeltaParams) -> np.ndarray:
    """The four poles (most unstable first) of J_mu at this lam_N."""
    eigs = np.linalg.eigvals(block_at(mu, lam_N, p0))
    return eigs[np.argsort(-eigs.real)]


def root_locus(mu: float, lam_N_values, p0: NotchDeltaParams) -> list[dict]:
    rows = []
    for lam_N in lam_N_values:
        J_mu = block_at(mu, lam_N, p0)
        eigs = np.linalg.eigvals(J_mu)
        eigs = eigs[np.argsort(-eigs.real)]
        g = G_exact(J_mu)
        for i, e in enumerate(eigs):
            rows.append({"lam_N": lam_N, "mu": mu, "pole_index": i,
                         "re": float(e.real), "im": float(e.imag), "G_mu": g})
    return rows


def critical_lam_N(mu: float, p0: NotchDeltaParams, bracket=(1.2, 1.4), xtol: float = 1e-13) -> dict:
    """lam_N where the leading pole crosses Re(s)=0 and where G(mu)=1,
    each from its own brentq."""
    lam_pole = brentq(lambda l: float(poles_at(mu, l, p0)[0].real), *bracket, xtol=xtol)
    lam_gain = brentq(lambda l: G_exact(block_at(mu, l, p0)) - 1.0, *bracket, xtol=xtol)
    return {"lam_N_pole": lam_pole, "lam_N_gain": lam_gain,
            "G_at_pole": G_exact(block_at(mu, lam_pole, p0)),
            "pole_at_gain": float(poles_at(mu, lam_gain, p0)[0].real)}
