"""Direct root-finding verification of the mu=-3 critical lam_N.

The value used elsewhere in this repository (docs/MANUSCRIPT_NUMBERS.md,
figures/data/exact_loop_gain.csv) comes from a 121-point geomspace lam_N
grid, so its trailing digits are a grid artifact. This script instead
brackets the same crossing directly with scipy.optimize.brentq and checks
that the pole-crossing (Re(s)=0) and unity-gain (G_exact(-3)=1) conditions
agree to machine precision, as the exact Schur-complement identity
(control/state_space.py, control/schur_verification.py) requires.

Reuses control.root_locus / control.state_space / model.* /
stability.mode_scan unchanged; adds no new numerical logic beyond the
brentq calls themselves.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.optimize import brentq
from tqdm import tqdm

from control.root_locus import poles_at
from control.state_space import G_exact
from model.fixed_point import find_fixed_points
from model.jacobian import jacobian_block
from model.params import P0
from stability.mode_scan import DOMINANT_PATTERN_MU, UNIFORM_MU

BRACKET = (1.2, 1.4)
XTOL = 1e-13


def _leading_pole_re(lam_N: float, mu: float, pbar: tqdm) -> float:
    pbar.update(1)
    return float(poles_at(mu, lam_N, P0)[0].real)


def _g_exact_minus_one(lam_N: float, mu: float, pbar: tqdm) -> float:
    pbar.update(1)
    p = dataclasses.replace(P0, lam_N=lam_N)
    fp = find_fixed_points(p)[0]
    return G_exact(jacobian_block(mu, fp, p)) - 1.0


def main():
    mu = DOMINANT_PATTERN_MU
    with tqdm(desc="pole crossing (Re(s)=0)", unit="eval") as pbar:
        lam_pole = brentq(_leading_pole_re, *BRACKET, args=(mu, pbar), xtol=XTOL)
    with tqdm(desc="unity-gain (G_exact-1=0)", unit="eval") as pbar:
        lam_gain = brentq(_g_exact_minus_one, *BRACKET, args=(mu, pbar), xtol=XTOL)

    p_pole = dataclasses.replace(P0, lam_N=lam_pole)
    fp_pole = find_fixed_points(p_pole)[0]
    g_at_pole = G_exact(jacobian_block(mu, fp_pole, p_pole))
    g_uniform_at_pole = G_exact(jacobian_block(UNIFORM_MU, fp_pole, p_pole))

    p_gain = dataclasses.replace(P0, lam_N=lam_gain)
    fp_gain = find_fixed_points(p_gain)[0]
    pole_at_gain = float(np.linalg.eigvals(jacobian_block(mu, fp_gain, p_gain)).real.max())

    print(f"pole crossing lam_N   = {lam_pole:.9f}")
    print(f"unity-gain    lam_N   = {lam_gain:.9f}  (difference {abs(lam_pole - lam_gain):.1e})")
    print(f"G_exact(-3) at crossing     = {g_at_pole:.12f}")
    print(f"leading pole at unity-gain  = {pole_at_gain:.1e}")
    print(f"G_exact(+6) at crossing     = {g_uniform_at_pole:.4f}")


if __name__ == "__main__":
    main()
