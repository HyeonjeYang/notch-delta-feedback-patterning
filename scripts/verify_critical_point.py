"""Locate the mu = -3 critical lam_N twice and check that the two locations agree.

Two separate brentq calls are made on [1.2, 1.4] with xtol = 1e-13: one on the
leading pole's real part (Re s = 0) and one on G(-3) - 1. They agree to
machine precision, as the static identity
det J_mu = det(A_mu) a_S [1 - G(mu)] requires (control/state_space.py).
Grid interpolation is not used.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from control.root_locus import block_at, critical_lam_N
from control.state_space import G_exact
from model.params import P0
from stability.mode_scan import DOMINANT_PATTERN_MU, UNIFORM_MU


def main() -> None:
    crit = critical_lam_N(DOMINANT_PATTERN_MU, P0)
    print(f"pole crossing lam_N   = {crit['lam_N_pole']:.9f}")
    print(f"unity-gain    lam_N   = {crit['lam_N_gain']:.9f}  "
          f"(difference {abs(crit['lam_N_pole'] - crit['lam_N_gain']):.1e})")
    print(f"G(-3) at crossing     = {crit['G_at_pole']:.12f}")
    print(f"leading pole at unity gain = {crit['pole_at_gain']:.1e}")
    print(f"G(+6) at crossing     = {G_exact(block_at(UNIFORM_MU, crit['lam_N_pole'], P0)):.4f}")


if __name__ == "__main__":
    main()
