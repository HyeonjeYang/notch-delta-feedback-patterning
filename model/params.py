"""D-N-A-S (Delta-Notch-activated Notch/NICD-proneural/ASC) kinetic
parameters, copied unmodified from the source repository's
``src/model1_notch_delta_stochastic/core.py::NotchDeltaParams`` and
``LITERATURE_PARAMS``. No value here was changed.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass


@dataclass(frozen=True)
class NotchDeltaParams:
    """The 13 Notch-Delta kinetic parameters."""

    lam_D:   float = 1.3   # lam_D  : max Delta production (Hill-modulated by S/ASC)
    lam_N:   float = 1.0   # lam_N  : constitutive Notch production
    lam_ASC: float = 1.0   # lam_S  : max ASC/S production (Hill-inhibited by A/NICD)
    Kd:      float = 0.3   # Kd     : shared EC50 for both Hill functions
    n_d:     float = 2.0   # n_d    : Hill coefficient, S -> D production
    n_ASC:   float = 2.0   # n_ASC  : Hill coefficient, A -> S inhibition
    d_D:     float = 0.1   # d_D    : Delta degradation rate
    d_N:     float = 0.1   # d_N    : Notch degradation rate
    d_A:     float = 0.1   # d_A    : NICD/A degradation rate
    d_ASC:   float = 0.1   # d_S    : ASC/S degradation rate
    f_D:     float = 0.6   # f_D    : trans Delta-Notch depletion rate
    f_N:     float = 0.6   # f_N    : trans Notch signalling / NICD production rate
    k_cis:   float = 0.5   # k_cis  : cis mutual inactivation rate


# Literature-inspired parameters producing a genuine Turing-like lateral
# inhibition regime (see docs/GEOMETRY_DERIVATION.md, docs/CONTROL_THEORY_DERIVATION.md).
LITERATURE_PARAMS = NotchDeltaParams(
    lam_D=5.0, lam_N=1.0, lam_ASC=1.0, Kd=0.3, n_d=3.0, n_ASC=3.0,
    d_D=0.1, d_N=0.1, d_A=1.0, d_ASC=0.1, f_D=0.05, f_N=0.6, k_cis=0.25,
)

# The validated continuous-time baseline p0 used throughout the final
# analysis: LITERATURE_PARAMS with k_cis overridden to 0.35 (the primary
# regime confirmed spatially stable before lam_N was identified as the
# relevant bifurcation parameter).
P0 = dataclasses.replace(LITERATURE_PARAMS, k_cis=0.35)

# The three validated points on the lam_N continuation axis (all other 12
# parameters fixed at P0), from the frozen patternability analysis.
LAM_N_BASELINE = 1.0    # spatially stable
LAM_N_P2 = 2.0          # robust validated pattern-forming point
LAM_N_P3 = 5.623        # deeper validated pattern-forming point
