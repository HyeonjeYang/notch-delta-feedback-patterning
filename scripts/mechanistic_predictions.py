"""Recompute compensated-turnover predictions stated in the Discussion."""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from scipy.optimize import brentq

from control.state_space import G_exact
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import P0
from model.reduced import reduce_S, scale_proneural_timescale


def main():
    base = dataclasses.replace(P0,lam_N=2)
    rows = []
    for alpha in (.25,1.,4.):
        p = scale_proneural_timescale(base,alpha)
        fp = homogeneous_fixed_point(p)
        J = jacobian_block(-3,fp,p)

        def threshold(lam_N):
            q = dataclasses.replace(p,lam_N=lam_N)
            return G_exact(jacobian_block(-3,homogeneous_fixed_point(q),q))-1

        rows.append(dict(alpha=alpha,lam_S=p.lam_ASC,d_S=p.d_ASC,
            S_q_prefactor=p.lam_ASC/p.d_ASC,S_star=fp[3],G_minus3=G_exact(J),
            stationary_lam_N=brentq(threshold,1.2,1.4,xtol=1e-13),
            growth_full=float(np.linalg.eigvals(J).real.max()),
            growth_DNA=float(np.linalg.eigvals(reduce_S(J)).real.max()),
            clamped_A_relaxation_time=1/p.d_ASC,clamped_A_half_power_frequency=p.d_ASC))
    pd.DataFrame(rows).to_csv(ROOT/"figures/data/mechanistic_predictions.csv",index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
