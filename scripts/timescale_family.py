"""Proneural-timescale family: (lam_S, d_S) -> alpha * (lam_S, d_S).

Every member of the family has the same homogeneous equilibrium and the same
static gains, so it also has the same G(mu) for every mu and the same
stationary boundary. Only the relaxation rate d_S of the proneural stage
changes. alpha = 1 is the explicit model. alpha -> infinity is the
matched DNA reduction (S slaved to A).

Writes figures/data/timescale_family.csv with, for each alpha:
  * the leading mu = -3 eigenvalue at P2 and P3 (full model);
  * the DNA- and DN-reduction growth rates, which do not depend on alpha;
  * the slow-proneural estimate d_S (G(-3) - 1);
  * G(-3), and the critical lam_N recomputed by brentq for a subset of
    alpha values, as a check that it does not depend on alpha.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from control.root_locus import critical_lam_N
from control.state_space import G_exact
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import LAM_N_P2, LAM_N_P3, P0
from model.reduced import reduce_AS, reduce_S, scale_proneural_timescale
from stability.mode_scan import DOMINANT_PATTERN_MU, lambda_pattern, leading_eigenvalue

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"
ALPHAS = np.geomspace(1e-2, 1e4, 61)
CRIT_CHECK_ALPHAS = (1e-2, 1e-1, 1.0, 1e1, 1e2, 1e4)


def main() -> None:
    rows = []
    for alpha in ALPHAS:
        row = {"alpha": alpha}
        for label, lam_N in (("P2", LAM_N_P2), ("P3", LAM_N_P3)):
            p = scale_proneural_timescale(dataclasses.replace(P0, lam_N=lam_N), alpha)
            fp = homogeneous_fixed_point(p)
            J = jacobian_block(DOMINANT_PATTERN_MU, fp, p)
            lead = leading_eigenvalue(J)
            scan = lambda_pattern(fp, p)
            G = G_exact(J)
            row.update({f"{label}_lambda": lead.real, f"{label}_lambda_im": abs(lead.imag),
                        f"{label}_best_mu": scan["nonuniform_best_mu"],
                        f"{label}_uniform_re": scan["uniform_re_lambda"],
                        f"{label}_G_minus3": G, f"{label}_slow_S_estimate": p.d_ASC * (G - 1.0),
                        f"{label}_lambda_DNA": leading_eigenvalue(reduce_S(J)).real,
                        f"{label}_lambda_DN": leading_eigenvalue(reduce_AS(J)).real,
                        f"{label}_D_star": fp[0]})
        if any(np.isclose(alpha, a) for a in CRIT_CHECK_ALPHAS):
            row["critical_lam_N"] = critical_lam_N(DOMINANT_PATTERN_MU, scale_proneural_timescale(P0, alpha),
                                                   bracket=(1.2, 1.4))["lam_N_pole"]
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(DATA_DIR / "timescale_family.csv", index=False)

    at1 = df.iloc[int(np.argmin(np.abs(np.log(df.alpha))))]
    print(f"alpha=1: P2 lambda {at1.P2_lambda:.6f}, DNA {at1.P2_lambda_DNA:.6f} "
          f"(x{at1.P2_lambda_DNA / at1.P2_lambda:.1f}), DN {at1.P2_lambda_DN:.4f} (x{at1.P2_lambda_DN / at1.P2_lambda:.0f}); "
          f"slow-S estimate {at1.P2_slow_S_estimate:.4f}")
    print(f"           P3 lambda {at1.P3_lambda:.6f}, DNA {at1.P3_lambda_DNA:.6f} (x{at1.P3_lambda_DNA / at1.P3_lambda:.1f})")
    print(f"P2 lambda range over alpha: {df.P2_lambda.min():.2e} .. {df.P2_lambda.max():.4f}")
    print(f"G(-3) spread over alpha (P2): {df.P2_G_minus3.max() - df.P2_G_minus3.min():.1e}; "
          f"D* spread {df.P2_D_star.max() - df.P2_D_star.min():.1e}")
    crit = df.critical_lam_N.dropna()
    print(f"critical lam_N over alpha in {CRIT_CHECK_ALPHAS}: {crit.min():.12f} .. {crit.max():.12f}")
    print(f"dominant mode always -3: {bool(np.allclose(df[['P2_best_mu', 'P3_best_mu']], -3))}; "
          f"oscillatory leading eigenvalue anywhere: {bool((df[['P2_lambda_im', 'P3_lambda_im']] > 1e-9).any().any())}; "
          f"uniform mode unstable anywhere: {bool((df[['P2_uniform_re', 'P3_uniform_re']] > 0).any().any())}")


if __name__ == "__main__":
    main()
