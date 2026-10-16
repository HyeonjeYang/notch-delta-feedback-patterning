"""Reproduces figures/data/geff_comparison.csv (and the 99.79%-vs-93.88%
classification numbers in docs/CONTROL_THEORY.md) directly from
figures/data/global_parameter_samples.csv -- the frozen 4000-point global
parameter sample, included in full so this comparison is reproducible
in-repo without any new search. Logic copied unmodified from the source
repository's ``preprint_final_analysis/scripts/explain_geff.py``; no
parameter search is (re)run here, only a recomputation of G_exact/G_eff
from the already-stored fixed points.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from control.loop_gain import G_eff as G_eff_fn
from control.state_space import G_exact as G_exact_fn
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"
MU_PROBE = -3.0

PARAM_NAMES = ["lam_D", "lam_N", "lam_ASC", "Kd", "n_d", "n_ASC",
              "d_D", "d_N", "d_A", "d_ASC", "f_D", "f_N", "k_cis"]


def main() -> None:
    df = pd.read_csv(DATA_DIR / "global_parameter_samples.csv")
    valid = df[df.classification != "no_positive_fixed_point"].copy()

    rows = []
    for _, r in valid.iterrows():
        p = NotchDeltaParams(*[r[name] for name in PARAM_NAMES])
        fp = np.array([r["D_star"], r["N_star"], r["A_star"], r["ASC_star"]])
        J_mu = jacobian_block(MU_PROBE, fp, p)
        rows.append({
            "sample_id": r["sample_id"], "lambda_pattern": r["nonuniform_max_re_lambda"],
            "G_exact_minus3": G_exact_fn(J_mu), "G_eff": G_eff_fn(J_mu),
            "unstable": bool(r["nonuniform_max_re_lambda"] > 0),
        })
    gdf = pd.DataFrame(rows)

    pred_exact = gdf["G_exact_minus3"].values > 1.0
    acc_exact = (pred_exact == gdf["unstable"].values).mean()

    thresholds = np.linspace(gdf["G_eff"].min(), gdf["G_eff"].max(), 300)
    y = gdf["unstable"].values
    best_acc, best_t = 0.0, None
    for t in thresholds:
        acc = ((gdf["G_eff"].values > t) == y).mean()
        if acc > best_acc:
            best_acc, best_t = acc, t

    print(f"n valid points: {len(gdf)}")
    print(f"corr(G_eff, G_exact(-3)) = {gdf[['G_eff','G_exact_minus3']].corr().iloc[0,1]:.4f}")
    print(f"G_exact(-3)>1 classification accuracy: {acc_exact:.4f}")
    print(f"G_eff best-threshold accuracy: {best_acc:.4f} at G_eff={best_t:.4f}")

    gdf.to_csv(DATA_DIR / "geff_comparison.csv", index=False)
    print(f"wrote {DATA_DIR / 'geff_comparison.csv'}")


if __name__ == "__main__":
    main()
