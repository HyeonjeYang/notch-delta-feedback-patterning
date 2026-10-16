"""Which kinetic parameters move the spatial stability margin (all about p0).

Writes to figures/data/:
  local_sensitivity.csv          d lambda_pattern / d ln p_i at p0, h = 0.01
  stage_gain_sensitivity.csv     d ln G(-3) / d ln p_i at p0, split into the
                                 kappa_D, T_-3 and kappa_S stage contributions
  one_dimensional_scans.csv      every parameter x [0.1, 10] (41 log-spaced
                                 multipliers); Hill exponents over [1, 6]
  two_dimensional_phase_maps.csv 25 x 25 grids for (lam_N, d_A),
                                 (n_ASC, Kd) and (lam_N, n_ASC). Same axes as v1,
                                 with no missing equilibria now.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from control.state_space import G_exact
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import P0
from stability.mode_scan import DOMINANT_PATTERN_MU, lambda_pattern
from stability.sensitivity import PARAM_NAMES, local_sensitivity, one_dimensional_scan, stage_gain_sensitivity

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"
HILL = ("n_d", "n_ASC")
GRIDS = {
    "lam_N": np.geomspace(0.2, 10.0, 25), "d_A": np.geomspace(0.05, 2.0, 25),
    "n_ASC": np.linspace(1.0, 6.0, 25), "Kd": np.geomspace(0.03, 1.0, 25),
}
PAIRS = [("lam_N", "d_A"), ("n_ASC", "Kd"), ("lam_N", "n_ASC")]


def main() -> None:
    sens = local_sensitivity(P0)
    sdf = pd.DataFrame({"parameter": list(sens), "baseline_value": [getattr(P0, k) for k in sens],
                        "S_lambda_pattern": list(sens.values())})
    sdf["rank"] = sdf.S_lambda_pattern.abs().rank(ascending=False).astype(int)
    sdf.sort_values("rank").to_csv(DATA_DIR / "local_sensitivity.csv", index=False)

    gdf = pd.DataFrame(stage_gain_sensitivity(P0))
    gdf["rank"] = gdf.dlnG_total.abs().rank(ascending=False).astype(int)
    gdf.sort_values("rank").to_csv(DATA_DIR / "stage_gain_sensitivity.csv", index=False)
    print("d lambda_pattern / d ln p (top 6):", sdf.sort_values("rank").head(6)[["parameter", "S_lambda_pattern"]].round(4).to_dict("records"))
    print("d ln G(-3) / d ln p (top 6):", gdf.sort_values("rank").head(6)[["parameter", "dlnG_total", "via_kappa_D", "via_T", "via_kappa_S"]].round(3).to_dict("records"))

    rows = []
    for name in PARAM_NAMES:
        base = getattr(P0, name)
        mults = np.linspace(1.0, 6.0, 21) / base if name in HILL else np.geomspace(0.1, 10.0, 41)
        rows += one_dimensional_scan(P0, name, mults)
    one = pd.DataFrame(rows)
    one.to_csv(DATA_DIR / "one_dimensional_scans.csv", index=False)
    crossing = one.groupby("parameter").apply(lambda g: bool((np.sign(g.nonuniform_max_re_lambda).diff().fillna(0) != 0).any()),
                                              include_groups=False)
    print("1-D scans that cross lambda_pattern = 0:", sorted(crossing[crossing].index))

    rows = []
    for px, py in PAIRS:
        for x in GRIDS[px]:
            for y in GRIDS[py]:
                p = dataclasses.replace(P0, **{px: x, py: y})
                fp = homogeneous_fixed_point(p)
                scan = lambda_pattern(fp, p)
                rows.append({"pair": f"{px}_x_{py}", "param_x": px, "param_y": py, "x": x, "y": y,
                             "G_minus3": G_exact(jacobian_block(DOMINANT_PATTERN_MU, fp, p)), **scan})
    two = pd.DataFrame(rows)
    two.to_csv(DATA_DIR / "two_dimensional_phase_maps.csv", index=False)
    unstable = two.nonuniform_max_re_lambda > 0
    print(f"2-D maps: {len(two)} points, unstable {unstable.sum()}, "
          f"G>1 agrees with instability at {((two.G_minus3 > 1) == unstable).mean():.4%}, "
          f"oscillatory leading modes {(two.nonuniform_oscillatory & unstable).sum()}")


if __name__ == "__main__":
    main()
