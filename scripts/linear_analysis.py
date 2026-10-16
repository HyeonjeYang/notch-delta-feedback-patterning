"""Linear analysis along the lam_N continuation (all other parameters at p0).

Writes to figures/data/:
  pole_continuation.csv    four poles of J_mu for mu = -3 and mu = 6 on 121
                           log-spaced lam_N values in [0.5, 6.5], with G(mu)
  critical_point.csv       critical lam_N from two independent brentq
                           calls, one on the pole crossing and one on G(-3) = 1
  loop_gain_spectrum.csv   G(mu) and the stage gains over every distinct graph
                           eigenvalue at lam_N = 1, critical, 2, 5.623
  schur_verification.csv   direct vs Schur-complement characteristic
                           polynomial (9 (lam_N, mu) combinations)
  graph_io_check.csv       det J_mu vs det(J_loc) det(I - mu K) on the same
                           combinations (control/graph_io.py)
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from control.graph_io import coupling_decomposition, static_io_gain
from control.root_locus import block_at, critical_lam_N, root_locus
from control.schur_verification import verify as schur_verify
from control.state_space import stage_gains
from model.fixed_point import homogeneous_fixed_point
from model.params import LAM_N_BASELINE, LAM_N_P2, LAM_N_P3, P0
from stability.mode_scan import DOMINANT_PATTERN_MU, UNIFORM_MU, UNIQUE_MU

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"


def main() -> None:
    lam_grid = np.geomspace(0.5, 6.5, 121)
    rows = root_locus(DOMINANT_PATTERN_MU, lam_grid, P0) + root_locus(UNIFORM_MU, lam_grid, P0)
    pd.DataFrame(rows).to_csv(DATA_DIR / "pole_continuation.csv", index=False)

    crit = critical_lam_N(DOMINANT_PATTERN_MU, P0)
    lam_c = crit["lam_N_pole"]
    crit["G_uniform_at_pole"] = float(stage_gains(block_at(UNIFORM_MU, lam_c, P0))["G"])
    crit["max_re_uniform_over_grid"] = float(max(r["re"] for r in rows if r["mu"] == UNIFORM_MU))
    pd.DataFrame([crit]).to_csv(DATA_DIR / "critical_point.csv", index=False, float_format="%.15g")
    print(f"critical lam_N: pole {crit['lam_N_pole']:.9f}, gain {crit['lam_N_gain']:.9f}, "
          f"|diff| {abs(crit['lam_N_pole'] - crit['lam_N_gain']):.1e}")
    print(f"G(-3) at pole crossing {crit['G_at_pole']:.12f}; G(+6) there {crit['G_uniform_at_pole']:.4f}; "
          f"max Re pole mu=6 over grid {crit['max_re_uniform_over_grid']:.6f}")

    spec = []
    for label, lam_N in [("baseline", LAM_N_BASELINE), ("critical", lam_c), ("P2", LAM_N_P2), ("P3", LAM_N_P3)]:
        for mu in UNIQUE_MU:
            spec.append({"label": label, "lam_N": lam_N, "mu": mu, **stage_gains(block_at(mu, lam_N, P0))})
    spec = pd.DataFrame(spec)
    spec.to_csv(DATA_DIR / "loop_gain_spectrum.csv", index=False)
    am = spec.loc[spec.groupby("label")["G"].idxmax(), ["label", "mu", "G"]]
    print("argmax_mu G(mu):", am.to_dict("records"))

    schur_rows, io_rows = [], []
    for lam_N in (LAM_N_BASELINE, LAM_N_P2, LAM_N_P3):
        p = dataclasses.replace(P0, lam_N=lam_N)
        fp = homogeneous_fixed_point(p)
        J_loc, _, _ = coupling_decomposition(fp, p)
        K = static_io_gain(fp, p)
        for mu in (UNIFORM_MU, DOMINANT_PATTERN_MU, 1.5):
            J = block_at(mu, lam_N, P0)
            res = schur_verify(J)
            schur_rows.append({"lam_N": lam_N, "mu": mu, "max_coeff_abs_diff": res["max_coeff_diff"],
                               "max_root_abs_mismatch": res["max_root_mismatch"]})
            lhs = np.linalg.det(J)
            rhs = np.linalg.det(J_loc) * np.linalg.det(np.eye(2) - mu * K)
            io_rows.append({"lam_N": lam_N, "mu": mu, "det_J": lhs, "det_Jloc_det_I_minus_muK": rhs,
                            "rel_diff": abs(lhs - rhs) / abs(lhs)})
    schur = pd.DataFrame(schur_rows); schur.to_csv(DATA_DIR / "schur_verification.csv", index=False)
    io = pd.DataFrame(io_rows); io.to_csv(DATA_DIR / "graph_io_check.csv", index=False)
    print(f"Schur: max coeff diff {schur.max_coeff_abs_diff.max():.2e}, max root mismatch {schur.max_root_abs_mismatch.max():.2e}")
    print(f"graph I/O identity: max relative difference {io.rel_diff.max():.2e}")


if __name__ == "__main__":
    main()
