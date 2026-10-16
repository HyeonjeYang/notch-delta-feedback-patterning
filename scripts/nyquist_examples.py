"""Data for the 'beyond the static gain' figure.

Writes to figures/data/:
  nyquist_examples.csv     L_-3(i omega), omega >= 0, for three points: p0 (stable),
                           P2 (stationary instability), and global-sample point
                           1730 (complex pair in the right half-plane while
                           G(-3) = 0.366)
  two_eigenvalue_cases.csv every global-sample point whose mu = -3 block has
                           exactly two right-half-plane eigenvalues, with those
                           eigenvalues and G(-3)

Requires figures/data/global_sample_results.csv (scripts/global_sample.py).
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from control.state_space import G_exact, nyquist_curve
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import LAM_N_P2, NotchDeltaParams, P0
from stability.mode_scan import DOMINANT_PATTERN_MU
from stability.sensitivity import PARAM_NAMES

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"
EXAMPLE_SAMPLE = 1730


def sample_params(params: pd.DataFrame, sample_id: int) -> NotchDeltaParams:
    r = params.set_index("sample_id").loc[sample_id]
    return NotchDeltaParams(*[float(r[n]) for n in PARAM_NAMES])


def main() -> None:
    params = pd.read_csv(DATA_DIR / "global_sample_parameters.csv", float_precision="round_trip")
    results = pd.read_csv(DATA_DIR / "global_sample_results.csv")
    examples = [("p0 (stable)", P0), ("P2 (stationary)", dataclasses.replace(P0, lam_N=LAM_N_P2)),
                (f"sample {EXAMPLE_SAMPLE} (complex pair)", sample_params(params, EXAMPLE_SAMPLE))]
    omega = np.concatenate([[0.0], np.geomspace(1e-4, 1e3, 3000)])
    rows = []
    for label, p in examples:
        J = jacobian_block(DOMINANT_PATTERN_MU, homogeneous_fixed_point(p), p)
        L = nyquist_curve(J, omega)
        ev = np.linalg.eigvals(J)
        rows.append(pd.DataFrame({"example": label, "omega": omega, "re_L": L.real, "im_L": L.imag,
                                  "G_minus3": G_exact(J), "n_rhp": int((ev.real > 0).sum())}))
    pd.concat(rows).to_csv(DATA_DIR / "nyquist_examples.csv", index=False)

    cases = []
    for sid in results.loc[results.n_rhp_minus3 == 2, "sample_id"]:
        p = sample_params(params, int(sid))
        J = jacobian_block(DOMINANT_PATTERN_MU, homogeneous_fixed_point(p), p)
        ev = np.linalg.eigvals(J)
        for e in ev[ev.real > 0]:
            cases.append({"sample_id": int(sid), "G_minus3": G_exact(J), "re": e.real, "im": e.imag,
                          "pair_type": "complex" if abs(e.imag) > 1e-12 else "real",
                          "v1_classification": results.set_index("sample_id").loc[sid, "v1_classification"]})
    cases = pd.DataFrame(cases)
    cases.to_csv(DATA_DIR / "two_eigenvalue_cases.csv", index=False)
    print(cases.groupby("sample_id").agg(G=("G_minus3", "first"), type=("pair_type", "first"),
                                         v1=("v1_classification", "first")).to_string())


if __name__ == "__main__":
    main()
