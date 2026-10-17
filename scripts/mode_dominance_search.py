"""Search for parameter sets in which a spatial mode other than mu_min = -3
dominates or destabilises first.

For stationary (zero-eigenvalue) onset, docs/THEORY.md (Proposition 3) shows
that no such parameter set exists. Non-stationary instabilities are not
covered by that argument, so this script samples new random parameter sets
and runs the full mode scan:

  * 20 000 points from the v1 sampling box (rates log-uniform, Hill
    exponents uniform in [1, 6]);
  * 20 000 points from a box widened 10-fold in each direction for every rate.

For each point it records whether a nonuniform mode mu != -3 is the most
unstable mode, or is unstable while mu = -3 is stable. Fixed seed.
Writes figures/data/mode_dominance_search.csv (one summary row per box).
Run time is about 5 minutes.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from tqdm import tqdm

from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams
from stability.mode_scan import UNIFORM_MU, UNIQUE_MU
from stability.sensitivity import PARAM_NAMES

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"
LO = dict(lam_D=0.5, lam_N=0.1, lam_ASC=0.1, Kd=0.03, n_d=1, n_ASC=1, d_D=0.01, d_N=0.01, d_A=0.1,
          d_ASC=0.01, f_D=0.005, f_N=0.06, k_cis=0.035)
HI = dict(lam_D=50, lam_N=10, lam_ASC=10, Kd=3, n_d=6, n_ASC=6, d_D=1, d_N=1, d_A=10, d_ASC=1,
          f_D=0.5, f_N=6, k_cis=3.5)
N_PER_BOX = 20_000
SEED = 20261016


def draw(rng: np.random.Generator, widen: float) -> NotchDeltaParams:
    vals = {}
    for name in PARAM_NAMES:
        if name in ("n_d", "n_ASC"):
            vals[name] = rng.uniform(LO[name], HI[name])
        else:
            vals[name] = float(np.exp(rng.uniform(np.log(LO[name] / widen), np.log(HI[name] * widen))))
    return NotchDeltaParams(**vals)


def main() -> None:
    # One continuous PCG64 stream: base box first, widened box second.
    rng = np.random.Generator(np.random.PCG64(SEED))
    nonuniform = UNIQUE_MU < UNIFORM_MU - 1e-9
    i_K = int(np.argmin(UNIQUE_MU))
    rows, samples = [], []
    for widen in (1.0, 10.0):
        unstable = other_dominant = other_only = oscillatory_K = 0
        for sample_id in tqdm(range(N_PER_BOX), desc=f"box x{widen:g}"):
            p = draw(rng, widen)
            fp = homogeneous_fixed_point(p)
            ev = np.linalg.eigvals(np.stack([jacobian_block(mu, fp, p) for mu in UNIQUE_MU]))
            lead = ev.real.max(axis=1)
            k = int(np.argmax(np.where(nonuniform, lead, -np.inf)))
            samples.append({"box_widening": widen, "sample_id": sample_id,
                            **{name: getattr(p, name) for name in PARAM_NAMES},
                            "lambda_pattern": float(lead[k]), "best_mu": float(UNIQUE_MU[k]),
                            "lambda_K": float(lead[i_K]), "uniform_re": float(lead[-1])})
            if lead[k] > 0:
                unstable += 1
                other_dominant += int(k != i_K and lead[k] > lead[i_K] + 1e-12)
                other_only += int(lead[i_K] <= 0)
                lead_K = ev[i_K][np.argmax(ev[i_K].real)]
                oscillatory_K += int(abs(lead_K.imag) > 1e-9)
        rows.append({"box_widening": widen, "n_samples": N_PER_BOX, "n_nonuniform_unstable": unstable,
                     "n_other_mode_dominant": other_dominant, "n_other_mode_unstable_K_stable": other_only,
                     "n_K_leading_oscillatory": oscillatory_K})
        print(rows[-1])
    sample_frame = pd.DataFrame(samples)
    snapshot = DATA_DIR / "mode_search_samples.csv.gz"
    if snapshot.exists():
        frozen = pd.read_csv(snapshot, float_precision="round_trip")
        cols = ["box_widening", "sample_id", *PARAM_NAMES]
        if frozen.shape[0] != sample_frame.shape[0] or not np.allclose(
                frozen[cols].values, sample_frame[cols].values, rtol=5e-14, atol=0):
            raise RuntimeError("mode-search inputs differ from the frozen parameter snapshot")
    pd.DataFrame(rows).to_csv(DATA_DIR / "mode_dominance_search.csv", index=False)
    # Freeze the actual inputs, rather than relying only on seed + summary.
    sample_frame.to_csv(snapshot, index=False, compression={"method": "gzip", "mtime": 0})


if __name__ == "__main__":
    main()
