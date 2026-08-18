"""Part A: structure factor of the already-completed stochastic sweep.

The original stochastic_check/run_stochastic_check.py saved scalar order
parameters (q_amp, std(D), P_mu=-3, Moran's I) per recorded time point, but
never persisted the raw final spatial fields. Since every run there was
fully deterministic given its (seed, sigma) pair (one master RNG per sigma
group, drawing the whole (n_replicates,24,24) noise array at every step),
this script REPLAYS the exact same 48 trajectories -- same code, same
seeds, same sigma values -- solely to recover the final D fields for
spectral analysis. This is not a new stochastic experiment: no sigma
value, seed, or parameter differs from the original sweep. Before trusting
the replay, every replica's final q_amp and final P_mu=-3 are checked
against the already-saved stochastic_check/per_replicate_summary.csv to
machine precision; the script aborts if they do not match.
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_STOCH = _HERE.parent
sys.path.insert(0, str(_STOCH))

import numpy as np
import pandas as pd

from run_stochastic_check import (  # noqa: E402
    COLS, DT_MAIN, LAM_N_P2, N_REPLICATES, P0, REL_AMP, ROWS, SIGMAS, T_FINAL,
    build_projection_basis, dominant_eigvec, find_fixed_points, mode_purity,
    order_parameters, perturbed_field, run_batch,
)
import dataclasses  # noqa: E402

OUT_DIR = _HERE
K0, N0 = 8, 16  # the (m,n) indices realizing mu=-3


def structure_factor(D: np.ndarray) -> np.ndarray:
    """S(m,n) = |FFT[D-mean(D)]|^2, raw (unnormalized)."""
    z = D - D.mean()
    Z = np.fft.fft2(z)
    return np.abs(Z) ** 2


def normalize_excluding_dc(S: np.ndarray) -> np.ndarray:
    S = S.copy()
    total = S.sum() - S[0, 0]  # exclude k=0 (identically ~0 for a mean-subtracted field)
    S[0, 0] = 0.0
    return S / total if total > 0 else S


def spectral_entropy(S_norm: np.ndarray) -> float:
    p = S_norm.ravel()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def peak_location(S_norm: np.ndarray) -> tuple[int, int, float]:
    S = S_norm.copy()
    S[0, 0] = -1  # exclude DC
    idx = np.unravel_index(np.argmax(S), S.shape)
    return int(idx[0]), int(idx[1]), float(S_norm[idx])


def k_point_power_fraction(S_norm: np.ndarray, k0: int = K0, n0: int = N0) -> float:
    return float(S_norm[k0, n0] + S_norm[n0, k0])


def main() -> None:
    p = dataclasses.replace(P0, lam_N=LAM_N_P2)
    fp = find_fixed_points(p)[0]
    D_star = float(fp[0])
    vec = dominant_eigvec(fp, p)
    D0, N0_, A0, S0 = perturbed_field(fp, vec, REL_AMP)
    Q = build_projection_basis()

    frozen_summary = pd.read_csv(_STOCH / "per_replicate_summary.csv")

    rows = []
    example_fields = {}
    seed_base = 1000
    for sigma_idx, sigma in enumerate(SIGMAS):
        result = run_batch(sigma, N_REPLICATES, DT_MAIN, T_FINAL, D0, N0_, A0, S0,
                           p, D_star, Q, seed=seed_base + sigma_idx)
        final_D = result["final_D"]  # (n_replicates, ROWS, COLS)

        # --- verification against the already-saved summary (replay integrity check) ---
        frozen_group = frozen_summary[frozen_summary.sigma == sigma].sort_values("replicate")
        replay_q_amp = np.array([order_parameters(final_D[i])["q_amp"] for i in range(N_REPLICATES)])
        replay_purity = np.array([mode_purity(final_D[i], Q) for i in range(N_REPLICATES)])
        max_q_amp_diff = float(np.max(np.abs(replay_q_amp - frozen_group["final_q_amp"].values)))
        max_purity_diff = float(np.max(np.abs(replay_purity - frozen_group["final_P_mu_minus3"].values)))
        print(f"sigma={sigma:g}: replay-vs-frozen max |delta q_amp|={max_q_amp_diff:.2e}, "
             f"max |delta P_mu-3|={max_purity_diff:.2e}")
        if max_q_amp_diff > 1e-8 or max_purity_diff > 1e-8:
            raise RuntimeError(f"Replay mismatch at sigma={sigma}: not a faithful reproduction "
                               f"of the frozen stochastic sweep -- STOPPING per task instruction.")

        for rep in range(N_REPLICATES):
            D = final_D[rep]
            S_raw = structure_factor(D)
            S_norm = normalize_excluding_dc(S_raw)
            k_frac = k_point_power_fraction(S_norm)
            purity = mode_purity(D, Q)
            k_diff = abs(k_frac - purity)

            m_peak, n_peak, peak_val = peak_location(S_norm)
            entropy = spectral_entropy(S_norm)

            rows.append({
                "sigma": sigma, "replicate": rep,
                "k_point_power_fraction": k_frac, "P_mu_minus3_recomputed": purity,
                "k_vs_purity_abs_diff": k_diff,
                "peak_m": m_peak, "peak_n": n_peak, "peak_value": peak_val,
                "peak_is_K_point": bool({(m_peak, n_peak)} <= {(K0, N0), (N0, K0)}),
                "spectral_entropy": entropy,
            })
            if rep == 0:
                example_fields[sigma] = S_norm

    df = pd.DataFrame(rows)
    max_verify_diff = df["k_vs_purity_abs_diff"].max()
    print(f"\nmax |k_point_power_fraction - P_mu_minus3| across all 48 replicates: {max_verify_diff:.3e}")
    if max_verify_diff > 1e-6:
        raise RuntimeError("Structure-factor K-point power does not reproduce P_mu=-3 to tolerance "
                           "-- STOPPING to diagnose normalization/indexing per task instruction.")
    print("VERIFIED: exact-K structure-factor power reproduces P_mu=-3 to numerical tolerance.")

    df.to_csv(OUT_DIR / "structure_factor_per_replicate.csv", index=False)

    summary_rows = []
    for sigma in SIGMAS:
        g = df[df.sigma == sigma]
        summary_rows.append({
            "sigma": sigma,
            "k_power_fraction_mean": g["k_point_power_fraction"].mean(),
            "k_power_fraction_sd": g["k_point_power_fraction"].std(),
            "k_power_fraction_median": g["k_point_power_fraction"].median(),
            "spectral_entropy_mean": g["spectral_entropy"].mean(),
            "spectral_entropy_sd": g["spectral_entropy"].std(),
            "spectral_entropy_median": g["spectral_entropy"].median(),
            "fraction_replicates_peak_at_K": g["peak_is_K_point"].mean(),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "structure_factor_summary.csv", index=False)
    print("\n" + summary_df.to_string(index=False))

    np.savez(OUT_DIR / "structure_factor_examples.npz",
            sigmas=np.array(SIGMAS),
            **{f"S_sigma_{i}": example_fields[s] for i, s in enumerate(SIGMAS)})

    print(f"\nwrote {OUT_DIR / 'structure_factor_summary.csv'}, "
         f"{OUT_DIR / 'structure_factor_per_replicate.csv'}, "
         f"{OUT_DIR / 'structure_factor_examples.npz'}")


if __name__ == "__main__":
    main()
