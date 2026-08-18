"""Part B+C: deterministic random-initial-condition control.

Tests whether the P_mu=-3 collapse seen under stochastic forcing is a
noise-specific effect, or largely an artifact of comparing against a
symmetry-protected pure-eigenmode deterministic initial condition. Uses
the already-validated deterministic adaptive solver (DOP853,
simulation.nonlinear_ode.integrate, rtol=1e-9, atol=1e-11) -- NO
stochastic forcing anywhere in this script. The only change from the
frozen P2 pure-eigenmode reference is the initial perturbation's SPATIAL
profile: the same 4-component biochemical direction (the dominant mu=-3
eigenvector) is preserved, but the shared spatial factor is replaced by a
mean-zero iid Gaussian field, renormalized so the total perturbation
norm exactly matches the frozen rel_amp=1e-6 pure-mode perturbation.
"""

from __future__ import annotations

import dataclasses
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_STOCH = _HERE.parent
_ROOT = _STOCH.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_STOCH))

import numpy as np
import pandas as pd

from run_stochastic_check import (  # noqa: E402
    COLS, LAM_N_P2, P0, REL_AMP, ROWS, build_projection_basis, dominant_eigvec,
    find_fixed_points, mode_purity, order_parameters, sender_receiver_ratio,
)
from simulation.nonlinear_ode import integrate, pack  # noqa: E402
from analyze_structure_factor import (  # noqa: E402
    k_point_power_fraction, normalize_excluding_dc, peak_location,
    spectral_entropy, structure_factor,
)

OUT_DIR = _HERE
N_SEEDS = 12
T_FINAL = 600.0
RECORD_INTERVAL = 1.0
RTOL, ATOL = 1e-9, 1e-11
T50_THRESHOLD = 0.50
T90_THRESHOLD = 0.90
M0, N0_MODE = 8, 16


def first_crossing(times: np.ndarray, series: np.ndarray, threshold: float) -> float:
    if not np.any(series >= threshold):
        return float("nan")
    return float(times[np.argmax(series >= threshold)])


def random_perturbed_field(fixed_point: np.ndarray, eigvec: np.ndarray, rel_amp: float,
                           seed: int, m0: int = M0, n0: int = N0_MODE,
                           rows: int = ROWS, cols: int = COLS):
    """Same construction as run_stochastic_check.perturbed_field, but with
    the shared spatial factor replaced by a normalized, mean-zero iid
    Gaussian field instead of cos(phi)."""
    D0, N0, A0, S0 = fixed_point
    r = np.arange(rows).reshape(-1, 1)
    c = np.arange(cols).reshape(1, -1)
    phi = 2 * np.pi * (m0 * r / rows + n0 * c / cols)
    cos_profile = np.cos(phi)
    target_norm2 = float(np.sum(cos_profile ** 2))

    rng = np.random.default_rng(seed)
    xi_raw = rng.standard_normal((rows, cols))
    xi_raw = xi_raw - xi_raw.mean()  # exactly mean-zero
    current_norm2 = float(np.sum(xi_raw ** 2))
    xi = xi_raw * np.sqrt(target_norm2 / current_norm2)  # sum(xi^2) == sum(cos_profile^2), exactly

    scale = rel_amp * np.array([D0, max(N0, 1e-6), max(A0, 1e-6), S0])
    D = np.clip(D0 + scale[0] * eigvec[0] * xi, 0, None)
    N = np.clip(N0 + scale[1] * eigvec[1] * xi, 0, None)
    A = np.clip(A0 + scale[2] * eigvec[2] * xi, 0, None)
    S = np.clip(S0 + scale[3] * eigvec[3] * xi, 0, None)

    # verify total perturbation norm matches the pure-mode construction exactly
    pure_scale_sq = np.sum((scale * eigvec) ** 2) * target_norm2
    random_scale_sq = np.sum((scale * eigvec) ** 2) * current_norm2 * (target_norm2 / current_norm2)
    assert abs(pure_scale_sq - random_scale_sq) < 1e-20 * max(pure_scale_sq, 1.0)
    return D, N, A, S, xi


def main() -> None:
    t0 = time.time()
    p = dataclasses.replace(P0, lam_N=LAM_N_P2)
    fp = find_fixed_points(p)[0]
    print(f"P2 fixed point: D*={fp[0]:.6f} N*={fp[1]:.6f} A*={fp[2]:.6f} S*={fp[3]:.6f}")
    vec = dominant_eigvec(fp, p)
    Q = build_projection_basis()

    all_traj_rows = []
    summary_rows = []
    max_clip_check = 0.0

    for seed in range(N_SEEDS):
        D0, N0f, A0, S0, xi = random_perturbed_field(fp, vec, REL_AMP, seed=seed)
        y0 = pack(D0, N0f, A0, S0)

        t1 = time.time()
        sol = integrate(y0, p, t_final=T_FINAL, record_interval=RECORD_INTERVAL, rtol=RTOL, atol=ATOL)
        assert sol.success, f"seed={seed}: DOP853 integration failed"
        elapsed = time.time() - t1

        min_val = sol.y.min()
        max_clip_check = min(max_clip_check, min_val)

        q_amp_list, std_list, purity_list, moran_list = [], [], [], []
        for i, t in enumerate(sol.t):
            D = sol.y[:576, i].reshape(ROWS, COLS)
            op = order_parameters(D)
            q_amp_list.append(op["q_amp"]); std_list.append(op["std"])
            moran_list.append(op["moran"]); purity_list.append(mode_purity(D, Q))
            for key, val in [("q_amp", op["q_amp"]), ("std_D", op["std"]),
                            ("P_mu_minus3", purity_list[-1]), ("moran", op["moran"])]:
                all_traj_rows.append({"seed": seed, "time": t, "metric": key, "value": val})

        q_amp_arr = np.array(q_amp_list)
        purity_arr = np.array(purity_list)
        t50 = first_crossing(sol.t, q_amp_arr, T50_THRESHOLD)
        t90 = first_crossing(sol.t, q_amp_arr, T90_THRESHOLD)

        D_final = sol.y[:576, -1].reshape(ROWS, COLS)
        S_norm = normalize_excluding_dc(structure_factor(D_final))
        k_frac = k_point_power_fraction(S_norm)
        m_peak, n_peak, peak_val = peak_location(S_norm)
        entropy = spectral_entropy(S_norm)

        summary_rows.append({
            "seed": seed, "t50": t50, "t90": t90,
            "final_q_amp": q_amp_arr[-1], "final_P_mu_minus3": purity_arr[-1],
            "final_moran": moran_list[-1],
            "final_sender_receiver_ratio": sender_receiver_ratio(D_final),
            "min_field_value": float(min_val),
            "k_point_power_fraction": k_frac, "spectral_entropy": entropy,
            "peak_m": m_peak, "peak_n": n_peak,
            "peak_is_K_point": bool({(m_peak, n_peak)} <= {(M0, N0_MODE), (N0_MODE, M0)}),
        })
        print(f"seed={seed}: t50={t50}, t90={t90}, final_q_amp={q_amp_arr[-1]:.4f}, "
             f"final_P_mu-3={purity_arr[-1]:.4f}, k_frac={k_frac:.4f}, "
             f"min_field={min_val:.3e} ({elapsed:.1f}s)")

    print(f"\nnumerical check: minimum field value across all 12 runs = {max_clip_check:.3e} "
         f"({'OK, no negative excursion' if max_clip_check >= -1e-9 else 'WARNING: negative value found'})")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "random_ic_summary.csv", index=False)

    traj_long_df = pd.DataFrame(all_traj_rows)
    traj_wide = traj_long_df.pivot_table(index=["seed", "time"], columns="metric", values="value").reset_index()

    n_t = len(np.arange(0.0, T_FINAL + 1e-9, RECORD_INTERVAL))
    q_amp_grid = np.full((N_SEEDS, n_t), np.nan)
    purity_grid = np.full((N_SEEDS, n_t), np.nan)
    std_grid = np.full((N_SEEDS, n_t), np.nan)
    times_ref = None
    for seed in range(N_SEEDS):
        sub = traj_wide[traj_wide.seed == seed].sort_values("time")
        times_ref = sub["time"].values
        q_amp_grid[seed, :len(sub)] = sub["q_amp"].values
        purity_grid[seed, :len(sub)] = sub["P_mu_minus3"].values
        std_grid[seed, :len(sub)] = sub["std_D"].values
    np.savez(OUT_DIR / "random_ic_trajectories.npz", times=times_ref,
            q_amp=q_amp_grid, P_mu_minus3=purity_grid, std_D=std_grid)

    print(f"\nwrote {OUT_DIR / 'random_ic_summary.csv'} and {OUT_DIR / 'random_ic_trajectories.npz'}")
    print(f"total elapsed: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
