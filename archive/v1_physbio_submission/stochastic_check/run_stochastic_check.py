"""One minimal, prespecified stochastic check at the validated P2 point
(lam_N=2.0). Uses the repository's existing frozen model equations
(model.equations.drift, unmodified) and fixed-point/mode machinery
(model.fixed_point, stability.mode_scan) -- no new deterministic result,
no parameter retuning, no parameter search.

Additive white noise on Delta only:
    dD_i = F_D,i(x) dt + sigma * D_star * dW_i
    dN_i = F_N,i(x) dt   (no noise)
    dA_i = F_A,i(x) dt   (no noise)
    dS_i = F_S,i(x) dt   (no noise)

Euler-Maruyama increment: sigma * D_star * sqrt(dt) * Normal(0,1) per cell
per step. Four prespecified sigma values (0, 1e-4, 1e-3, 1e-2), 12
replicates each, dt=0.005, T=600 -- plus one dt=0.0025 convergence control
at sigma=1e-3, replicate labels 0-3 (compared distributionally against the
matched replicate labels at dt=0.005, not pathwise -- independent noise
draws at each dt, same convention as the repository's earlier
numerical-convergence work).

Batched across replicates (one vectorized Euler-Maruyama loop per sigma
group, shape (n_replicates, 24, 24)) for speed; the drift function itself
(model.equations.drift) and the neighbour-sum construction are otherwise
byte-for-byte the same operations the frozen deterministic code uses.
"""

from __future__ import annotations

import dataclasses
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy import stats

from geometry.lattice import COLS, ROWS
from model.equations import drift
from model.fixed_point import find_fixed_points
from model.metrics import morans_i, order_parameters
from model.params import LAM_N_P2, P0
from simulation.perturbation import dominant_eigvec, perturbed_field
from stability.mode_scan import DOMINANT_PATTERN_MU

OUT_DIR = Path(__file__).resolve().parent
REL_AMP = 1e-6
SIGMAS = [0.0, 1e-4, 1e-3, 1e-2]
N_REPLICATES = 12
T_FINAL = 600.0
DT_MAIN = 0.005
DT_CONTROL = 0.0025
CONTROL_SIGMA = 1e-3
CONTROL_REPLICATES = 4
RECORD_INTERVAL = 1.0
T50_THRESHOLD = 0.50
T90_THRESHOLD = 0.90


# ---------------------------------------------------------------- geometry
def hex_adj_batched(arr: np.ndarray) -> np.ndarray:
    """Six-neighbour contact sum for a (batch, rows, cols) stack -- exactly
    geometry.lattice.hex_adj's offsets, batched over axis 0."""
    r1 = np.roll(arr, 1, axis=1)
    rm1 = np.roll(arr, -1, axis=1)
    return (
        r1 + np.roll(r1, -1, axis=2) +
        np.roll(arr, 1, axis=2) + np.roll(arr, -1, axis=2) +
        rm1 + np.roll(rm1, 1, axis=2)
    )


def build_projection_basis(m0: int = 8, n0: int = 16, rows: int = ROWS, cols: int = COLS) -> np.ndarray:
    """Orthonormal (via QR) basis for the 2-D mu=-3 real eigenspace
    (cos/sin of the (m0,n0) plane wave)."""
    r = np.arange(rows).reshape(-1, 1)
    c = np.arange(cols).reshape(1, -1)
    phi = 2 * np.pi * (m0 * r / rows + n0 * c / cols)
    B = np.column_stack([np.cos(phi).ravel(), np.sin(phi).ravel()])
    Q, _ = np.linalg.qr(B)
    return Q  # (rows*cols, 2)


def mode_purity(D: np.ndarray, Q: np.ndarray) -> float:
    """P_mu=-3(t) = |Pi_{-3}(D-mean(D))|^2 / |D-mean(D)|^2."""
    z = (D - D.mean()).ravel()
    denom = float(z @ z)
    if denom <= 0.0:
        return float("nan")
    proj_coeffs = Q.T @ z
    return float((proj_coeffs @ proj_coeffs) / denom)


def sender_receiver_ratio(D: np.ndarray) -> float:
    """Same high/low split convention as model.metrics.order_parameters's
    relative-contrast calculation (D > mean -> sender-like, else receiver-like)."""
    mean = D.mean()
    n_high = int(np.sum(D > mean))
    n_low = int(np.sum(D <= mean))
    return float(n_high) / float(n_low) if n_low > 0 else float("nan")


# ------------------------------------------------------------ integration
def run_batch(sigma: float, n_replicates: int, dt: float, t_final: float,
             D0: np.ndarray, N0: np.ndarray, A0: np.ndarray, S0: np.ndarray,
             p, D_star: float, Q: np.ndarray, seed: int) -> dict:
    """Vectorized Euler-Maruyama for n_replicates independent copies,
    starting from the identical deterministic perturbed field D0/N0/A0/S0
    (broadcast to every replicate)."""
    n_steps = int(round(t_final / dt))
    record_every = max(1, int(round(RECORD_INTERVAL / dt)))
    rng = np.random.default_rng(seed)

    D = np.broadcast_to(D0, (n_replicates, ROWS, COLS)).copy()
    N = np.broadcast_to(N0, (n_replicates, ROWS, COLS)).copy()
    A = np.broadcast_to(A0, (n_replicates, ROWS, COLS)).copy()
    S = np.broadcast_to(S0, (n_replicates, ROWS, COLS)).copy()

    clip_counts = np.zeros(n_replicates, dtype=np.int64)
    times, q_amp_series, std_series, purity_series, moran_series = [], [], [], [], []

    for step in range(n_steps + 1):
        if step % record_every == 0 or step == n_steps:
            t = step * dt
            times.append(t)
            q_amp_row, std_row, purity_row, moran_row = [], [], [], []
            for i in range(n_replicates):
                op = order_parameters(D[i])
                q_amp_row.append(op["q_amp"])
                std_row.append(op["std"])
                moran_row.append(op["moran"])
                purity_row.append(mode_purity(D[i], Q))
            q_amp_series.append(q_amp_row)
            std_series.append(std_row)
            purity_series.append(purity_row)
            moran_series.append(moran_row)
        if step == n_steps:
            break
        s_D = hex_adj_batched(D)
        s_N = hex_adj_batched(N)
        dD, dN, dA, dS = drift(D, N, A, S, p, s_D, s_N)
        noise = sigma * D_star * np.sqrt(dt) * rng.standard_normal(D.shape)
        D_new = D + dD * dt + noise
        clip_counts += np.sum(D_new < 0.0, axis=(1, 2))
        D = np.clip(D_new, 0, None)
        N = np.clip(N + dN * dt, 0, None)
        A = np.clip(A + dA * dt, 0, None)
        S = np.clip(S + dS * dt, 0, None)

    total_cell_steps = ROWS * COLS * n_steps
    clip_fraction = clip_counts / total_cell_steps

    return {
        "times": np.array(times),
        "q_amp": np.array(q_amp_series),      # (n_recorded, n_replicates)
        "std_D": np.array(std_series),
        "P_mu_minus3": np.array(purity_series),
        "moran": np.array(moran_series),
        "clip_fraction": clip_fraction,        # (n_replicates,)
        "final_D": D,                          # (n_replicates, ROWS, COLS)
        "final_sender_receiver": np.array([sender_receiver_ratio(D[i]) for i in range(n_replicates)]),
    }


def first_crossing(times: np.ndarray, series: np.ndarray, threshold: float) -> float:
    idx = np.argmax(series >= threshold) if np.any(series >= threshold) else None
    if idx is None or not np.any(series >= threshold):
        return float("nan")
    return float(times[idx])


def main() -> None:
    t0 = time.time()
    p = dataclasses.replace(P0, lam_N=LAM_N_P2)
    fps = find_fixed_points(p)
    assert len(fps) == 1, f"expected 1 P2 fixed point, found {len(fps)}"
    fp = fps[0]
    D_star = float(fp[0])
    print(f"P2 fixed point: D*={fp[0]:.6f} N*={fp[1]:.6f} A*={fp[2]:.6f} S*={fp[3]:.6f}")

    from stability.mode_scan import lambda_pattern
    stab = lambda_pattern(fp, p)
    print(f"deterministic mu=-3 growth rate: {stab['nonuniform_max_re_lambda']:.5f} "
         f"(expected ~0.10307, mode={stab['nonuniform_best_mu']})")
    assert abs(stab["nonuniform_max_re_lambda"] - 0.10307) < 1e-3
    assert stab["nonuniform_best_mu"] == -3.0

    vec = dominant_eigvec(fp, p, mu=DOMINANT_PATTERN_MU)
    D0, N0, A0, S0 = perturbed_field(fp, vec, REL_AMP)
    Q = build_projection_basis()

    all_rows = []
    summary_rows = []
    seed_base = 1000

    for sigma_idx, sigma in enumerate(SIGMAS):
        t1 = time.time()
        result = run_batch(sigma, N_REPLICATES, DT_MAIN, T_FINAL, D0, N0, A0, S0,
                           p, D_star, Q, seed=seed_base + sigma_idx)
        elapsed = time.time() - t1
        print(f"sigma={sigma:g}: {N_REPLICATES} replicates, dt={DT_MAIN}, T={T_FINAL} "
             f"-> {elapsed:.1f}s, mean clip_fraction={result['clip_fraction'].mean():.3e}")

        for rep in range(N_REPLICATES):
            for k, t in enumerate(result["times"]):
                all_rows.append({
                    "sigma": sigma, "dt": DT_MAIN, "replicate": rep, "time": t,
                    "q_amp": result["q_amp"][k, rep], "std_D": result["std_D"][k, rep],
                    "P_mu_minus3": result["P_mu_minus3"][k, rep], "moran": result["moran"][k, rep],
                })
            t50 = first_crossing(result["times"], result["q_amp"][:, rep], T50_THRESHOLD)
            t90 = first_crossing(result["times"], result["q_amp"][:, rep], T90_THRESHOLD)
            summary_rows.append({
                "sigma": sigma, "dt": DT_MAIN, "replicate": rep,
                "t50": t50, "t90": t90,
                "final_q_amp": result["q_amp"][-1, rep],
                "final_P_mu_minus3": result["P_mu_minus3"][-1, rep],
                "final_moran": result["moran"][-1, rep],
                "final_sender_receiver_ratio": result["final_sender_receiver"][rep],
                "clip_fraction": result["clip_fraction"][rep],
            })

        if sigma == 0.0:
            det_final_q_amp = float(result["q_amp"][-1, 0])
            print(f"  deterministic (sigma=0) reference final q_amp={det_final_q_amp:.6f} "
                 f"(compare to frozen P2 value ~0.9896)")

        if sigma_idx == len(SIGMAS) - 1:
            final_fields = {"sigma0": None, "sigma1e-3": None, "sigma1e-2": None}

    # representative fields for the figure: rerun-free, captured from the
    # main sweep's first replicate at sigma=0, 1e-3, 1e-2
    rep_fields = {}
    for sigma in [0.0, 1e-3, 1e-2]:
        idx = SIGMAS.index(sigma)
        result = run_batch(sigma, 1, DT_MAIN, T_FINAL, D0, N0, A0, S0, p, D_star, Q,
                           seed=seed_base + idx)
        rep_fields[sigma] = result["final_D"][0]

    trajectories_df = pd.DataFrame(all_rows)
    trajectories_df.to_csv(OUT_DIR / "trajectories.csv", index=False)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "per_replicate_summary.csv", index=False)

    # ---- dt convergence control: sigma=1e-3, replicate labels 0-3, dt=0.0025
    t1 = time.time()
    control_result = run_batch(CONTROL_SIGMA, CONTROL_REPLICATES, DT_CONTROL, T_FINAL,
                               D0, N0, A0, S0, p, D_star, Q, seed=9999)
    print(f"dt-convergence control: sigma={CONTROL_SIGMA:g}, {CONTROL_REPLICATES} replicates, "
         f"dt={DT_CONTROL} -> {time.time()-t1:.1f}s")

    control_rows = []
    for rep in range(CONTROL_REPLICATES):
        t50 = first_crossing(control_result["times"], control_result["q_amp"][:, rep], T50_THRESHOLD)
        t90 = first_crossing(control_result["times"], control_result["q_amp"][:, rep], T90_THRESHOLD)
        control_rows.append({
            "sigma": CONTROL_SIGMA, "dt": DT_CONTROL, "replicate": rep,
            "t50": t50, "t90": t90,
            "final_q_amp": control_result["q_amp"][-1, rep],
            "final_P_mu_minus3": control_result["P_mu_minus3"][-1, rep],
            "clip_fraction": control_result["clip_fraction"][rep],
        })
    control_df = pd.DataFrame(control_rows)

    main_sigma1e3 = summary_df[(summary_df.sigma == CONTROL_SIGMA) & (summary_df.replicate < CONTROL_REPLICATES)]
    dt_compare = main_sigma1e3[["replicate", "t50", "t90", "final_q_amp", "final_P_mu_minus3", "clip_fraction"]].copy()
    dt_compare["dt"] = DT_MAIN
    dt_convergence_df = pd.concat([dt_compare, control_df], ignore_index=True)
    dt_convergence_df.to_csv(OUT_DIR / "dt_convergence.csv", index=False)

    # ---- statistical summary
    def summarize(group: pd.DataFrame, col: str) -> dict:
        vals = group[col].dropna().values
        if len(vals) == 0:
            return {"mean": np.nan, "sd": np.nan, "median": np.nan, "iqr_lo": np.nan, "iqr_hi": np.nan, "n": 0}
        return {"mean": float(np.mean(vals)), "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
               "median": float(np.median(vals)),
               "iqr_lo": float(np.percentile(vals, 25)), "iqr_hi": float(np.percentile(vals, 75)),
               "n": len(vals)}

    stat_rows = []
    for sigma in SIGMAS:
        group = summary_df[summary_df.sigma == sigma]
        for col in ["t50", "t90", "final_q_amp", "final_P_mu_minus3", "clip_fraction"]:
            s = summarize(group, col)
            stat_rows.append({"sigma": sigma, "metric": col, **s})
    stat_df = pd.DataFrame(stat_rows)
    stat_df.to_csv(OUT_DIR / "summary.csv", index=False)

    # Kruskal-Wallis, descriptive only
    kw_results = {}
    for col in ["t50", "t90", "final_q_amp", "final_P_mu_minus3"]:
        groups = [summary_df[summary_df.sigma == s][col].dropna().values for s in SIGMAS]
        groups = [g for g in groups if len(g) > 0]
        if len(groups) >= 2 and all(len(g) > 0 for g in groups):
            try:
                stat_val, p_val = stats.kruskal(*groups)
                kw_results[col] = (float(stat_val), float(p_val))
            except ValueError:
                kw_results[col] = (float("nan"), float("nan"))

    print("\nKruskal-Wallis (descriptive only):")
    for col, (stat_val, p_val) in kw_results.items():
        print(f"  {col}: H={stat_val:.3f}, p={p_val:.4f}")

    # ---- figure
    make_figure(trajectories_df, summary_df, rep_fields)

    print(f"\nAll data written to {OUT_DIR}")
    print(f"total elapsed: {time.time()-t0:.1f}s")

    return summary_df, stat_df, kw_results, det_final_q_amp


def make_figure(trajectories_df: pd.DataFrame, summary_df: pd.DataFrame, rep_fields: dict) -> None:
    import matplotlib.pyplot as plt

    sys.path.insert(0, str(OUT_DIR.parents[0]))
    from figures.hex_plot import HEX_CAPTION_NOTE, draw_hex_field

    plt.rcParams.update({"font.size": 10.5, "figure.facecolor": "white",
                         "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300})
    colors = {0.0: "#2E86AB", 1e-4: "#16A085", 1e-3: "#D68910", 1e-2: "#C0392B"}

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1])

    ax_a = fig.add_subplot(gs[0, 0])
    for sigma in SIGMAS:
        sub = trajectories_df[trajectories_df.sigma == sigma]
        piv = sub.pivot(index="time", columns="replicate", values="q_amp")
        mean = piv.mean(axis=1); sd = piv.std(axis=1)
        ax_a.plot(piv.index, mean, color=colors[sigma], label=f"sigma={sigma:g}")
        ax_a.fill_between(piv.index, np.maximum(mean - sd, 1e-30), mean + sd, color=colors[sigma], alpha=0.15)
    ax_a.set_yscale("log")
    ax_a.set_xlabel("time"); ax_a.set_ylabel(r"$q_{amp}$")
    ax_a.legend(fontsize=7.5, frameon=False)
    ax_a.set_title("A. Pattern amplitude (mean +/- SD, 12 replicates)", fontsize=10.5)

    ax_b = fig.add_subplot(gs[0, 1])
    data_t50 = [summary_df[summary_df.sigma == s]["t50"].dropna().values for s in SIGMAS]
    positions = range(len(SIGMAS))
    bp = ax_b.boxplot(data_t50, positions=positions, widths=0.5, patch_artist=True)
    for patch, sigma in zip(bp["boxes"], SIGMAS):
        patch.set_facecolor(colors[sigma]); patch.set_alpha(0.5)
    ax_b.set_xticks(positions); ax_b.set_xticklabels([f"{s:g}" for s in SIGMAS])
    ax_b.set_xlabel("sigma"); ax_b.set_ylabel("t50 (time to q_amp>=0.5)")
    ax_b.set_title("B. Formation time (t50) vs. sigma", fontsize=10.5)

    ax_c = fig.add_subplot(gs[0, 2])
    for sigma in SIGMAS:
        vals = summary_df[summary_df.sigma == sigma]["final_P_mu_minus3"].dropna().values
        ax_c.scatter([sigma] * len(vals), vals, color=colors[sigma], alpha=0.7, s=25)
    ax_c.set_xscale("symlog", linthresh=1e-4)
    ax_c.set_xlabel("sigma"); ax_c.set_ylabel(r"final $P_{\mu=-3}$")
    ax_c.set_title("C. Final mode purity vs. sigma", fontsize=10.5)

    gs_row = gs[1, :].subgridspec(1, 3)
    fields = [rep_fields[0.0], rep_fields[1e-3], rep_fields[1e-2]]
    labels = ["sigma=0", "sigma=1e-3", "sigma=1e-2"]
    vmin = min(f.min() for f in fields); vmax = max(f.max() for f in fields)
    axes_d = [fig.add_subplot(gs_row[i]) for i in range(3)]
    sm = None
    for ax, field, label in zip(axes_d, fields, labels):
        sm = draw_hex_field(ax, field, vmin, vmax, cmap="viridis")
        ax.set_title(label, fontsize=9.5)
    cax = fig.add_axes([0.93, 0.1, 0.012, 0.22])
    fig.colorbar(sm, cax=cax, label="D (common scale)")
    axes_d[0].text(-0.05, 0.5, "D. Final D field\n(same color scale)", fontsize=10.5,
                  transform=axes_d[0].transAxes, ha="right", va="center", rotation=90)

    fig.suptitle("Stochastic check at P2 (lam_N=2.0): weak-to-moderate Delta noise", y=1.02, fontsize=14)
    fig.text(0.35, 0.0, HEX_CAPTION_NOTE, ha="center", fontsize=8, style="italic", color="#555")
    fig.savefig(OUT_DIR / "Figure_stochastic_check.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Figure_stochastic_check.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'Figure_stochastic_check.png'}")


if __name__ == "__main__":
    main()
