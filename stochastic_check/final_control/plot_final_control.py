"""Figure_final_stochastic_control: the one required figure for the final
control task. Reads only already-computed CSV/NPZ files -- no simulation."""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_STOCH = _HERE.parent
_ROOT = _STOCH.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_STOCH))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_DIR = _HERE
COLORS = {"pure_det": "#2E86AB", "random_det": "#27AE60", "stoch_1e-3": "#C0392B"}
LABELS = {"pure_det": "det.\n(pure mu=-3)", "random_det": "det.\n(random IC)",
         "stoch_1e-3": "stochastic\n(sigma=1e-3)"}
ORDER = ["pure_det", "random_det", "stoch_1e-3"]

plt.rcParams.update({"font.size": 10.5, "figure.facecolor": "white",
                      "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300})


def boxplot_panel(ax, groups, col, title, ylabel):
    data = [groups[g][col].dropna().values for g in ORDER]
    bp = ax.boxplot(data, positions=range(3), widths=0.5, patch_artist=True)
    for patch, g in zip(bp["boxes"], ORDER):
        patch.set_facecolor(COLORS[g]); patch.set_alpha(0.6)
    ax.set_xticks(range(3)); ax.set_xticklabels([LABELS[g] for g in ORDER], fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10.5)


def main():
    per_rep = pd.read_csv(_STOCH / "per_replicate_summary.csv")
    pure_det = per_rep[per_rep.sigma == 0.0]
    stoch = per_rep[per_rep.sigma == 1e-3]
    random_det = pd.read_csv(OUT_DIR / "random_ic_summary.csv")
    sf_per_rep = pd.read_csv(OUT_DIR / "structure_factor_per_replicate.csv")
    sf_pure = sf_per_rep[sf_per_rep.sigma == 0.0]
    sf_stoch = sf_per_rep[sf_per_rep.sigma == 1e-3]
    groups = {"pure_det": pure_det, "random_det": random_det, "stoch_1e-3": stoch}

    npz = np.load(OUT_DIR / "structure_factor_examples.npz")
    S_pure = npz["S_sigma_0"]
    S_stoch = npz["S_sigma_2"]  # SIGMAS = [0, 1e-4, 1e-3, 1e-2] -> index 2 == 1e-3
    S_random_det = np.load(OUT_DIR / "random_det_seed0_field.npz")["S_norm"]

    fig = plt.figure(figsize=(14.5, 7.5))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.1])

    boxplot_panel(fig.add_subplot(gs[0, 0]), groups, "t50", "A. Formation time (t50)", "t50")
    boxplot_panel(fig.add_subplot(gs[0, 1]), groups, "final_q_amp", "B. Final pattern amplitude", "final q_amp")
    boxplot_panel(fig.add_subplot(gs[0, 2]), groups, "final_P_mu_minus3", "C. Final mode purity", "final P_mu=-3")

    ax_d = fig.add_subplot(gs[0, 3])
    ax_d.hist(random_det["spectral_entropy"], bins=8, color=COLORS["random_det"], alpha=0.65,
             label="random det", range=(0, 4.2))
    ax_d.hist(sf_stoch["spectral_entropy"], bins=8, color=COLORS["stoch_1e-3"], alpha=0.55,
             label="stochastic 1e-3", range=(0, 4.2))
    ax_d.axvline(sf_pure["spectral_entropy"].iloc[0], color=COLORS["pure_det"], lw=2.5, label="pure det (exact, ln2)")
    ax_d.set_xlabel("spectral entropy"); ax_d.set_ylabel("count")
    ax_d.legend(fontsize=7, frameon=False)
    ax_d.set_title("D. Spectral entropy distributions", fontsize=10.5)

    gs_row2 = gs[1, :3].subgridspec(1, 3)
    axes_sk = [fig.add_subplot(gs_row2[i]) for i in range(3)]
    vmax = max(S_pure.max(), S_stoch.max(), S_random_det.max())
    for ax, S, label in zip(axes_sk, [S_pure, S_random_det, S_stoch],
                            ["pure det (seed n/a)", "random det (seed 0)", "stochastic 1e-3 (rep 0)"]):
        im = ax.imshow(np.fft.fftshift(S), cmap="viridis", vmin=0, vmax=vmax)
        ax.set_title(label, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=axes_sk[-1], fraction=0.046, label="S_norm(k)")

    ax_inset = fig.add_subplot(gs[1, 3])
    ax_inset.axis("off")
    summary_text = (
        "mean (t50 / final q_amp / final P_mu=-3):\n\n"
        f"pure det: {pure_det['t50'].mean():.0f} / {pure_det['final_q_amp'].mean():.4f} / "
        f"{pure_det['final_P_mu_minus3'].mean():.3f}\n\n"
        f"random det: {random_det['t50'].mean():.0f} / {random_det['final_q_amp'].mean():.4f} / "
        f"{random_det['final_P_mu_minus3'].mean():.3f}\n\n"
        f"stochastic 1e-3: {stoch['t50'].mean():.0f} / {stoch['final_q_amp'].mean():.4f} / "
        f"{stoch['final_P_mu_minus3'].mean():.3f}\n\n"
        "Mann-Whitney p (purity,\nrandom det vs. stoch): 0.11\n"
        "Mann-Whitney p (t50,\nrandom det vs. stoch): 3.4e-5"
    )
    ax_inset.text(0.0, 1.0, summary_text, ha="left", va="top", fontsize=8, transform=ax_inset.transAxes,
                 bbox=dict(fc="#F4F6F6", ec="#888", pad=6))

    fig.suptitle("Final stochastic control: pure-mode det. vs. random-IC det. vs. stochastic (P2, lam_N=2.0)",
                y=1.02, fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Figure_final_stochastic_control.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Figure_final_stochastic_control.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'Figure_final_stochastic_control.png'}")


if __name__ == "__main__":
    main()
