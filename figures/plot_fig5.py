"""Figure 5 -- Linear instability becomes a nonlinear pattern (P2, lam_N=2.0)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures.hex_plot import HEX_CAPTION_NOTE, draw_hex_field

DATA_DIR = Path(__file__).resolve().parent / "data"
OUT_DIR = Path(__file__).resolve().parent / "main"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 10.5, "figure.facecolor": "white",
                      "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300})

PREDICTED_RE_LAMBDA = 0.10307  # mu=-3 leading pole at lam_N=2.0 (control/root_locus.py)


TIME_LABEL = r"simulation time $t$ (dimensionless)"


def panel_A(ax, df):
    ax.plot(df["time"], np.maximum(df["q_amp"], 1e-30), color="#2E86AB", lw=1.8)
    ax.axhline(1e-5, color="#333", ls="--", lw=1, label="finite-amplitude gate")
    ax.set_yscale("log")
    ax.set_xlabel(TIME_LABEL); ax.set_ylabel(r"$q_{amp}$ (dimensionless)")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title(r"A. Amplitude vs. time (P2, $\lambda_N=2.0$)", fontsize=10.5)


def panel_B(ax, df):
    sub = df[df.time <= 40]
    amp = np.maximum(sub["std_D"].values, 1e-30)
    t = sub["time"].values
    ax.semilogy(t, amp, "o", ms=3, color="#2E86AB", label=r"observed $\mathrm{std}(D)$")
    fit_mask = (t >= 5) & (t <= 20)
    slope, intercept = np.polyfit(t[fit_mask], np.log(amp[fit_mask]), 1)
    ax.semilogy(t, np.exp(intercept + slope * t), "-", color="#333", lw=1.2, label=fr"fit slope$={slope:.4f}$")
    ax.semilogy(t, amp[0] * np.exp(PREDICTED_RE_LAMBDA * (t - t[0])), "--", color="#C0392B", lw=1.5,
               label=fr"predicted pole$={PREDICTED_RE_LAMBDA:.4f}$")
    ax.set_xlabel(TIME_LABEL); ax.set_ylabel(r"$\mathrm{std}(D)$ (dimensionless)")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("B. Linear growth window: observed vs. predicted", fontsize=10.5)


def panel_C(fig, gs_row, npz):
    times_to_show = [0.0, 40.0, 600.0]
    labels = [r"initial ($t=0$)", r"growth ($t=40$)", r"saturated ($t=600$)"]
    all_t = npz["times"]
    D = npz["D"]
    idxs = [int(np.argmin(np.abs(all_t - t))) for t in times_to_show]
    fields = [D[:, :, i] for i in idxs]
    vmin = min(f.min() for f in fields); vmax = max(f.max() for f in fields)
    axes = [fig.add_subplot(gs_row[i]) for i in range(3)]
    sm = None
    for ax, field, label in zip(axes, fields, labels):
        sm = draw_hex_field(ax, field, vmin, vmax, cmap="viridis")
        ax.set_title(label, fontsize=9.5)
    cax = fig.add_axes([0.93, 0.1, 0.012, 0.22])
    fig.colorbar(sm, cax=cax, label=r"$D$ (dimensionless, common scale)")
    axes[0].text(-0.05, 0.5, "C. Real hexagonal-cell snapshots\n(same color scale)", fontsize=10.5,
               transform=axes[0].transAxes, ha="right", va="center", rotation=90)


def main():
    df = pd.read_csv(DATA_DIR / "nonlinear_figure_data.csv")
    npz = np.load(DATA_DIR / "p2_field_snapshots.npz")

    fig = plt.figure(figsize=(13, 8.5))
    # wspace/hspace give the rotated y-axis labels (e.g. panel B's std(D))
    # room so they don't run into the neighbouring panel's tick labels.
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], wspace=0.45, hspace=0.5)
    ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1])
    panel_A(ax_a, df); panel_B(ax_b, df)

    ax_d = fig.add_subplot(gs[0, 2])
    D_final = npz["D"][:, :, -1].ravel()
    ax_d.hist(D_final, bins=30, color="#8E44AD", alpha=0.8)
    ax_d.axvline(D_final.mean(), color="black", ls="--", lw=1)
    ax_d.set_xlabel(r"$D$ (dimensionless, final, $t=600$)"); ax_d.set_ylabel("cell count")
    ax_d.set_title("D. Bimodal state separation\nat saturation (sender/receiver)", fontsize=10.5)

    gs_row = gs[1, :].subgridspec(1, 3, wspace=0.1)
    panel_C(fig, gs_row, npz)

    fig.suptitle(r"Figure 5 -- Linear instability becomes a nonlinear pattern (P2, $\lambda_N=2.0$)",
                y=1.02, fontsize=14)
    fig.text(0.35, 0.0, HEX_CAPTION_NOTE, ha="center", fontsize=8, style="italic", color="#555")
    fig.savefig(OUT_DIR / "Fig5_nonlinear_pattern.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Fig5_nonlinear_pattern.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'Fig5_nonlinear_pattern.png'}")


if __name__ == "__main__":
    main()
