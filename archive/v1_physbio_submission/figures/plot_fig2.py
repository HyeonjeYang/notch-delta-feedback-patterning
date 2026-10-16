"""Figure 2 -- What controls patternability? Reads only the small,
already-reduced CSVs in figures/data/ (no new computation)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent / "data"
OUT_DIR = Path(__file__).resolve().parent / "main"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 10.5, "figure.facecolor": "white",
                      "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300})

STRONG = ["d_A", "lam_N", "Kd", "n_ASC", "d_ASC", "lam_ASC"]
WEAK = ["k_cis", "lam_D", "f_D", "f_N", "d_D", "d_N"]

PARAM_TEX = {
    "d_A": r"$d_A$", "lam_N": r"$\lambda_N$", "Kd": r"$K_d$", "n_ASC": r"$n_{ASC}$",
    "d_ASC": r"$d_{ASC}$", "lam_ASC": r"$\lambda_{ASC}$", "k_cis": r"$k_{cis}$",
    "lam_D": r"$\lambda_D$", "f_D": r"$f_D$", "f_N": r"$f_N$", "d_D": r"$d_D$", "d_N": r"$d_N$",
}


def main():
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    df = pd.read_csv(DATA_DIR / "local_sensitivity.csv").sort_values("S_nonuniform", key=abs)
    colors = ["#C0392B" if v > 0 else "#2E86AB" for v in df["S_nonuniform"]]
    ax.barh(df["parameter"], df["S_nonuniform"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(r"$S_i = d\lambda_{pattern} / d\log(p_i)$")
    ax.set_title("A. Ranked local sensitivity at baseline")

    ax = axes[0, 1]
    df1d = pd.read_csv(DATA_DIR / "one_dimensional_scans.csv")
    for name, color in [("lam_N", "#2E86AB"), ("d_A", "#C0392B")]:
        sub = df1d[(df1d.parameter == name) & (df1d.tier.isin(["tier1", "tier2"]))].sort_values("multiplier")
        ax.plot(sub["multiplier"], sub["nonuniform_max_re_lambda"], "o-", ms=3, color=color, label=PARAM_TEX[name])
    sub = df1d[(df1d.parameter == "n_ASC") & (df1d.tier == "hill_set")].sort_values("value")
    ax2 = ax.twiny()
    ax2.plot(sub["value"], sub["nonuniform_max_re_lambda"], "s--", ms=4, color="#16A085",
            label=r"$n_{ASC}$ (top axis)")
    ax2.set_xlabel(r"$n_{ASC}$ value", color="#16A085")
    ax.axhline(0, color="#999", lw=1, ls="--")
    ax.axvline(1.0, color="#888", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"multiplier of baseline ($\lambda_N$, $d_A$)")
    ax.set_ylabel(r"nonuniform max Re($\lambda$)")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, frameon=False, loc="upper left")
    ax.set_title("B. Representative 1-D continuations")

    ax = axes[1, 0]
    df2d = pd.read_csv(DATA_DIR / "two_dimensional_phase_maps.csv")
    sub = df2d[df2d.pair == "lam_N_x_d_A"]
    piv = sub.pivot(index="y", columns="x", values="nonuniform_max_re_lambda")
    X, Y = np.meshgrid(piv.columns.values, piv.index.values)
    vmax = np.nanmax(np.abs(piv.values))
    pc = ax.pcolormesh(X, Y, piv.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="nearest")
    ax.contour(X, Y, piv.values, levels=[0], colors="black", linewidths=1.5)
    ax.plot(1.0, 1.0, "o", color="black", ms=9, mfc="white", mew=2, label="baseline")
    ax.plot(2.0, 1.0, "^", color="black", ms=9, mfc="gold", mew=1.5, label="P2")
    ax.plot(5.623, 1.0, "s", color="black", ms=9, mfc="orange", mew=1.5, label="P3")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\lambda_N$"); ax.set_ylabel(r"$d_A$")
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    fig.colorbar(pc, ax=ax, label=r"nonuniform max Re($\lambda$)", fraction=0.046)
    ax.set_title(r"C. Phase diagram: $\lambda_N \times d_A$")

    ax = axes[1, 1]
    strong_vals = df.set_index("parameter").loc[STRONG, "S_nonuniform"].abs()
    weak_vals = df.set_index("parameter").loc[WEAK, "S_nonuniform"].abs()
    positions = list(range(len(STRONG))) + list(range(len(STRONG) + 1, len(STRONG) + 1 + len(WEAK)))
    labels = [PARAM_TEX[p] for p in STRONG + WEAK]
    colors_bar = ["#D68910"] * len(STRONG) + ["#95A5A6"] * len(WEAK)
    ax.bar(positions, list(strong_vals.values) + list(weak_vals.values), color=colors_bar)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(r"$|S_i|$ (log-sensitivity magnitude)")
    ax.set_title("D. Downstream regulatory (orange) vs.\nclassical cis/trans (grey) drivers")

    fig.suptitle("Figure 2 -- What controls patternability?", y=1.02, fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Fig2_patternability_drivers.png", dpi=300)
    fig.savefig(OUT_DIR / "Fig2_patternability_drivers.pdf")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'Fig2_patternability_drivers.png'}")


if __name__ == "__main__":
    main()
