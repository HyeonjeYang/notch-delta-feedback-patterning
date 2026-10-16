"""Shared figure style: one palette and one set of rcParams for every figure.

Categorical slots follow a fixed order (blue, orange, aqua, yellow). This order
was checked with the dataviz palette validator against a white surface:
adjacent and all-pairs CVD separation pass. Aqua and yellow sit below 3:1
contrast, so they always appear with a legend or a direct label and a
distinct marker shape. Diverging maps run blue -> neutral grey -> red.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

FIG_DIR = Path(__file__).resolve().parent
DATA_DIR = FIG_DIR / "data"
MAIN_DIR = FIG_DIR / "main"
SUPP_DIR = FIG_DIR / "supplementary"

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SERIES = [BLUE, ORANGE, AQUA, YELLOW]
INK, INK_2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
CONTEXT = "#b9b8b1"  # neutral grey for context points (not a series colour)

DIVERGING = LinearSegmentedColormap.from_list(
    "blue_grey_red", ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f2a3a2", "#e34948", "#9c2223"])

PARAM_TEX = {
    "lam_D": r"$\lambda_D$", "lam_N": r"$\lambda_N$", "lam_ASC": r"$\lambda_S$", "Kd": r"$K_d$",
    "n_d": r"$n_D$", "n_ASC": r"$n_S$", "d_D": r"$d_D$", "d_N": r"$d_N$", "d_A": r"$d_A$",
    "d_ASC": r"$d_S$", "f_D": r"$f_D$", "f_N": r"$f_N$", "k_cis": r"$k_c$",
}

TIME_LABEL = r"time $t$ (dimensionless)"


def apply_style() -> None:
    plt.rcParams.update({
        "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 7.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK, "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "lines.linewidth": 1.6, "lines.markersize": 4,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300,
        "legend.frameon": False, "mathtext.fontset": "dejavusans",
    })


def panel_label(ax, letter: str, x: float = -0.16, y: float = 1.06) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left", color=INK)


def save(fig, stem: str, supplementary: bool = False) -> None:
    out = SUPP_DIR if supplementary else MAIN_DIR
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{stem}.png")
    fig.savefig(out / f"{stem}.pdf")
    plt.close(fig)
    print(f"wrote {out / stem}.png/.pdf")
