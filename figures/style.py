"""Publication style with clear mid-saturation colours and redundant markers.

Colour names retain their API for existing plots. No perceptual-validator
or colour-vision certification is implied by this palette.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

FIG_DIR = Path(__file__).resolve().parent
DATA_DIR = FIG_DIR / "data"
MAIN_DIR = FIG_DIR / "main"
SUPP_DIR = FIG_DIR / "supplementary"

BLUE, ORANGE, AQUA, YELLOW = "#2E78B7", "#D16B42", "#209781", "#CB9A2D"
SERIES = [BLUE, ORANGE, AQUA, YELLOW]
INK, INK_2, MUTED, GRID, AXIS = "#25323C", "#52606A", "#88929A", "#E9ECEF", "#BBC3C9"
CONTEXT = "#A8AFB4"

DIVERGING = LinearSegmentedColormap.from_list(
    "blue_neutral_coral", ["#2F73AD", "#99BED9", "#F8F7F2", "#E5A389", "#C95743"])

SEQUENTIAL = LinearSegmentedColormap.from_list("cream_blue", ["#F3E4AF", "#9AC9C7", "#4A9FAA", "#226F9D"])

PARAM_TEX = {
    "lam_D": r"$\lambda_D$", "lam_N": r"$\lambda_N$", "lam_ASC": r"$\lambda_S$", "Kd": r"$K_d$",
    "n_d": r"$n_D$", "n_ASC": r"$n_S$", "d_D": r"$d_D$", "d_N": r"$d_N$", "d_A": r"$d_A$",
    "d_ASC": r"$d_S$", "f_D": r"$f_D$", "f_N": r"$f_N$", "k_cis": r"$k_c$",
}

TIME_LABEL = r"time $t$ (dimensionless)"


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
        "axes.titlesize": 9.5, "axes.titlepad": 12, "axes.labelsize": 9, "legend.fontsize": 7.5,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.labelpad": 6,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK, "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
        "axes.grid": True, "axes.grid.axis": "y", "axes.axisbelow": True,
        "grid.color": GRID, "grid.linewidth": 0.55,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.minor.visible": False, "ytick.minor.visible": False,
        "lines.linewidth": 1.8, "lines.markersize": 4,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300,
        "legend.frameon": False, "legend.handlelength": 1.7,
        "mathtext.fontset": "dejavusans", "pdf.fonttype": 42, "ps.fonttype": 42,
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
