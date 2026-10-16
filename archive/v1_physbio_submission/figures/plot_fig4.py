"""Figure 4 -- Geometry selects the spatial mode."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures.hex_plot import HEX_CAPTION_NOTE, draw_hex_field
from geometry.lattice import COLS, ROWS

DATA_DIR = Path(__file__).resolve().parent / "data"
OUT_DIR = Path(__file__).resolve().parent / "main"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 10.5, "figure.facecolor": "white",
                      "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300})


def panel_A(ax):
    ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.2, 2.2); ax.set_aspect("equal")
    theta = np.linspace(0, 2 * np.pi, 7)
    R = 4 * np.pi / 3 / np.sqrt(3)
    bz_x = R * np.cos(theta + np.pi / 6)
    bz_y = R * np.sin(theta + np.pi / 6)
    ax.plot(bz_x, bz_y, "-", color="#555", lw=1.5)
    ax.plot(0, 0, "o", color="black", ms=10)
    ax.annotate(r"$\Gamma$" "\n" r"($\mu=6$)", (0, 0), textcoords="offset points", xytext=(10, 8), fontsize=9)
    K_angles = np.pi / 6 + np.arange(6) * np.pi / 3
    for a in K_angles:
        ax.plot(R * np.cos(a), R * np.sin(a), "s", color="#C0392B", ms=8)
    ax.annotate("K-like\n" r"($\mu=-3$)", (R * np.cos(K_angles[0]), R * np.sin(K_angles[0])),
               textcoords="offset points", xytext=(10, 6), fontsize=9, color="#C0392B")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("A. Brillouin zone of the triangular\ncontact lattice (schematic)", fontsize=10.5)


def panel_B(ax):
    df = pd.read_csv(DATA_DIR / "graph_dispersion.csv")
    piv = df.pivot(index="m", columns="n", values="mu_formula")
    im = ax.imshow(piv.values, cmap="RdBu_r", origin="lower", extent=[0, COLS, 0, ROWS])
    ax.plot(0, 0, "o", color="black", ms=8)
    ax.annotate(r"$\Gamma$ ($\mu=6$)", (0, 0), textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.plot(16, 8, "s", color="lime", ms=8, mec="black")
    ax.annotate(r"K-like ($\mu=-3$)", (16, 8), textcoords="offset points", xytext=(-70, 6), fontsize=8)
    ax.set_xlabel(r"$n$ ($k_2$ index)"); ax.set_ylabel(r"$m$ ($k_1$ index)")
    plt.colorbar(im, ax=ax, label=r"$\mu(k_1,k_2)$", fraction=0.046)
    ax.set_title(r"B. Graph dispersion $\mu(k_1,k_2)$" "\n" r"(numerical, matches closed form to $2.8\times10^{-14}$)",
                fontsize=10.5)


def panel_C(ax):
    df = pd.read_csv(DATA_DIR / "reciprocal_space_growth.csv")
    piv = df.pivot(index="m", columns="n", values="lambda_max")
    vmax = np.nanmax(np.abs(piv.values))
    im = ax.imshow(piv.values, cmap="RdBu_r", origin="lower", extent=[0, COLS, 0, ROWS], vmin=-vmax, vmax=vmax)
    ax.plot(0, 0, "o", color="black", ms=8)
    ax.plot(16, 8, "s", color="lime", ms=8, mec="black")
    ax.set_xlabel("$n$"); ax.set_ylabel("$m$")
    plt.colorbar(im, ax=ax, label=r"$\lambda_{max}(k)$", fraction=0.046)
    ax.set_title(r"C. Growth-rate surface at $\lambda_N=2.0$ (P2)" "\n" r"$k=0$ stable, K-like strongest growth",
                fontsize=10.5)


def panel_D(ax):
    r = np.arange(ROWS).reshape(-1, 1)
    c = np.arange(COLS).reshape(1, -1)
    phi = 2 * np.pi * (8 * r / ROWS + 16 * c / COLS)
    field = np.cos(phi * np.ones((ROWS, COLS)))
    sm = draw_hex_field(ax, field, -1, 1, cmap="RdBu_r")
    plt.colorbar(sm, ax=ax, fraction=0.046, label="mode amplitude")
    ax.set_title(r"D. $\mu=-3$ eigenmode $(m,n)=(8,16)$" "\n" "on the actual hex contact graph", fontsize=10.5)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(12, 10.5))
    panel_A(axes[0, 0]); panel_B(axes[0, 1]); panel_C(axes[1, 0]); panel_D(axes[1, 1])
    fig.suptitle("Figure 4 -- Geometry selects the spatial mode", y=1.02, fontsize=14)
    fig.text(0.5, -0.01, HEX_CAPTION_NOTE, ha="center", fontsize=8, style="italic", color="#555")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Fig4_geometry_selects_mode.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Fig4_geometry_selects_mode.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'Fig4_geometry_selects_mode.png'}")


if __name__ == "__main__":
    main()
