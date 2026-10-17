"""Figure 1: circuit, contact geometry, and the graph spectrum.

Panel A is intentionally blank for the author's separately prepared circuit
schematic. Equations remain in the main text. Panel D defines the graph eigenvalue mu
and the plane-wave indices (m, n) used throughout.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from figures.hex_plot import HEX_CAPTION_NOTE, draw_hex_cells, hexagon_offsets, sublattice_label
from figures.style import DIVERGING, INK, ORANGE, apply_style, panel_label, save
from geometry.graph_spectrum import circulant_eigenvalue
from geometry.lattice import COLS, NEIGHBOR_OFFSETS, ROWS

CELL_GREY, NEIGHBOUR = "#E6EAEB", "#F1C3A7"


def panel_A(ax):
    ax.axis("off")
    ax.set_title("Feedback circuit (cell $i$, neighbour $j$)", loc="left")


def panel_B(ax):
    offsets = hexagon_offsets(3)
    colours = [ORANGE if o == (0, 0) else (NEIGHBOUR if o in NEIGHBOR_OFFSETS else CELL_GREY) for o in offsets]
    draw_hex_cells(ax, offsets, colours)
    ax.set_title("Six contact neighbours per cell", loc="left")


def panel_C(ax):
    offsets = hexagon_offsets(4)
    colours = {0: "#6DA3CE", 1: "#E59B73", 2: "#DEC05D"}
    draw_hex_cells(ax, offsets, [colours[sublattice_label(dr, dc)] for dr, dc in offsets])
    handles = [plt.Rectangle((0, 0), 1, 1, color=colours[i]) for i in range(3)]
    ax.legend(handles, ["sublattice 1", "sublattice 2", "sublattice 3"], loc="upper center",
              bbox_to_anchor=(0.5, 0.02), ncol=3, handlelength=1.0, columnspacing=1.0)
    ax.set_title(r"Three-sublattice mode ($\mu=-3$)", loc="left")


def panel_D(ax):
    mu = np.array([[circulant_eigenvalue(m, n) for n in range(COLS)] for m in range(ROWS)])
    im = ax.imshow(mu, cmap=DIVERGING.reversed(), norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=6),
                   origin="lower", extent=[-0.5, COLS - 0.5, -0.5, ROWS - 0.5])
    ax.grid(False)
    box = dict(fc="white", ec="none", alpha=0.9, pad=1.2)
    arrow = dict(arrowstyle="-", color=INK, lw=0.7)
    ax.plot(0, 0, "o", ms=6, mfc="white", mec=INK, mew=1.2, clip_on=False, zorder=5)
    ax.annotate(r"$\mu=6$", (0, 0), xytext=(3.0, 3.0), fontsize=7.5, color=INK, bbox=box, arrowprops=arrow)
    for m, n in [(8, 16), (16, 8)]:
        ax.plot(n, m, "s", ms=6, mfc="white", mec=INK, mew=1.2, zorder=5)
    ax.annotate(r"$\mu=-3$", (16, 8), xytext=(18.0, 3.0), fontsize=7.5, color=INK, bbox=box, arrowprops=arrow)
    ax.annotate(r"$\mu=-3$", (8, 16), xytext=(2.0, 20.0), fontsize=7.5, color=INK, bbox=box, arrowprops=arrow)
    ax.set_xlabel(r"$n$"); ax.set_ylabel(r"$m$")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03, ticks=[-3, 0, 3, 6])
    cb.set_label(r"graph eigenvalue $\mu(m,n)$")
    cb.ax.grid(False)
    ax.set_title(r"Contact spectrum ($24\times24$)", loc="left")


def main():
    apply_style()
    fig, grid = plt.subplots(2, 2, figsize=(8.0, 7.0), gridspec_kw={"wspace": 0.35, "hspace": 0.40})
    axes = list(grid.flat)
    panel_A(axes[0]); panel_B(axes[1]); panel_C(axes[2]); panel_D(axes[3])
    for ax, letter in zip(axes, "ABCD"):
        panel_label(ax, letter, x=-0.12)
    fig.text(0.5, -0.01, HEX_CAPTION_NOTE, ha="center", fontsize=7, style="italic", color="#52514e")
    save(fig, "Fig1_circuit_and_geometry")


if __name__ == "__main__":
    main()
