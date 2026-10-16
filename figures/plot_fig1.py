"""Figure 1: circuit, contact geometry, and the graph spectrum.

Panel A is left empty on purpose: the author draws the circuit by hand
(Notch-Delta-Figure-1.pptx) and lays it over this frame. The equations are
in the main text, not in the figure. Panel D defines the graph eigenvalue mu
and the plane-wave indices (m, n) used throughout.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import RegularPolygon

from figures.hex_plot import HEX_CAPTION_NOTE, sublattice_labels
from figures.style import BLUE, DIVERGING, INK, ORANGE, YELLOW, apply_style, panel_label, save
from geometry.graph_spectrum import circulant_eigenvalue
from geometry.lattice import COLS, NEIGHBOR_OFFSETS, ROWS, hex_centers

CELL_GREY, NEIGHBOUR = "#d6d5ce", "#f6c3ad"


def _hex_axes(ax, x, y):
    ax.set_xlim(x.min() - 1, x.max() + 1); ax.set_ylim(y.min() - 1, y.max() + 1)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_A(ax):
    ax.axis("off")
    ax.set_title("Feedback circuit (cell $i$, neighbour $j$)", loc="left")


def panel_B(ax):
    rows, cols, cr, cc = 7, 7, 3, 3
    x, y = hex_centers(rows, cols)
    for r in range(rows):
        for c in range(cols):
            colour = ORANGE if (r, c) == (cr, cc) else (
                NEIGHBOUR if (r - cr, c - cc) in NEIGHBOR_OFFSETS else CELL_GREY)
            ax.add_patch(RegularPolygon((x[r, c], y[r, c]), 6, radius=0.58, facecolor=colour,
                                        edgecolor="white", linewidth=1))
    _hex_axes(ax, x, y)
    ax.set_title("Six contact neighbours per cell", loc="left")


def panel_C(ax):
    rows, cols = 9, 9
    labels = sublattice_labels(rows, cols)
    x, y = hex_centers(rows, cols)
    colours = {0: BLUE, 1: ORANGE, 2: YELLOW}
    for r in range(rows):
        for c in range(cols):
            ax.add_patch(RegularPolygon((x[r, c], y[r, c]), 6, radius=0.58, facecolor=colours[int(labels[r, c])],
                                        edgecolor="white", linewidth=1))
    _hex_axes(ax, x, y)
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
    ax.set_title(r"Spectrum of the contact graph ($24\times24$)", loc="left")


def main():
    apply_style()
    fig, grid = plt.subplots(2, 2, figsize=(7.0, 6.6), gridspec_kw={"wspace": 0.3, "hspace": 0.35})
    axes = list(grid.flat)
    panel_A(axes[0]); panel_B(axes[1]); panel_C(axes[2]); panel_D(axes[3])
    for ax, letter in zip(axes, "ABCD"):
        panel_label(ax, letter, x=-0.12)
    fig.text(0.5, -0.01, HEX_CAPTION_NOTE, ha="center", fontsize=7, style="italic", color="#52514e")
    save(fig, "Fig1_circuit_and_geometry")


if __name__ == "__main__":
    main()
