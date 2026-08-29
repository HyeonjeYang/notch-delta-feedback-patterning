"""Figure 1 -- Biological circuit and tissue geometry."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

from figures.hex_plot import HEX_CAPTION_NOTE, sublattice_labels
from geometry.lattice import hex_centers

OUT_DIR = Path(__file__).resolve().parent / "main"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 10.5, "figure.facecolor": "white",
                      "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300})


def panel_A(ax):
    # Left blank intentionally -- hand-drawn by the author; keep the
    # coordinate frame so a hand-drawn overlay lines up with panels B-D.
    ax.set_xlim(-2.6, 8.6); ax.set_ylim(-1.2, 4.2); ax.axis("off")
    ax.set_title(r"A. Feedback circuit (cell $i$, neighbour $j$)", fontsize=10.5)


def panel_B(ax):
    rows, cols = 7, 7
    x, y = hex_centers(rows, cols)
    cr, cc = 3, 3
    offsets = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, -1)]
    for r in range(rows):
        for c in range(cols):
            is_center = (r == cr and c == cc)
            is_neighbor = any(r == cr + dr and c == cc + dc for dr, dc in offsets)
            color = "#C0392B" if is_center else ("#F5CBA7" if is_neighbor else "#D6DBDF")
            hexagon = RegularPolygon((x[r, c], y[r, c]), numVertices=6, radius=0.58,
                                     orientation=0.0, facecolor=color, edgecolor="white", linewidth=1)
            ax.add_patch(hexagon)
    ax.set_xlim(x.min() - 1, x.max() + 1); ax.set_ylim(y.min() - 1, y.max() + 1)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("B. Hexagonal tiling: one cell's\nsix contact neighbours", fontsize=10.5)


def panel_C(ax):
    rows, cols = 9, 9
    labels = sublattice_labels(rows, cols)
    x, y = hex_centers(rows, cols)
    colors3 = {0: "#2E86AB", 1: "#C0392B", 2: "#F1C40F"}
    for r in range(rows):
        for c in range(cols):
            hexagon = RegularPolygon((x[r, c], y[r, c]), numVertices=6, radius=0.58,
                                     orientation=0.0, facecolor=colors3[int(labels[r, c])],
                                     edgecolor="white", linewidth=1)
            ax.add_patch(hexagon)
    ax.set_xlim(x.min() - 1, x.max() + 1); ax.set_ylim(y.min() - 1, y.max() + 1)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors3[i]) for i in range(3)]
    ax.legend(handles, ["sublattice A", "sublattice B", "sublattice C"], fontsize=7.5,
             frameon=False, loc="upper right", bbox_to_anchor=(1.15, 1.05))
    ax.set_title("C. Three-sublattice (A/B/C)\ndominant mode structure", fontsize=10.5)


def panel_D(ax):
    ax.axis("off")
    text = (
        r"$\dot D = \lambda_D\,h_D(S) - d_D D - f_D D\,s_N - k_{cis} N D$" "\n"
        r"$\dot N = \lambda_N - d_N N - f_N N\,s_D - k_{cis} N D$" "\n"
        r"$\dot A = -d_A A + f_N N\,s_D$" "\n"
        r"$\dot S = -d_S S + \lambda_S\,h_S(A)$" "\n\n"
        r"$s_D=\sum_{j\sim i}D_j,\ \ s_N=\sum_{j\sim i}N_j$" "\n\n"
        "Small-signal block preview:\n"
        r"$\dot z = A_\mu z + bS,\ \ \dot S = c^{\!\top}z + a_S S$"
    )
    txt = ax.text(0.02, 0.5, text, ha="left", va="center", fontsize=10.5, transform=ax.transAxes,
                 bbox=dict(fc="#F4F6F6", ec="#888", pad=8))
    ax.set_title("D. Governing equations (exact, model.equations)\n+ block preview", fontsize=10.5)
    # Equation lines stay left-aligned relative to each other, but the block
    # as a whole is re-centred under the "D" title (its width isn't known
    # until after a render pass).
    ax.figure.canvas.draw()
    bbox_axes = txt.get_window_extent(renderer=ax.figure.canvas.get_renderer()).transformed(ax.transAxes.inverted())
    txt.set_x(0.5 - (bbox_axes.x1 - bbox_axes.x0) / 2)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 10))
    panel_A(axes[0, 0]); panel_B(axes[0, 1]); panel_C(axes[1, 0]); panel_D(axes[1, 1])
    fig.suptitle("Figure 1 -- Biological circuit and tissue geometry", y=1.02, fontsize=14)
    fig.text(0.5, -0.01, HEX_CAPTION_NOTE, ha="center", fontsize=8, style="italic", color="#555")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Fig1_circuit_and_geometry.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Fig1_circuit_and_geometry.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'Fig1_circuit_and_geometry.png'}")


if __name__ == "__main__":
    main()
