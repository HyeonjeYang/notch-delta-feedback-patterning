"""Real hexagonal-cell visualization utility (not a square-matrix
heatmap), copied unmodified in spirit from the source repository's
``preprint_final_analysis/scripts/hex_plot.py``.

Caption note (attach verbatim near any figure using this): "Polygons
visualize the six-neighbour contact topology used by the simulation, not
measured cell shape."
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.patches import RegularPolygon

from geometry.lattice import hex_centers

HEX_CAPTION_NOTE = ("Polygons visualize the six-neighbour contact topology used by the "
                    "simulation, not measured cell shape.")


def draw_hex_field(ax, field: np.ndarray, vmin: float, vmax: float, cmap: str = "viridis",
                   radius: float = 0.62):
    rows, cols = field.shape
    x, y = hex_centers(rows, cols)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cmap_obj = plt.get_cmap(cmap)
    for r in range(rows):
        for c in range(cols):
            color = cmap_obj(norm(field[r, c]))
            hexagon = RegularPolygon((x[r, c], y[r, c]), numVertices=6, radius=radius,
                                     orientation=0.0, facecolor=color, edgecolor="none")
            ax.add_patch(hexagon)
    ax.set_xlim(x.min() - 1, x.max() + 1)
    ax.set_ylim(y.min() - 1, y.max() + 1)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    return sm


def sublattice_labels(rows: int, cols: int, m0: int = 8, n0: int = 16) -> np.ndarray:
    """Exact three-sublattice class (0/1/2) for the dominant mu=-3 mode."""
    r = np.arange(rows).reshape(-1, 1)
    c = np.arange(cols).reshape(1, -1)
    return (m0 * r + n0 * c) % 3
