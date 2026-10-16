"""Hexagonal-cell drawing (polygons, not a square-matrix heatmap).

Cells are drawn in a regular hexagonal window: every cell within hex
distance ``window`` of a centre cell, with lattice indices wrapped
periodically. On the 24 x 24 periodic lattice a window of 11 shows 397
distinct cells and no cell twice. Cell (r, c) sits at x = c + r/2,
y = r sqrt(3)/2, the embedding that matches the six neighbour offsets in
geometry/lattice.py.

Caption note (attach near any figure using this): "Polygons visualize the
six-neighbour contact topology used by the simulation, not measured cell
shape."
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.patches import RegularPolygon

from geometry.lattice import SQRT3_2

HEX_CAPTION_NOTE = ("Polygons visualize the six-neighbour contact topology used by the "
                    "simulation, not measured cell shape.")


def hexagon_offsets(window: int) -> list[tuple[int, int]]:
    """(dr, dc) lattice offsets within hex distance ``window`` of a centre cell."""
    return [(dr, dc) for dr in range(-window, window + 1) for dc in range(-window, window + 1)
            if max(abs(dr), abs(dc), abs(dr + dc)) <= window]


def draw_hex_cells(ax, offsets, colours, radius: float = 0.58, edgecolor: str = "white", linewidth: float = 1.0):
    """Draw one hexagon per (dr, dc) offset, in a regular hexagonal window."""
    xs, ys = [], []
    for (dr, dc), colour in zip(offsets, colours):
        x, y = dc + 0.5 * dr, dr * SQRT3_2
        ax.add_patch(RegularPolygon((x, y), numVertices=6, radius=radius, orientation=0.0,
                                    facecolor=colour, edgecolor=edgecolor, linewidth=linewidth))
        xs.append(x); ys.append(y)
    ax.set_xlim(min(xs) - 1, max(xs) + 1)
    ax.set_ylim(min(ys) - 1, max(ys) + 1)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_hex_field(ax, field: np.ndarray, vmin: float, vmax: float, cmap: str = "viridis",
                   window: int = 11, radius: float = 0.6):
    """Field on a periodic lattice, shown in a regular hexagonal window around the centre."""
    rows, cols = field.shape
    assert 2 * window < min(rows, cols), "window would show some cells twice"
    r0, c0 = rows // 2, cols // 2
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cmap_obj = plt.get_cmap(cmap)
    offsets = hexagon_offsets(window)
    colours = [cmap_obj(norm(field[(r0 + dr) % rows, (c0 + dc) % cols])) for dr, dc in offsets]
    draw_hex_cells(ax, offsets, colours, radius=radius, edgecolor="none", linewidth=0)
    return ScalarMappable(norm=norm, cmap=cmap)


def sublattice_label(dr: int, dc: int, m0: int = 8, n0: int = 16) -> int:
    """Three-sublattice class (0/1/2) of the mu = -3 mode at lattice offset (dr, dc)."""
    return (m0 * dr + n0 * dc) % 3
