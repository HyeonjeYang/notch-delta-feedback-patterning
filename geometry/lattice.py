"""The periodic six-neighbour triangular contact graph corresponding to a
hexagonal cellular tiling, and its real-space embedding.

``hex_adj`` is copied unmodified (only the six fixed neighbour offsets) from
the source repository's ``core.py::hex_adj`` / ``_neighbor_input`` (the
``contact_sum`` mode -- the only mode used anywhere in the final analysis).
``hex_centers`` is copied unmodified from ``src/utils/hex_lattice.py``.
"""

from __future__ import annotations

import numpy as np

ROWS = 24
COLS = 24

SQRT3_2 = float(np.sqrt(3.0) / 2.0)

# The six offsets hex_adj sums, derived from its np.roll composition (see
# docs/GEOMETRY_DERIVATION.md): cell (r,c) receives contributions from
# (r+dr, c+dc) at each of these six offsets.
NEIGHBOR_OFFSETS = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, -1)]


def hex_adj(arr: np.ndarray) -> np.ndarray:
    """Reciprocal six-neighbour sum on the periodic contact graph (raw,
    unnormalized sum -- not an average)."""
    r1 = np.roll(arr, 1, axis=0)
    rm1 = np.roll(arr, -1, axis=0)
    return (
        r1 + np.roll(r1, -1, axis=1) +
        np.roll(arr, 1, axis=1) + np.roll(arr, -1, axis=1) +
        rm1 + np.roll(rm1, 1, axis=1)
    )


def hex_centers(rows: int = ROWS, cols: int = COLS) -> tuple[np.ndarray, np.ndarray]:
    """Real (x, y) cell-centre coordinates matching hex_adj's fixed
    neighbour offsets (lattice constant a=1), for plotting."""
    r_idx, c_idx = np.indices((rows, cols), dtype=float)
    x = c_idx + 0.5 * r_idx
    y = r_idx * SQRT3_2
    return x, y
