"""Graph spectrum / dispersion relation of the periodic six-neighbour
triangular contact lattice. Logic copied unmodified from the source
repository's ``continuous_time_validation/scripts/graph_spectrum.py``.

For the plane wave exp(i(k1*r + k2*c)) on the periodic (rows x cols) index
array (k1 = 2*pi*m/rows, k2 = 2*pi*n/cols), applying the six neighbour
offsets and summing gives, exactly:

    mu(k1, k2) = 2*cos(k1) + 2*cos(k2) + 2*cos(k1 - k2)

Verified numerically against the explicit adjacency matrix to 2.8e-14
(machine precision). mu(Gamma, k=0) = 6 (uniform mode, the graph maximum).
The dominant three-sublattice ("K-like") points, at (m,n) multiples of
rows/3, give mu = -3 exactly (the graph minimum): setting
k1=2*pi/3, k2=4*pi/3 makes all three cosine terms equal -1/2.
"""

from __future__ import annotations

import numpy as np

from geometry.lattice import COLS, NEIGHBOR_OFFSETS, ROWS


def circulant_eigenvalue(m: int, n: int, rows: int = ROWS, cols: int = COLS) -> float:
    k1 = 2 * np.pi * m / rows
    k2 = 2 * np.pi * n / cols
    return 2 * np.cos(k1) + 2 * np.cos(k2) + 2 * np.cos(k1 - k2)


def graph_eigenvalues(rows: int = ROWS, cols: int = COLS) -> np.ndarray:
    """All rows*cols graph eigenvalues, sorted descending (uniform mode
    mu=6 first)."""
    vals = [circulant_eigenvalue(m, n, rows, cols) for m in range(rows) for n in range(cols)]
    return np.array(sorted(vals, reverse=True))


def build_adjacency_matrix(rows: int = ROWS, cols: int = COLS) -> np.ndarray:
    """Explicit (rows*cols)x(rows*cols) adjacency matrix, for direct
    verification against the closed-form dispersion relation."""
    n = rows * cols

    def idx(r, c):
        return (r % rows) * cols + (c % cols)

    A = np.zeros((n, n))
    for r in range(rows):
        for c in range(cols):
            i = idx(r, c)
            for dr, dc in NEIGHBOR_OFFSETS:
                A[i, idx(r + dr, c + dc)] += 1.0
    return A


def verify_dispersion(rows: int = ROWS, cols: int = COLS) -> float:
    """Return the max abs deviation between the closed-form dispersion
    relation and the explicit adjacency matrix's eigenvalues."""
    eig_matrix = np.sort(np.linalg.eigvalsh(build_adjacency_matrix(rows, cols)))[::-1]
    eig_formula = graph_eigenvalues(rows, cols)
    return float(np.max(np.abs(eig_matrix - eig_formula)))
