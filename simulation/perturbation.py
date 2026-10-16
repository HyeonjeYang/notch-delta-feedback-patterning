"""Initial conditions and amplitude read-out for the growth-rate validation.

The perturbation is a single spatial mode with an exact eigenvector
direction:

    x_i(0) = x* + eps * v * cos(phi_i),   phi_i = 2 pi (m r_i / R + n c_i / C)

Here v is the real eigenvector of the leading eigenvalue of J_mu
(mu = -3, (m, n) = (8, 16) on the 24 x 24 lattice). v is scaled as a whole by
one scalar eps = rel_amp * D*/v_D, so the Delta component has relative
amplitude rel_amp and the 4-species direction is exactly v. Linear dynamics then
grow the K-mode amplitude as exp(lambda t) from t = 0, with no transient.

The v1 analysis multiplied v componentwise by the equilibrium values
(x* + rel_amp * diag(x*) v cos phi). That vector is not an eigenvector of J_mu
and excites the other three poles of the block. ``v1_componentwise_state``
reproduces that protocol for the comparison in the Supplementary figure only.
"""

from __future__ import annotations

import numpy as np

from geometry.lattice import COLS, ROWS

MODE_M, MODE_N = 8, 16  # plane-wave indices (m, n) with mu = -3 on the 24 x 24 grid


def mode_phase(rows: int = ROWS, cols: int = COLS, m: int = MODE_M, n: int = MODE_N) -> np.ndarray:
    r = np.arange(rows).reshape(-1, 1)
    c = np.arange(cols).reshape(1, -1)
    return 2 * np.pi * (m * r / rows + n * c / cols)


def leading_real_eigvec(J: np.ndarray) -> tuple[float, np.ndarray]:
    """(lambda, v) for the leading eigenvalue, which must be real."""
    eigvals, eigvecs = np.linalg.eig(J)
    k = int(np.argmax(eigvals.real))
    if abs(eigvals[k].imag) > 1e-12:
        raise ValueError("leading eigenvalue is complex; no single real eigenvector")
    return float(eigvals[k].real), eigvecs[:, k].real


def eigenmode_state(fixed_point: np.ndarray, eigvec: np.ndarray, rel_amp: float) -> list[np.ndarray]:
    """x* + eps v cos(phi): exact eigenvector direction, Delta relative amplitude rel_amp."""
    spatial = np.cos(mode_phase())
    eps = rel_amp * fixed_point[0] / eigvec[0]
    return [x0 + eps * v * spatial for x0, v in zip(fixed_point, eigvec)]


def v1_componentwise_state(fixed_point: np.ndarray, eigvec: np.ndarray, rel_amp: float) -> list[np.ndarray]:
    """v1 protocol (comparison only): x* + rel_amp * diag(x*) * (v/|v|) * cos(phi)."""
    spatial = np.cos(mode_phase())
    v = eigvec / np.linalg.norm(eigvec)
    return [np.clip(x0 + rel_amp * max(x0, 1e-6) * vi * spatial, 0, None)
            for x0, vi in zip(fixed_point, v)]


def k_mode_amplitude(field: np.ndarray) -> float:
    """Amplitude of the (m, n) = (8, 16) plane-wave component of one field
    (equals eps * v_X for the initial state above)."""
    dev = field - field.mean()
    return float(2.0 / field.size * np.abs(np.sum(dev * np.exp(-1j * mode_phase(*field.shape)))))
