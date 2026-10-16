"""Deterministic integration of the full D-N-A-S lattice and of the matched
3-variable DNA reduction (model/reduced.py) on the periodic 24 x 24
six-neighbour contact lattice.

The drift is forward-invariant on the nonnegative orthant, so no clipping
is applied. The integrator is adaptive DOP853 with rtol = 1e-9 and atol = 1e-11.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from geometry.lattice import COLS, ROWS, hex_adj
from model.equations import drift
from model.params import NotchDeltaParams
from model.reduced import drift_DNA


def pack(*fields: np.ndarray) -> np.ndarray:
    return np.concatenate([f.ravel() for f in fields])


def unpack(y: np.ndarray, n_species: int = 4):
    return [f.reshape(ROWS, COLS) for f in np.split(y, n_species)]


def rhs_lattice(t: float, y: np.ndarray, p: NotchDeltaParams) -> np.ndarray:
    D, N, A, S = unpack(y, 4)
    return pack(*drift(D, N, A, S, p, hex_adj(D), hex_adj(N)))


def rhs_lattice_DNA(t: float, y: np.ndarray, p: NotchDeltaParams) -> np.ndarray:
    D, N, A = unpack(y, 3)
    return pack(*drift_DNA(D, N, A, p, hex_adj(D), hex_adj(N)))


def integrate(y0: np.ndarray, p: NotchDeltaParams, t_final: float = 600.0,
              record_interval: float = 0.5, model: str = "full",
              rtol: float = 1e-9, atol: float = 1e-11):
    rhs = {"full": rhs_lattice, "DNA": rhs_lattice_DNA}[model]
    t_eval = np.arange(0.0, t_final + 1e-9, record_interval)
    return solve_ivp(rhs, (0.0, t_final), y0, args=(p,), method="DOP853",
                     rtol=rtol, atol=atol, t_eval=t_eval)
