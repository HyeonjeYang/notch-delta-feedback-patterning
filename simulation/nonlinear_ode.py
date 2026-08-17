"""Validated nonlinear deterministic simulation: full 24x24-lattice
right-hand side for ``scipy.integrate.solve_ivp``, plus the initial-
condition convention used throughout the final analysis. Copied (logic
unchanged) from the source repository's
``continuous_time_validation/scripts/ode_model.py``.

Positivity note: the drift equations are forward-invariant on the
nonnegative orthant by construction (see docs/CONTROL_THEORY_DERIVATION.md
/ docs/GEOMETRY_DERIVATION.md provenance notes) -- no clipping is applied
here, unlike the (separate, not reused) stochastic Euler-Maruyama code path
in the source repository.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from geometry.lattice import COLS, ROWS, hex_adj
from model.equations import drift
from model.params import NotchDeltaParams


def pack(D: np.ndarray, N: np.ndarray, A: np.ndarray, S: np.ndarray) -> np.ndarray:
    return np.concatenate([D.ravel(), N.ravel(), A.ravel(), S.ravel()])


def unpack(y: np.ndarray):
    D, N, A, S = np.split(y, 4)
    return (D.reshape(ROWS, COLS), N.reshape(ROWS, COLS),
            A.reshape(ROWS, COLS), S.reshape(ROWS, COLS))


def rhs_lattice(t: float, y: np.ndarray, p: NotchDeltaParams) -> np.ndarray:
    """Vectorized full-lattice RHS for solve_ivp. y has shape (4*576,)."""
    D, N, A, S = unpack(y)
    s_D = hex_adj(D)
    s_N = hex_adj(N)
    dD, dN, dA, dS = drift(D, N, A, S, p, s_D, s_N)
    return pack(dD, dN, dA, dS)


def generate_ic(seed: int, rows: int = ROWS, cols: int = COLS, ic_amp: float = 0.04):
    """D=0.4, N=0.5, A=0, S=0.3 plus Normal(0, ic_amp) noise, clipped
    nonnegative -- the production initial-condition convention."""
    rng = np.random.default_rng(seed)
    shape = (rows, cols)
    D = np.clip(0.4 + rng.normal(0, ic_amp, shape), 0, None)
    N = np.clip(0.5 + rng.normal(0, ic_amp, shape), 0, None)
    A = np.clip(0.0 + rng.normal(0, min(ic_amp, 0.01), shape), 0, None)
    S = np.clip(0.3 + rng.normal(0, ic_amp, shape), 0, None)
    return D, N, A, S


def integrate(y0: np.ndarray, p: NotchDeltaParams, t_final: float = 600.0,
             record_interval: float = 0.5, rtol: float = 1e-9, atol: float = 1e-11):
    """Adaptive high-accuracy integration (DOP853) -- the solver validated
    throughout the final analysis (agrees with Euler dt=0.01/dt=0.005 to
    6 significant figures at the validated pattern-forming points)."""
    t_eval = np.arange(0.0, t_final + 1e-9, record_interval)
    return solve_ivp(rhs_lattice, (0.0, t_final), y0, args=(p,), method="DOP853",
                     rtol=rtol, atol=atol, t_eval=t_eval)
