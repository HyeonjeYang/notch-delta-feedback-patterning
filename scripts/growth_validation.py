"""Nonlinear growth-rate validation on the 24 x 24 lattice.

Runs (adaptive DOP853, rtol=1e-9, atol=1e-11):
  full_P2        full model, lam_N = 2.0,   exact-eigenvector initial state
  full_P2_v1ic   full model, lam_N = 2.0,   v1 componentwise-scaled state (comparison)
  DNA_P2         DNA reduction, lam_N = 2.0, exact eigenvector of its own 3x3 block
  full_P3        full model, lam_N = 5.623, exact-eigenvector initial state
  full_P3_v1ic   full model, lam_N = 5.623, v1 componentwise-scaled state

The exact-eigenvector runs start from relative Delta amplitude 1e-6 in the
(m, n) = (8, 16) mu = -3 plane wave. The legacy comparison retains its
original componentwise scaling, whose Delta amplitude is generally different.

Fitting window, fixed before any run: from t = 0 to the first recorded time
at which the K-mode amplitude of D reaches 1e-3 * D*. That is three decades
of growth for the corrected initial state. The endpoint is trajectory-dependent;
the amplitude rule is prespecified. The slope of ln a_D(t) over the
window is compared with the leading eigenvalue of the corresponding
Jacobian block.

Writes figures/data/growth_timeseries.csv, growth_validation.csv,
pattern_final_fields.npz. Run time is a few minutes.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import LAM_N_P2, LAM_N_P3, P0
from model.reduced import reduce_S
from simulation.nonlinear_ode import integrate, pack, unpack
from simulation.perturbation import (eigenmode_state, k_mode_amplitude, leading_real_eigvec,
                                     v1_componentwise_state)
from stability.mode_scan import DOMINANT_PATTERN_MU

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"
REL_AMP = 1e-6
WINDOW_END_REL_AMP = 1e-3
SPECIES = ("D", "N", "A", "S")

RUNS = [
    ("full_P2", "full", LAM_N_P2, "eigenmode", 400.0),
    ("full_P2_v1ic", "full", LAM_N_P2, "v1", 600.0),
    ("DNA_P2", "DNA", LAM_N_P2, "eigenmode", 400.0),
    ("full_P3", "full", LAM_N_P3, "eigenmode", 200.0),
    ("full_P3_v1ic", "full", LAM_N_P3, "v1", 300.0),
]


def run(name: str, model: str, lam_N: float, ic: str, t_final: float):
    p = dataclasses.replace(P0, lam_N=lam_N)
    fp = homogeneous_fixed_point(p)
    J = jacobian_block(DOMINANT_PATTERN_MU, fp, p)
    n_sp = 4 if model == "full" else 3
    if model == "DNA":
        J = reduce_S(J)
    lam_pred, v = leading_real_eigvec(J)
    fp_m = fp[:n_sp]
    fields = (eigenmode_state if ic == "eigenmode" else v1_componentwise_state)(fp_m, v, REL_AMP)
    sol = integrate(pack(*fields), p, t_final=t_final, record_interval=0.25, model=model)
    if not sol.success:
        raise RuntimeError(sol.message)

    amps = np.array([[k_mode_amplitude(f) for f in unpack(y, n_sp)] for y in sol.y.T])
    ts = pd.DataFrame({"run": name, "time": sol.t, **{f"a_{s}": amps[:, i] for i, s in enumerate(SPECIES[:n_sp])}})

    rel = amps[:, 0] / fp[0]
    reached = np.flatnonzero(rel >= WINDOW_END_REL_AMP)
    if len(reached) == 0:
        raise RuntimeError(f"{name}: prespecified growth-window threshold was not reached")
    end = int(reached[0])
    if end < 2 or np.any(amps[:end + 1, 0] <= 0):
        raise RuntimeError(f"{name}: too few valid amplitude observations for fitting")
    win = slice(0, end + 1)
    slope = float(np.polyfit(sol.t[win], np.log(amps[win, 0]), 1)[0])
    pure = amps[0, 0] * np.exp(lam_pred * sol.t[win])
    early = sol.t <= 8.0
    t_pattern = float(sol.t[int(np.argmax(rel >= 0.1))]) if np.any(rel >= 0.1) else np.nan
    final = unpack(sol.y[:, -1], n_sp)
    summary = {"run": name, "model": model, "lam_N": lam_N, "initial_state": ic, "rel_amp": REL_AMP,
               "predicted_lambda": lam_pred, "window_start": 0.0, "window_end": float(sol.t[end]),
               "fitted_slope": slope, "relative_slope_error": abs(slope - lam_pred) / lam_pred,
               "max_rel_dev_from_pure_exponential": float(np.max(np.abs(amps[win, 0] / pure - 1))),
               "early_slope_t0_8": float(np.polyfit(sol.t[early], np.log(amps[early, 0]), 1)[0]),
               "t_rel_amp_0.1": t_pattern, "t_final": t_final,
               "final_D_min": float(final[0].min()), "final_D_max": float(final[0].max())}
    return ts, summary, final[0]


def main() -> None:
    series, summaries, finals = [], [], {}
    for spec in RUNS:
        ts, summary, final_D = run(*spec)
        series.append(ts); summaries.append(summary); finals[spec[0]] = final_D
        print(f"{summary['run']:13s} predicted {summary['predicted_lambda']:.6f}  fitted {summary['fitted_slope']:.6f} "
              f"on [0, {summary['window_end']:.2f}]  (rel err {summary['relative_slope_error']:.1e}; "
              f"early [0,8] slope {summary['early_slope_t0_8']:.4f}; t(rel amp 0.1) {summary['t_rel_amp_0.1']:.1f})")
    pd.concat(series).to_csv(DATA_DIR / "growth_timeseries.csv", index=False)
    out = pd.DataFrame(summaries)
    d_full, d_dna = finals["full_P2"], finals["DNA_P2"]
    out["final_D_maxabs_diff_full_vs_DNA_P2"] = float(np.max(np.abs(d_full - d_dna)))
    out.to_csv(DATA_DIR / "growth_validation.csv", index=False)
    np.savez_compressed(DATA_DIR / "pattern_final_fields.npz", D_full_P2=d_full, D_DNA_P2=d_dna)
    print(f"final P2 pattern, full vs DNA: max |D_full - D_DNA| = {np.max(np.abs(d_full - d_dna)):.2e} "
          f"(D range {d_full.min():.4f}..{d_full.max():.4f})")


if __name__ == "__main__":
    main()
