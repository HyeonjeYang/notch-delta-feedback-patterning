"""Targeted numerical checks for reported classifications and growth validation.

This is robustness to numerical settings, not robustness to a new parameter
sampling seed. It reads the published inputs and writes a separate report.
"""
from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from control.state_space import nyquist_rhp_count
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams, P0
from model.reduced import reduce_S
from simulation.nonlinear_ode import integrate, pack, unpack
from simulation.perturbation import eigenmode_state, k_mode_amplitude, leading_real_eigvec
from stability.sensitivity import PARAM_NAMES, local_sensitivity, stage_gain_sensitivity


def main():
    data = ROOT / "figures/data"
    params = pd.read_csv(data / "global_sample_parameters.csv", float_precision="round_trip").set_index("sample_id")
    results = pd.read_csv(data / "global_sample_results.csv").set_index("sample_id")
    # All gain mismatches and the 16 points nearest G=1; union avoids duplicates.
    ids = sorted(set(results.index[(results.G_minus3 > 1) != (results.lambda_pattern > 0)]) |
                 set((results.G_minus3 - 1).abs().nsmallest(16).index))
    nyquist = []
    for sample_id in ids:
        p = NotchDeltaParams(**{k: float(params.loc[sample_id, k]) for k in PARAM_NAMES})
        J = jacobian_block(-3, homogeneous_fixed_point(p), p)
        expected = int((np.linalg.eigvals(J).real > 0).sum())
        counts = {str(n): nyquist_rhp_count(J, n=n) for n in (2000, 4000, 8000)}
        if any(c != expected for c in counts.values()):
            raise RuntimeError(f"Nyquist grid does not resolve sample {sample_id}")
        nyquist.append({"sample_id": int(sample_id), "expected_rhp": expected, "counts": counts})

    base = pd.read_csv(data / "growth_validation.csv").set_index("run")
    growth = []
    for name, model, lam_N, end in (("full_P2", "full", 2, 70), ("DNA_P2", "DNA", 2, 5),
                                     ("full_P3", "full", 5.623, 34)):
        p = dataclasses.replace(P0, lam_N=lam_N)
        fp = homogeneous_fixed_point(p)
        J = jacobian_block(-3, fp, p)
        n_sp = 4 if model == "full" else 3
        predicted, v = leading_real_eigvec(J if model == "full" else reduce_S(J))
        fields = eigenmode_state(fp[:n_sp], v, 1e-6)
        # Tighten both tolerances and halve the recording interval.
        sol = integrate(pack(*fields), p, t_final=end, record_interval=0.125,
                        model=model, rtol=1e-11, atol=1e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        amps = np.array([k_mode_amplitude(unpack(y, n_sp)[0]) for y in sol.y.T])
        hit = np.flatnonzero(amps >= 1e-3 * fp[0])
        if len(hit) == 0 or np.any(amps <= 0):
            raise RuntimeError(f"growth window not reached or amplitude invalid: {name}")
        stop = int(hit[0]) + 1
        slope = float(np.polyfit(sol.t[:stop], np.log(amps[:stop]), 1)[0])
        error = abs(slope / predicted - 1)
        if error >= 5e-4:
            raise RuntimeError(f"growth validation failed at tighter tolerance: {name}, {error}")
        growth.append({"run": name, "predicted_lambda": predicted,
                       "baseline_fitted_slope": float(base.loc[name, "fitted_slope"]),
                       "tighter_fitted_slope": slope, "relative_error_vs_prediction": error,
                       "window_end": float(sol.t[stop - 1]), "rtol": 1e-11, "atol": 1e-14,
                       "record_interval": 0.125})
    sensitivity = {}
    for name, function, h in (("growth", local_sensitivity, 0.01), ("static", stage_gain_sensitivity, 1e-4)):
        def vector(step):
            out = function(P0, h=step)
            return np.array([out[k] for k in PARAM_NAMES]) if isinstance(out, dict) else np.array([r["dlnG_total"] for r in out])
        original = vector(h)
        sensitivity[name] = {str(step): float(np.max(np.abs(vector(step) - original)) / np.max(np.abs(original)))
                             for step in (h / 2, h * 2)}
    report = {"nyquist_n_positive_frequency_points": [2000, 4000, 8000], "nyquist": nyquist,
              "growth_tolerance_checks": growth, "sensitivity_max_change_relative_to_largest_component": sensitivity,
              "global_min_abs_growth": float(results.lambda_pattern.abs().min()),
              "global_unstable_by_absolute_threshold": {
                  str(t): int((results.lambda_pattern > t).sum()) for t in (0, 1e-10, 1e-8, 1e-6)}}
    (data / "numerical_robustness.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "nyquist"}, indent=2))
    print(f"Nyquist grid checks passed at {len(ids)} points")


if __name__ == "__main__":
    main()
