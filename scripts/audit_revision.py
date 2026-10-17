"""Independently audit stored statistics and retain the complete per-mode scan.

Unlike the dominant-mode summaries, this records whether subdominant modes
are also unstable. Run with ``python scripts/audit_revision.py``.
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

from control.state_space import G_exact, nyquist_rhp_count
from model.fixed_point import homogeneous_fixed_point, rhs_single_cell
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams, P0
from model.reduced import reduce_AS, reduce_S
from stability.mode_scan import UNIQUE_MU
from stability.sensitivity import PARAM_NAMES


def main():
    data = ROOT / "figures" / "data"
    params = pd.read_csv(data / "global_sample_parameters.csv", float_precision="round_trip")
    stored = pd.read_csv(data / "global_sample_results.csv").set_index("sample_id")
    nonuniform = UNIQUE_MU < 6 - 1e-9
    other = nonuniform & (UNIQUE_MU > -3 + 1e-9)
    rows, errors, residuals = [], [], []
    counts = dict(unstable=0, gain_mismatches=0, nyquist_mismatches=0,
                  other_nonuniform_unstable=0, uniform_unstable=0,
                  DNA_only=0, DN_only=0, DNA_missed=0, DN_missed=0)
    ratios = []
    for _, row in params.iterrows():
        p = NotchDeltaParams(**{k: float(row[k]) for k in PARAM_NAMES})
        fp = homogeneous_fixed_point(p)
        blocks = np.stack([jacobian_block(mu, fp, p) for mu in UNIQUE_MU])
        ev = np.linalg.eigvals(blocks)
        lead = ev.real.max(axis=1)
        full = float(lead[nonuniform].max())
        dna = float(np.linalg.eigvals(np.stack([reduce_S(J) for J in blocks])).real[nonuniform].max())
        dn = float(np.linalg.eigvals(np.stack([reduce_AS(J) for J in blocks])).real[nonuniform].max())
        gain = G_exact(blocks[0])
        rhp = int((ev[0].real > 0).sum())
        winding = nyquist_rhp_count(blocks[0])
        unstable = full > 0
        counts["unstable"] += int(unstable)
        counts["gain_mismatches"] += int((gain > 1) != unstable)
        counts["nyquist_mismatches"] += int(winding != rhp)
        counts["other_nonuniform_unstable"] += int(np.any(lead[other] > 0))
        counts["uniform_unstable"] += int(lead[-1] > 0)
        for name, value in (("DNA", dna), ("DN", dn)):
            counts[f"{name}_only"] += int(value > 0 and not unstable)
            counts[f"{name}_missed"] += int(value <= 0 and unstable)
        if unstable and dna > 0:
            ratios.append(dna / full)
        saved = stored.loc[int(row.sample_id)]
        errors.append(max(abs(full - saved.lambda_pattern), abs(dna - saved.lambda_DNA),
                          abs(dn - saved.lambda_DN), abs(gain - saved.G_minus3)))
        residuals.append(float(np.max(np.abs(rhs_single_cell(fp, p))) / max(1, np.max(fp))))
        for mu, rate in zip(UNIQUE_MU, lead):
            rows.append((int(row.sample_id), mu, rate))
    p2 = dataclasses.replace(P0, lam_N=2)
    fp2 = homogeneous_fixed_point(p2)
    p2_modes = [(float(mu), float(np.linalg.eigvals(jacobian_block(mu, fp2, p2)).real.max()))
                for mu in UNIQUE_MU]
    report = {"n_samples": len(params), **counts,
              "unstable_percent": 100 * counts["unstable"] / len(params),
              "max_stored_statistic_difference": max(errors),
              "max_relative_equilibrium_residual": max(residuals),
              "DNA_full_growth_ratio_median": float(np.median(ratios)),
              "DNA_full_growth_ratio_min": min(ratios),
              "min_abs_stored_lambda_pattern": float(stored.lambda_pattern.abs().min()),
              "min_abs_stored_G_minus_one": float((stored.G_minus3 - 1).abs().min()),
              "unstable_counts_by_absolute_tolerance": {
                  str(tol): int((stored.lambda_pattern > tol).sum()) for tol in (0, 1e-10, 1e-8, 1e-6)},
              "P2_unstable_distinct_mu": [mu for mu, rate in p2_modes if rate > 0]}
    pd.DataFrame(rows, columns=["sample_id", "mu", "leading_re"]).to_csv(data / "audit_mode_scan.csv", index=False)
    (data / "revision_audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if max(errors) >= 1e-9 or counts["nyquist_mismatches"]:
        raise RuntimeError("stored statistics or Nyquist counts differ from recalculation")


if __name__ == "__main__":
    main()
