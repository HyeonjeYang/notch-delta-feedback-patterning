"""Replicated LHS study, sample-size checks and a wider-box stress test.

The protocol is written before computation. Inputs and results are retained
per replicate; numerical/biological probabilities are not conflated.
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
import sys
from pathlib import Path

# Small matrix problems: avoid BLAS thread oversubscription in this process.
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from control.state_space import nyquist_rhp_count
from model.equations import hill_ASC, hill_D
from model.fixed_point import homogeneous_fixed_point, rhs_single_cell
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams, P0
from stability.mode_scan import UNIQUE_MU
from stability.sensitivity import PARAM_NAMES

OUT = ROOT / "figures/data/robustness"
PRIMARY_SEEDS = tuple(range(2026101701, 2026101711))
WIDE_SEEDS = tuple(range(2026101801, 2026101806))
TOL = 1e-8
CORE_PATHS = [*sorted((ROOT / "model").glob("*.py")),
              ROOT / "geometry/graph_spectrum.py", ROOT / "geometry/lattice.py",
              ROOT / "stability/mode_scan.py", ROOT / "stability/sensitivity.py",
              ROOT / "control/state_space.py", Path(__file__).resolve()]


def kernel_digest(paths=None, root=ROOT):
    """Canonical text hashing survives a Git LF/CRLF checkout conversion."""
    if paths is None:
        paths = CORE_PATHS
    chunks = [p.relative_to(root).as_posix().encode() + b"\0" + p.read_text(encoding="utf-8").encode()
              for p in paths]
    return hashlib.sha256(b"\0".join(chunks)).hexdigest()


def protocol() -> dict:
    return {
        "protocol_id": "independent-lhs-robustness-v1",
        "primary": {"seeds": list(PRIMARY_SEEDS), "n": 4000, "bounds": [0.1, 10]},
        "size_checks": {"seeds": list(PRIMARY_SEEDS[:5]), "n": [1000, 2000, 8000], "bounds": [0.1, 10]},
        "wider_box": {"seeds": list(WIDE_SEEDS), "n": 4000, "bounds": [0.01, 100]},
        "sampler": {"class": "scipy.stats.qmc.LatinHypercube", "legacy_seed_api": True,
                    "scramble": True, "strength": 1, "optimization": None,
                    "rate_measure": "log-uniform", "Hill_measure": "uniform [1,6]"},
        "parameter_order": PARAM_NAMES,
        "classification": {"main": "strict real part > 0", "marginal_band": TOL,
                           "sensitivity_thresholds": [0, 1e-10, 1e-8, 1e-6],
                           "mode_dominance_excess_tolerance": 1e-9},
        "nyquist": "all static-gain mismatches plus 16 closest-to-unit-gain points per replicate; grid 4000, refine on disagreement",
        "intervals": "approximate 95% Student-t intervals across independent replicate fractions, not iid-binomial intervals across LHS points",
        "size_dependence": "same first five seeds reused across sizes; separately generated non-nested designs; do not pool sizes for inference",
        "counterexamples": "retain all; validate up to five per category selected by largest signed margin, then sample ID, with 80-digit arithmetic",
        "total_new_points": 115000,
    }


def designs():
    for seed in PRIMARY_SEEDS:
        yield "primary", 4000, seed, 0.1, 10.0
    for n in (1000, 2000, 8000):
        for seed in PRIMARY_SEEDS[:5]:
            yield "size", n, seed, 0.1, 10.0
    for seed in WIDE_SEEDS:
        yield "wide", 4000, seed, 0.01, 100.0


def parameters(seed: int, n: int, low: float = 0.1, high: float = 10) -> np.ndarray:
    unit = qmc.LatinHypercube(d=13, scramble=True, strength=1, optimization=None, seed=seed).random(n)
    columns = []
    for i, name in enumerate(PARAM_NAMES):
        if name in ("n_d", "n_ASC"):
            columns.append(1 + 5 * unit[:, i])
        else:
            lo, hi = low * getattr(P0, name), high * getattr(P0, name)
            columns.append(np.exp(np.log(lo) + unit[:, i] * (np.log(hi) - np.log(lo))))
    return np.column_stack(columns)


def mode_blocks(fp, p, mus=UNIQUE_MU):
    """Exact original Jacobian entries, evaluated for all modes together."""
    D, N = fp[:2]
    J = np.broadcast_to(jacobian_block(0, fp, p), (len(mus), 4, 4)).copy()
    J[:, 0, 1] = -p.k_cis * D - p.f_D * D * mus
    J[:, 1, 0] = -p.k_cis * N - p.f_N * N * mus
    J[:, 2, 0] = p.f_N * N * mus
    return J


def analyse_fast(p: NotchDeltaParams) -> dict:
    fp = homogeneous_fixed_point(p)
    J = mode_blocks(fp, p)
    ev = np.linalg.eigvals(J)
    dna = J[:, :3, :3] - (J[:, :3, 3, None] * J[:, None, 3, :3]) / J[:, 3, 3, None, None]
    dn = J[:, :2, :2] - J[:, :2, 2:] @ np.linalg.solve(J[:, 2:, 2:], J[:, 2:, :2])
    edna, edn = np.linalg.eigvals(dna), np.linalg.eigvals(dn)
    lead = ev.real.max(axis=1)
    kd = np.argmax(ev.real, axis=1)
    best = int(np.argmax(lead[:-1]))
    rhs = rhs_single_cell(fp, p)
    D, N, A, S = fp
    production = np.array([p.lam_D * hill_D(S, p), p.lam_N, 6 * p.f_N * N * D, p.lam_ASC * hill_ASC(A, p)])
    balance = float(np.max(np.abs(rhs) / np.maximum(production + np.abs(production - rhs), np.finfo(float).tiny)))
    solved = np.linalg.solve(J[:, :3, :3], J[:, :3, 3, None])[..., 0]
    gains = np.einsum("mi,mi->m", J[:, 3, :3], solved) / J[:, 3, 3]
    return {
        "D_star": D, "N_star": N, "A_star": A, "S_star": S,
        "eq_residual": float(np.max(np.abs(rhs)) / max(1, np.max(fp))),
        "balance_residual": balance,
        "lambda_pattern": float(lead[best]), "best_mu": float(UNIQUE_MU[best]),
        "lambda_pattern_im": float(abs(ev[best, kd[best]].imag)),
        "lambda_K": float(lead[0]), "uniform_re": float(lead[-1]),
        "n_rhp_K": int((ev[0].real > 0).sum()),
        "n_complex_rhp_K": int(((ev[0].real > 0) & (abs(ev[0].imag) > 1e-10)).sum()),
        "min_abs_K_real": float(np.min(np.abs(ev[0].real))),
        "G_minus3": float(gains[0]), "G_argmax_mu": float(UNIQUE_MU[np.argmax(gains)]),
        "dominance_excess": float(lead[best] - lead[0]),
        "plant_max_re": float(np.linalg.eigvals(J[:, :3, :3]).real.max()),
        "other_nonuniform_unstable": bool(np.any(lead[1:-1] > 0)),
        "lambda_DNA": float(edna[:-1].real.max()), "lambda_DN": float(edn[:-1].real.max()),
        "uniform_DNA": float(edna[-1].real.max()), "uniform_DN": float(edn[-1].real.max()),
    }


def summarize(df: pd.DataFrame, scenario: str, seed: int, n: int) -> dict:
    u, udna, udn = df.lambda_pattern > 0, df.lambda_DNA > 0, df.lambda_DN > 0
    mismatch = u != (df.G_minus3 > 1)
    out = {"scenario": scenario, "seed": seed, "n": n, "n_unstable": int(u.sum()),
           "unstable_fraction": float(u.mean()), "n_gain_mismatch": int(mismatch.sum()),
           "gain_mismatch_fraction": float(mismatch.mean()),
           "n_other_mode_dominant": int((u & (df.dominance_excess > 1e-9)).sum()),
           "n_other_unstable_K_stable": int((u & (df.lambda_K <= 0)).sum()),
           "n_gain_max_not_K": int((df.G_argmax_mu != -3).sum()),
           "n_plant_not_hurwitz": int((df.plant_max_re >= 0).sum()),
           "n_near_marginal_full": int((df.lambda_pattern.abs() <= TOL).sum()),
           "min_abs_growth": float(df.lambda_pattern.abs().min()),
           "max_balance_residual": float(df.balance_residual.max()),
           "n_mismatch_not_two_rhp": int((mismatch & (df.n_rhp_K != 2)).sum())}
    for key, uk in (("DNA", udna), ("DN", udn)):
        both = u & uk
        ratio = df.loc[both, f"lambda_{key}"] / df.loc[both, "lambda_pattern"]
        out.update({f"{key}_only": int((uk & ~u).sum()), f"{key}_missed": int((~uk & u).sum()),
                    f"{key}_ratio_median": float(ratio.median()),
                    f"{key}_ratio_q25": float(ratio.quantile(.25)), f"{key}_ratio_q75": float(ratio.quantile(.75)),
                    f"{key}_ratio_min": float(ratio.min()),
                    f"{key}_not_faster": int((ratio <= 1).sum()),
                    f"{key}_slower_confirmable": int((both & (df[f"lambda_{key}"] < df.lambda_pattern - 1e-9)).sum())})
    for tol in (1e-10, 1e-8, 1e-6):
        out[f"n_unstable_tol_{tol}"] = int((df.lambda_pattern > tol).sum())
    return out


def check_nyquist(df: pd.DataFrame, scenario: str, seed: int, n: int):
    mismatch = (df.lambda_pattern > 0) != (df.G_minus3 > 1)
    ids = sorted(set(df.index[mismatch]) | set((df.G_minus3 - 1).abs().nsmallest(16).index))
    rows = []
    for sid in ids:
        r = df.loc[sid]
        p = NotchDeltaParams(**{k: float(r[k]) for k in PARAM_NAMES})
        J = jacobian_block(-3, homogeneous_fixed_point(p), p)
        expected = int((np.linalg.eigvals(J).real > 0).sum())
        checks = {}
        for grid in (4000, 8000, 16000, 32000):
            checks[str(grid)] = nyquist_rhp_count(J, n=grid)
            if checks[str(grid)] == expected:
                break
        rows.append({"scenario": scenario, "seed": seed, "n": n, "sample_id": int(sid),
                     "expected_rhp": expected, "computed_rhp": checks[str(grid)], "grid": grid,
                     "gain_mismatch": bool(mismatch.loc[sid]), "counts": json.dumps(checks)})
    return rows


def aggregate(summary: pd.DataFrame) -> dict:
    groups = []
    for (scenario, n), g in summary.groupby(["scenario", "n"], sort=False):
        values = g.unstable_fraction.values
        mean, sem = float(values.mean()), float(stats.sem(values))
        half = float(stats.t.ppf(.975, len(g) - 1) * sem)
        groups.append({"scenario": scenario, "n": int(n), "replicates": len(g),
                       "mean_unstable_fraction": mean, "mc95_low": mean - half, "mc95_high": mean + half,
                       "min_unstable_fraction": float(values.min()), "max_unstable_fraction": float(values.max()),
                       "DNA_replicate_median_min": float(g.DNA_ratio_median.min()),
                       "DNA_replicate_median_max": float(g.DNA_ratio_median.max()),
                       "DN_replicate_median_min": float(g.DN_ratio_median.min()),
                       "DN_replicate_median_max": float(g.DN_ratio_median.max()),
                       "gain_mismatches_total": int(g.n_gain_mismatch.sum()),
                       "DNA_only_total": int(g.DNA_only.sum()), "DNA_missed_total": int(g.DNA_missed.sum()),
                       "DNA_not_faster_total": int(g.DNA_not_faster.sum()),
                       "DN_only_total": int(g.DN_only.sum()), "DN_missed_total": int(g.DN_missed.sum()),
                       "DN_not_faster_total": int(g.DN_not_faster.sum()),
                       "other_mode_dominant_total": int(g.n_other_mode_dominant.sum()),
                       "near_marginal_total": int(g.n_near_marginal_full.sum())})
    return {"groups": groups, "total_points": int(summary.n.sum()),
            "note": "Intervals describe integration uncertainty under the specified sampling measure, not biological variability. Sizes share some seeds; do not pool their intervals."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol-only", action="store_true")
    parser.add_argument("--reuse", action="store_true", help="reuse checksum-verified replicates from this exact kernel")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    plan = protocol()
    plan_path = OUT / "protocol.json"
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise RuntimeError("Protocol differs from the frozen study; use a new study directory")
    plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    if args.protocol_only:
        print(f"wrote protocol before computation: {plan_path}")
        return
    kernel_hash = kernel_digest()
    summaries, nyquist, failures = [], [], []
    for scenario, n, seed, low, high in designs():
        name = f"{scenario}-n{n}-seed{seed}"
        path, meta_path = OUT / f"{name}.csv.gz", OUT / f"{name}.json"
        inputs = parameters(seed, n, low, high)
        input_hash = hashlib.sha256(inputs.tobytes()).hexdigest()
        if path.exists():
            frozen_inputs = pd.read_csv(path, usecols=PARAM_NAMES, float_precision="round_trip")
            np.testing.assert_allclose(frozen_inputs.values, inputs, rtol=5e-14, atol=0)
        if args.reuse and path.exists() and meta_path.exists():
            meta = json.loads(meta_path.read_text())
            if meta["kernel_sha256"] != kernel_hash or meta["parameter_array_sha256"] != input_hash:
                raise RuntimeError(f"Cached study has different code or inputs: {name}")
            if hashlib.sha256(path.read_bytes()).hexdigest() != meta["csv_sha256"]:
                raise RuntimeError(f"Cached CSV checksum differs: {name}")
            df = pd.read_csv(path, float_precision="round_trip")
            np.testing.assert_allclose(df[PARAM_NAMES].values, inputs, rtol=5e-14, atol=0)
            print(f"reusing {name}", flush=True)
        else:
            records = []
            for sid, values in enumerate(inputs):
                p = NotchDeltaParams(**dict(zip(PARAM_NAMES, values)))
                try:
                    record = analyse_fast(p)
                    if not np.isfinite(list(record.values())).all():
                        raise RuntimeError("non-finite analysis result")
                    records.append({"sample_id": sid, **dict(zip(PARAM_NAMES, values)), **record})
                except Exception as exc:
                    failures.append({"scenario": scenario, "n": n, "seed": seed, "sample_id": sid,
                                     **dict(zip(PARAM_NAMES, values)), "error": str(exc)})
                if (sid + 1) % 1000 == 0:
                    print(f"{name}: {sid + 1}/{n}", flush=True)
            if len(records) != n:
                pd.DataFrame(failures).to_csv(OUT / "failures.csv", index=False)
                raise RuntimeError(f"{name}: solver/analysis failure retained; not labelled no-equilibrium")
            df = pd.DataFrame(records)
            df.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
            meta = {"kernel_sha256": kernel_hash, "parameter_array_sha256": input_hash,
                    "kernel_hash_format": "path-delimited UTF-8 text with universal newlines",
                    "scenario": scenario, "n": n, "seed": seed, "bounds": [low, high],
                    "csv_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        df = df.set_index("sample_id")
        item = summarize(df, scenario, seed, n)
        summaries.append(item)
        nyquist.extend(check_nyquist(df, scenario, seed, n))
        pd.DataFrame(summaries).to_csv(OUT / "replicates.csv", index=False)
        pd.DataFrame(nyquist).to_csv(OUT / "nyquist_checks.csv", index=False)
        print(f"{name}: unstable {item['unstable_fraction']:.3%}, DNA median {item['DNA_ratio_median']:.3f}, "
              f"DNA missed {item['DNA_missed']}, DNA not faster {item['DNA_not_faster']}", flush=True)
    summary = pd.DataFrame(summaries)
    report = aggregate(summary)
    report["nyquist_checked_points"] = len(nyquist)
    report["nyquist_disagreements"] = sum(r["expected_rhp"] != r["computed_rhp"] for r in nyquist)
    report["max_balance_residual"] = float(summary.max_balance_residual.max())
    report["plant_failures"] = int(summary.n_plant_not_hurwitz.sum())
    report["gain_order_failures"] = int(summary.n_gain_max_not_K.sum())
    report["status"] = "complete"
    (OUT / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if report["nyquist_disagreements"]:
        raise RuntimeError("Nyquist count disagreement remains unresolved; retained for review")


if __name__ == "__main__":
    main()
