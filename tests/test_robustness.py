"""Check the new computation against the original unbatched analysis."""
import numpy as np
import pandas as pd

from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import P0, NotchDeltaParams
from scripts.global_sample import analyse, latin_hypercube_parameters
from scripts.robustness_study import (PRIMARY_SEEDS, aggregate, analyse_fast, designs,
                                      kernel_digest, mode_blocks, parameters, protocol, summarize)
from stability.mode_scan import UNIQUE_MU
from stability.sensitivity import PARAM_NAMES
from scripts.validate_robustness_cases import select_cases


def test_protocol_designs_are_unique_and_total_is_correct():
    runs = list(designs())
    assert len(set((s, n, seed) for s, n, seed, _, _ in runs)) == len(runs) == 30
    assert sum(n for _, n, _, _, _ in runs) == protocol()["total_new_points"] == 115000
    assert all(seed != 20260817 for _, _, seed, _, _ in runs)


def test_new_sampler_preserves_original_measure_and_legacy_api():
    np.testing.assert_array_equal(parameters(PRIMARY_SEEDS[0], 32),
                                  latin_hypercube_parameters(PRIMARY_SEEDS[0], 32))
    wide = parameters(PRIMARY_SEEDS[0], 32, .01, 100)
    for j, name in enumerate(PARAM_NAMES):
        lo, hi = (1, 6) if name in ("n_d", "n_ASC") else (.01 * getattr(P0, name), 100 * getattr(P0, name))
        assert np.all((wide[:, j] > lo) & (wide[:, j] < hi))


def test_batched_blocks_and_reductions_match_reference_analysis():
    values = np.vstack([parameters(1730, 5), parameters(810, 5, .01, 100)])
    for v in values:
        p = NotchDeltaParams(**dict(zip(PARAM_NAMES, v)))
        fp = homogeneous_fixed_point(p)
        blocks = np.stack([jacobian_block(mu, fp, p) for mu in UNIQUE_MU])
        np.testing.assert_array_equal(mode_blocks(fp, p), blocks)
        actual, expected = analyse_fast(p), analyse(p)
        for key in ("lambda_pattern", "lambda_DNA", "lambda_DN", "G_minus3", "uniform_re", "D_star", "S_star"):
            np.testing.assert_allclose(actual[key], expected[key], rtol=1e-9, atol=1e-12)


def test_monte_carlo_interval_uses_replicates_not_individual_lhs_points():
    rows = []
    for seed in PRIMARY_SEEDS[:3]:
        values = parameters(seed, 8)
        df = pd.DataFrame([analyse_fast(NotchDeltaParams(**dict(zip(PARAM_NAMES, v)))) for v in values])
        rows.append(summarize(df, "primary", seed, 8))
    frame = pd.DataFrame(rows)
    first = aggregate(frame)["groups"][0]
    frame["n"] = 8000
    second = aggregate(frame)["groups"][0]
    assert first["mc95_low"] == second["mc95_low"]
    assert first["mc95_high"] == second["mc95_high"]
    assert first["replicates"] == second["replicates"] == 3


def test_no_counterexamples_does_not_select_unrelated_rows():
    df = pd.DataFrame({"lambda_pattern": [.1,.2], "lambda_DNA": [1,2], "lambda_DN": [5,10],
                       "G_minus3": [2,3], "n_complex_rhp_K": [0,0], "sample_id": [0,1],
                       "dominance_excess": [0,0]})
    assert select_cases(df) == []


def test_kernel_fingerprint_survives_windows_line_ending_conversion(tmp_path):
    path = tmp_path / "model.py"
    path.write_bytes(b"x = 1\ny = 2\n")
    before = kernel_digest([path], root=tmp_path)
    path.write_bytes(b"x = 1\r\ny = 2\r\n")
    assert kernel_digest([path], root=tmp_path) == before
    path.write_bytes(b"x = 2\r\ny = 2\r\n")
    assert kernel_digest([path], root=tmp_path) != before
