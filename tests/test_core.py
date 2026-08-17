"""Minimal validation suite: checks that this repository's refactored code
reproduces the exact numbers from the frozen source analysis (not a
re-derivation -- a fidelity check)."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from control.root_locus import find_critical_lam_N
from control.schur_verification import verify as schur_verify
from control.state_space import G_exact
from geometry.graph_spectrum import graph_eigenvalues, verify_dispersion
from model.fixed_point import find_fixed_points
from model.jacobian import jacobian_block
from model.params import LAM_N_P2, LAM_N_P3, P0
from stability.mode_scan import DOMINANT_PATTERN_MU, UNIFORM_MU, lambda_pattern


def test_baseline_fixed_point_matches_frozen_value():
    fps = find_fixed_points(P0)
    assert len(fps) == 1
    D_star = fps[0][0]
    assert D_star == pytest.approx(28.548283450528967, rel=1e-6)


def test_p2_dominant_mode_growth_rate_matches_frozen_value():
    p = dataclasses.replace(P0, lam_N=LAM_N_P2)
    fp = find_fixed_points(p)[0]
    stab = lambda_pattern(fp, p)
    assert stab["nonuniform_max_re_lambda"] == pytest.approx(0.10307, abs=1e-4)
    assert stab["nonuniform_best_mu"] == pytest.approx(-3.0, abs=1e-6)


def test_p3_dominant_mode_growth_rate_matches_frozen_value():
    p = dataclasses.replace(P0, lam_N=LAM_N_P3)
    fp = find_fixed_points(p)[0]
    stab = lambda_pattern(fp, p)
    assert stab["nonuniform_max_re_lambda"] == pytest.approx(0.2163, abs=1e-4)


def test_baseline_is_spatially_stable():
    fp = find_fixed_points(P0)[0]
    stab = lambda_pattern(fp, P0)
    assert stab["nonuniform_max_re_lambda"] < 0
    assert stab["uniform_re_lambda"] < 0


def test_graph_dispersion_matches_adjacency_matrix():
    assert verify_dispersion() < 1e-10


def test_uniform_and_dominant_mode_eigenvalues():
    eigs = graph_eigenvalues()
    assert eigs.max() == pytest.approx(6.0, abs=1e-9)
    assert eigs.min() == pytest.approx(-3.0, abs=1e-9)


def test_schur_complement_matches_direct_jacobian_eigenvalues():
    fp = find_fixed_points(P0)[0]
    for mu in (UNIFORM_MU, DOMINANT_PATTERN_MU):
        J_mu = jacobian_block(mu, fp, P0)
        result = schur_verify(J_mu)
        assert result["max_coeff_diff"] < 1e-10
        assert result["max_root_mismatch"] < 1e-10


def test_critical_lam_n_and_unity_loop_gain():
    lam_N_values = np.geomspace(0.5, 6.5, 61)
    lam_crit = find_critical_lam_N(DOMINANT_PATTERN_MU, lam_N_values, P0)
    assert lam_crit == pytest.approx(1.305, abs=0.05)

    p = dataclasses.replace(P0, lam_N=lam_crit)
    fp = find_fixed_points(p)[0]
    J_mu = jacobian_block(DOMINANT_PATTERN_MU, fp, p)
    assert G_exact(J_mu) == pytest.approx(1.0, abs=0.01)
