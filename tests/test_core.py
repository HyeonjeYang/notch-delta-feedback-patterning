"""Validation suite.

The first group checks that the revised code reproduces the v1 numbers that
were correct in v1: the baseline equilibrium, the P2/P3 growth rates, the
dispersion relation, the Schur identity, and the critical lam_N. The second
group checks each analytic statement used in the revised manuscript on
random parameter sets.
"""

from __future__ import annotations

import dataclasses

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from control.graph_io import coupling_decomposition, static_io_gain
from control.root_locus import critical_lam_N
from control.schur_verification import verify as schur_verify
from control.state_space import G_exact, nyquist_rhp_count, partition, stage_gains
from geometry.graph_spectrum import graph_eigenvalues, verify_dispersion
from model.fixed_point import _scalar_residual, homogeneous_fixed_point, rhs_single_cell
from model.jacobian import jacobian_block
from model.params import LAM_N_P2, LAM_N_P3, NotchDeltaParams, P0
from model.reduced import reduce_AS, reduce_S, scale_proneural_timescale
from simulation.perturbation import eigenmode_state, k_mode_amplitude, leading_real_eigvec
from stability.mode_scan import DOMINANT_PATTERN_MU, UNIFORM_MU, UNIQUE_MU, lambda_pattern
from stability.sensitivity import PARAM_NAMES

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"

# Bounds of the v1 Latin-hypercube box (rates log-uniform, Hill exponents uniform).
LO = dict(lam_D=0.5, lam_N=0.1, lam_ASC=0.1, Kd=0.03, n_d=1, n_ASC=1, d_D=0.01, d_N=0.01, d_A=0.1,
          d_ASC=0.01, f_D=0.005, f_N=0.06, k_cis=0.035)
HI = dict(lam_D=50, lam_N=10, lam_ASC=10, Kd=3, n_d=6, n_ASC=6, d_D=1, d_N=1, d_A=10, d_ASC=1,
          f_D=0.5, f_N=6, k_cis=3.5)


def random_params(n: int, seed: int = 7, widen: float = 1.0) -> list[NotchDeltaParams]:
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        kw = {}
        for k in LO:
            if k in ("n_d", "n_ASC"):
                kw[k] = rng.uniform(LO[k], HI[k])
            else:
                kw[k] = float(np.exp(rng.uniform(np.log(LO[k] / widen), np.log(HI[k] * widen))))
        out.append(NotchDeltaParams(**kw))
    return out


RANDOM = random_params(60) + random_params(30, seed=11, widen=10.0)


# --- reproduction of correct v1 numbers -------------------------------------

def test_baseline_fixed_point_matches_v1():
    assert homogeneous_fixed_point(P0)[0] == pytest.approx(28.548283450528967, rel=1e-9)


def test_p2_p3_growth_rates_match_v1():
    for lam_N, expected in ((LAM_N_P2, 0.10307), (LAM_N_P3, 0.2163)):
        p = dataclasses.replace(P0, lam_N=lam_N)
        stab = lambda_pattern(homogeneous_fixed_point(p), p)
        assert stab["nonuniform_max_re_lambda"] == pytest.approx(expected, abs=1e-4)
        assert stab["nonuniform_best_mu"] == pytest.approx(DOMINANT_PATTERN_MU)


def test_baseline_is_spatially_stable():
    stab = lambda_pattern(homogeneous_fixed_point(P0), P0)
    assert stab["nonuniform_max_re_lambda"] < 0 and stab["uniform_re_lambda"] < 0


def test_graph_dispersion_matches_adjacency_matrix():
    assert verify_dispersion() < 1e-10
    eigs = graph_eigenvalues()
    assert eigs.max() == pytest.approx(6.0) and eigs.min() == pytest.approx(-3.0)


def test_schur_complement_matches_direct_characteristic_polynomial():
    fp = homogeneous_fixed_point(P0)
    for mu in (UNIFORM_MU, DOMINANT_PATTERN_MU):
        res = schur_verify(jacobian_block(mu, fp, P0))
        assert res["max_coeff_diff"] < 1e-10 and res["max_root_mismatch"] < 1e-10


def test_critical_lam_n_pole_crossing_equals_unity_gain():
    crit = critical_lam_N(DOMINANT_PATTERN_MU, P0)
    assert crit["lam_N_pole"] == pytest.approx(1.304815343, abs=1e-8)
    assert abs(crit["lam_N_pole"] - crit["lam_N_gain"]) < 1e-12
    assert crit["G_at_pole"] == pytest.approx(1.0, abs=1e-10)


# --- statements of the revised analysis -------------------------------------

@pytest.mark.parametrize("p", RANDOM)
def test_unique_positive_equilibrium(p):
    """F(D) is strictly decreasing with F(0) > 0 > F(lam_D/d_D), and the root solves the full system."""
    D = np.linspace(0.0, p.lam_D / p.d_D, 400)
    F = np.array([_scalar_residual(d, p, 6.0) for d in D])
    assert F[0] > 0 > F[-1]
    assert np.all(np.diff(F) < 0)
    fp = homogeneous_fixed_point(p)
    assert np.all(fp > 0)
    assert np.max(np.abs(rhs_single_cell(fp, p))) <= 1e-12 * max(1.0, np.max(fp))


def test_v1_solver_failures_have_equilibria():
    """All 733 global-sample points that v1 labelled 'no_positive_fixed_point'."""
    params = pd.read_csv(DATA_DIR / "global_sample_parameters.csv", float_precision="round_trip")
    failed = params[params.v1_classification == "no_positive_fixed_point"]
    assert len(failed) == 733
    for _, r in failed.iterrows():
        p = NotchDeltaParams(*[float(r[n]) for n in PARAM_NAMES])
        fp = homogeneous_fixed_point(p)
        assert np.all(fp > 0) and np.max(np.abs(rhs_single_cell(fp, p))) <= 1e-12 * max(1.0, np.max(fp))


@pytest.mark.parametrize("p", RANDOM)
def test_open_loop_hurwitz_and_static_identity(p):
    """A_mu is Hurwitz for every graph mode, and det J = det(A) a_S (1 - G)."""
    fp = homogeneous_fixed_point(p)
    for mu in UNIQUE_MU:
        J = jacobian_block(mu, fp, p)
        A_mu, _, _, a_S = partition(J)
        assert np.max(np.linalg.eigvals(A_mu).real) < 0
        lhs, rhs = np.linalg.det(J), np.linalg.det(A_mu) * a_S * (1.0 - G_exact(J))
        assert lhs == pytest.approx(rhs, rel=1e-8, abs=1e-14 * max(1.0, abs(np.linalg.det(A_mu) * a_S)))


@pytest.mark.parametrize("p", RANDOM)
def test_static_gain_is_maximal_at_most_negative_mode(p):
    fp = homogeneous_fixed_point(p)
    G = np.array([G_exact(jacobian_block(mu, fp, p)) for mu in UNIQUE_MU])
    assert UNIQUE_MU[int(np.argmax(G))] == pytest.approx(DOMINANT_PATTERN_MU)
    assert G[-1] < 0  # uniform mode mu = 6: never a stationary instability


@pytest.mark.parametrize("p", RANDOM[:20])
def test_stage_gains_multiply_to_G(p):
    fp = homogeneous_fixed_point(p)
    for mu in (DOMINANT_PATTERN_MU, 0.0, UNIFORM_MU):
        g = stage_gains(jacobian_block(mu, fp, p))
        assert g["G"] == pytest.approx(G_exact(jacobian_block(mu, fp, p)), rel=1e-9, abs=1e-300)


@pytest.mark.parametrize("p", RANDOM[:20])
def test_graph_io_identity(p):
    fp = homogeneous_fixed_point(p)
    J_loc, _, _ = coupling_decomposition(fp, p)
    K = static_io_gain(fp, p)
    for mu in (DOMINANT_PATTERN_MU, 1.5, UNIFORM_MU):
        J = jacobian_block(mu, fp, p)
        assert np.linalg.det(J) == pytest.approx(np.linalg.det(J_loc) * np.linalg.det(np.eye(2) - mu * K),
                                                 rel=1e-7, abs=1e-12 * abs(np.linalg.det(J_loc)))


@pytest.mark.parametrize("p", RANDOM[:20])
def test_reductions_share_the_zero_eigenvalue_condition(p):
    """det J_full = a_S det J_DNA = det(J_AS-block) det J_DN, so the stationary boundary is shared."""
    fp = homogeneous_fixed_point(p)
    J = jacobian_block(DOMINANT_PATTERN_MU, fp, p)
    a_S = J[3, 3]
    assert np.linalg.det(J) == pytest.approx(a_S * np.linalg.det(reduce_S(J)), rel=1e-8, abs=1e-300)
    assert np.linalg.det(J) == pytest.approx(np.linalg.det(J[2:, 2:]) * np.linalg.det(reduce_AS(J)),
                                             rel=1e-8, abs=1e-300)


def test_proneural_timescale_family_keeps_static_quantities():
    p2 = dataclasses.replace(P0, lam_N=LAM_N_P2)
    fp, G = homogeneous_fixed_point(p2), G_exact(jacobian_block(DOMINANT_PATTERN_MU, homogeneous_fixed_point(p2), p2))
    rates = []
    for alpha in (0.1, 1.0, 10.0):
        q = scale_proneural_timescale(p2, alpha)
        fq = homogeneous_fixed_point(q)
        Jq = jacobian_block(DOMINANT_PATTERN_MU, fq, q)
        assert np.allclose(fq, fp, rtol=1e-12) and G_exact(Jq) == pytest.approx(G, rel=1e-12)
        rates.append(np.max(np.linalg.eigvals(Jq).real))
    assert rates[0] < rates[1] < rates[2]


def test_nyquist_count_matches_eigenvalues():
    p2 = dataclasses.replace(P0, lam_N=LAM_N_P2)
    for p, expected in ((P0, 0), (p2, 1)):
        J = jacobian_block(DOMINANT_PATTERN_MU, homogeneous_fixed_point(p), p)
        assert nyquist_rhp_count(J) == expected == int(np.sum(np.linalg.eigvals(J).real > 0))


def test_eigenmode_initial_state_has_exact_direction():
    p2 = dataclasses.replace(P0, lam_N=LAM_N_P2)
    fp = homogeneous_fixed_point(p2)
    _, v = leading_real_eigvec(jacobian_block(DOMINANT_PATTERN_MU, fp, p2))
    fields = eigenmode_state(fp, v, 1e-6)
    amps = np.array([k_mode_amplitude(f) for f in fields])
    assert amps[0] == pytest.approx(1e-6 * fp[0], rel=1e-9)
    assert np.allclose(amps / amps[0], np.abs(v / v[0]), rtol=1e-8)


def test_global_sample_regenerates_from_seed():
    import importlib.util
    spec = importlib.util.spec_from_file_location("global_sample", Path(__file__).resolve().parents[1] / "scripts" / "global_sample.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    stored = pd.read_csv(DATA_DIR / "global_sample_parameters.csv", float_precision="round_trip")
    # log/exp transformations can differ in their last bits across platforms.
    np.testing.assert_allclose(module.latin_hypercube_parameters(), stored[PARAM_NAMES].values,
                               rtol=5e-14, atol=0)


def test_mode_search_regenerates_both_boxes_from_continuous_stream():
    from scripts.mode_dominance_search import N_PER_BOX, SEED, draw
    stored = pd.read_csv(DATA_DIR / "mode_search_samples.csv.gz", float_precision="round_trip")
    assert len(stored) == 2 * N_PER_BOX
    rng = np.random.Generator(np.random.PCG64(SEED))
    generated = []
    for widening in (1.0, 10.0):
        for sample_id in range(N_PER_BOX):
            p = draw(rng, widening)
            generated.append([widening, sample_id, *[getattr(p, k) for k in PARAM_NAMES]])
    np.testing.assert_allclose(stored[["box_widening", "sample_id", *PARAM_NAMES]].values,
                               generated, rtol=5e-14, atol=0)
