"""Check pairing, positivity, spatial alignment and all-mode propagation."""
import dataclasses

import numpy as np
import pytest

from geometry.graph_spectrum import circulant_eigenvalue
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import P0
from model.reduced import reduce_S
from scripts.initial_condition_robustness import (SEEDS, initial_fields, linear_delta,
    jacobian_sparsity, protocol, random_fields, spatial_metrics, specifications, translation_comparison)
from simulation.nonlinear_ode import rhs_lattice, rhs_lattice_DNA, pack


def test_protocol_has_unique_paired_designs():
    specs = list(specifications())
    assert len(set(specs)) == len(specs) == protocol()["total_primary_runs"] == 112
    assert sum(s[4] == "free" for s in specs) == 16


def test_initial_states_match_shared_species_and_remain_positive():
    p = dataclasses.replace(P0,lam_N=2)
    fp = homogeneous_fixed_point(p)
    xi = random_fields(SEEDS[0])
    np.testing.assert_allclose(xi.mean(axis=(1,2)),0,atol=1e-16)
    np.testing.assert_allclose(np.max(np.abs(xi),axis=(1,2)),1)
    qss = initial_fields(fp,p,SEEDS[0],.01,"qss")
    free = initial_fields(fp,p,SEEDS[0],.01,"free")
    np.testing.assert_array_equal(qss[:3],free[:3])
    np.testing.assert_allclose(free/fp[:,None,None]-1,.01*xi,atol=2e-16)
    assert qss.min() > 0 and free.min() > 0


def test_translation_alignment_removes_only_periodic_phase():
    a = random_fields(SEEDS[0])[0]+2
    b = np.roll(a,(7,11),axis=(0,1))
    corr, distance, shift = translation_comparison(a,b)
    assert corr == pytest.approx(1)
    assert distance < 1e-15
    np.testing.assert_array_equal(np.roll(b,shift,axis=(0,1)),a)
    assert translation_comparison(a,random_fields(SEEDS[1])[0]+2)[1] > .1


@pytest.mark.parametrize("model",["full","DNA"])
def test_all_mode_propagator_matches_lattice_derivative_and_initial_state(model):
    p = dataclasses.replace(P0,lam_N=2)
    fp = homogeneous_fixed_point(p)
    ns = 4 if model == "full" else 3
    initial = initial_fields(fp,p,SEEDS[0],1e-6)[:ns]
    dt = 1e-5
    pred = linear_delta(initial,fp,p,model,[0,dt])
    np.testing.assert_allclose(pred[0],initial[0]-fp[0],atol=1e-20)
    rhs = rhs_lattice if model == "full" else rhs_lattice_DNA
    actual = rhs(0,pack(*initial),p)[:24*24].reshape(24,24)
    np.testing.assert_allclose((pred[1]-pred[0])/dt,actual,atol=2e-10)


def test_pure_k_eigenmode_propagates_at_its_growth_rate():
    fp = homogeneous_fixed_point(P0)
    J = jacobian_block(-3,fp,P0)
    values,vectors = np.linalg.eig(J)
    k = np.argmax(values.real)
    r,c = np.indices((24,24))
    spatial = np.cos(2*np.pi*(8*r+16*c)/24)
    initial = fp[:,None,None]+1e-6*vectors[:,k].real[:,None,None]*spatial
    pred = linear_delta(initial,fp,P0,"full",[0,2])
    np.testing.assert_allclose(pred[1],pred[0]*np.exp(values[k].real*2),atol=1e-15)
    assert spatial_metrics(initial[0])[1] == pytest.approx(1)


@pytest.mark.parametrize("model",["full","DNA"])
def test_sparse_solver_pattern_contains_all_nonlocal_derivatives(model):
    p=dataclasses.replace(P0,lam_N=2)
    fp=homogeneous_fixed_point(p)
    ns=4 if model=="full" else 3
    y=pack(*initial_fields(fp,p,SEEDS[0],.01)[:ns])
    rhs=rhs_lattice if model=="full" else rhs_lattice_DNA
    sparsity=jacobian_sparsity(model)
    for species in range(ns):
        for cell in (0,1,23,24,200,575):
            column=species*576+cell
            step=1e-5*y[column]
            plus,minus=y.copy(),y.copy()
            plus[column]+=step
            minus[column]-=step
            derivative=(rhs(0,plus,p)-rhs(0,minus,p))/(2*step)
            assert np.all(sparsity[:,column].toarray().ravel()[np.abs(derivative)>1e-7])
