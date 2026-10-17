"""Independent checks for the extra biological and geometry controls."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.linalg import block_diag

from geometry.lattice import hex_adj
from geometry.graph_spectrum import circulant_eigenvalue
from model.params import P0
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from scripts.reviewer_models import adjacency,lattice_rhs,lattice_jacobian,production_lag_block,lag_similarity,random_initial


def test_variable_graph_matches_original_contacts_and_open_graph_degree():
    x=np.arange(24*24).reshape(24,24)
    np.testing.assert_array_equal((adjacency(x.shape)@x.ravel()).reshape(x.shape),hex_adj(x))
    W=adjacency((7,8),False)
    np.testing.assert_array_equal(W.toarray(),W.T.toarray())
    assert np.min(W.sum(axis=1))==2
    assert np.max(W.sum(axis=1))==6


@pytest.mark.parametrize("model",["full","DNA"])
@pytest.mark.parametrize("periodic",[True,False])
def test_patterned_analytic_jacobian_matches_directional_derivatives(model,periodic):
    p=replace(P0,lam_N=2)
    fp=homogeneous_fixed_point(p)
    ns=4 if model=="full" else 3
    fields=random_initial(fp,p,(5,6),123,.1)[:ns]
    y=fields.ravel()
    W=adjacency((5,6),periodic)
    v=np.random.default_rng(4242).normal(size=len(y))*y
    J=lattice_jacobian(y,p,W,model)
    h=1e-5
    fd=(lattice_rhs(0,y+h*v,p,W,model)-lattice_rhs(0,y-h*v,p,W,model))/(2*h)
    np.testing.assert_allclose(J@v,fd,rtol=2e-7,atol=1e-8)


def test_full_sparse_jacobian_has_the_expected_plane_wave_blocks():
    p=replace(P0,lam_N=2)
    fp=homogeneous_fixed_point(p)
    shape=(6,9);r,c=np.indices(shape)
    wave=np.exp(2j*np.pi*(2*r/6+4*c/9)).ravel()
    q=np.array([1.,-.3,.5,.2])
    y=np.broadcast_to(fp[:,None],(4,54)).ravel()
    actual=lattice_jacobian(y,p,adjacency(shape))@np.outer(q,wave).ravel()
    expected=np.outer(jacobian_block(circulant_eigenvalue(2,4,6,9),fp,p)@q,wave).ravel()
    np.testing.assert_allclose(actual,expected,atol=1e-13)


def test_production_lag_is_a_linear_similarity_for_every_mode():
    for lam in (1.,2.,5.623):
        p=replace(P0,lam_N=lam)
        fp=homogeneous_fixed_point(p)
        T=lag_similarity(fp,p)
        for mu in (-3.,-2.,0.,6.):
            J=jacobian_block(mu,fp,p)
            np.testing.assert_allclose(production_lag_block(mu,fp,p),T@J@np.linalg.inv(T),atol=2e-13)
