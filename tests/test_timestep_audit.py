"""Independent accuracy and grid checks for the fixed-step audit integrator."""
import numpy as np
import pytest

from scripts.timestep_audit import CONFIGS, lattice_cases, protocol, rk4_integrate


def test_rk4_observed_order_against_exact_nonlinear_solution():
    # y'=y^2, y(0)=1 has y(t)=1/(1-t), rather than the linear
    # matrix problem used in the scientific audit.
    errors=[]
    for steps in (8,16,32):
        times,states,_=rk4_integrate(lambda t,y:y*y,[1.],.5/steps,steps)
        assert times[-1]==pytest.approx(.5)
        errors.append(abs(states[-1,0]-2))
    assert 13<errors[0]/errors[1]<19
    assert 13<errors[1]/errors[2]<19


def test_nested_output_grid_preserves_the_same_fixed_steps():
    fun=lambda t,y:np.array([y[1],-y[0]])
    times,fine,_=rk4_integrate(fun,[1.,0.],.01,100)
    coarse_times,coarse,_=rk4_integrate(fun,[1.,0.],.01,100,10)
    np.testing.assert_array_equal(coarse,fine[::10])
    np.testing.assert_allclose(coarse_times,times[::10],atol=1e-15)
    np.testing.assert_allclose(coarse[-1],[np.cos(1),-np.sin(1)],atol=1e-9)


def test_audit_design_contains_all_amplitudes_models_and_downstream_controls():
    cases=list(lattice_cases())
    assert len(set(cases))==len(cases)==14
    assert len(cases)*len(CONFIGS)==protocol()["adaptive_total_runs"]==70
    assert sum(s[-1]=="free" for s in cases)==2


@pytest.mark.parametrize("dt,n,every",[(0,10,1),(-.1,10,1),(.1,0,1),(.1,10,3)])
def test_invalid_fixed_grids_are_rejected(dt,n,every):
    with pytest.raises(ValueError):
        rk4_integrate(lambda t,y:y,[1.],dt,n,every)
