"""The numerical drift and analytic Jacobian must describe the same model."""
import dataclasses

import numpy as np
import pytest

from model.equations import dhillD_dS, dhillS_dA, hill_ASC, hill_D
from model.params import P0


@pytest.mark.parametrize("value", [1e-10, 1e6])
def test_hill_functions_match_unclipped_model_and_derivatives(value):
    p = dataclasses.replace(P0, n_d=2.5, n_ASC=2.5)
    ratio = value / p.Kd
    assert hill_D(np.array(value), p) == pytest.approx(ratio**p.n_d / (1 + ratio**p.n_d), rel=1e-13, abs=0)
    assert hill_ASC(np.array(value), p) == pytest.approx(1 / (1 + ratio**p.n_ASC), rel=1e-13, abs=0)
    step = value * 1e-4
    # Use the unsaturated function for each end of the concentration range.
    f, derivative = (hill_D, dhillD_dS) if value < p.Kd else (hill_ASC, dhillS_dA)
    fd = (f(np.array(value + step), p) - f(np.array(value - step), p)) / (2 * step)
    assert derivative(value, p) == pytest.approx(fd, rel=1e-7, abs=0)


def test_hill_endpoints_and_extreme_arguments_are_finite():
    values = np.array([0, 1e-250, 1e250])
    assert np.array_equal(hill_D(values, P0), [0, 0, 1])
    assert np.array_equal(hill_ASC(values, P0), [1, 1, 0])
