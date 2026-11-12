import numpy as np
import pytest
from fa_metabolism_bistability.analysis import (
    classify_stability,
    finite_difference_jacobian,
)
from fa_metabolism_bistability.model import default_parameters, rhs
from fa_metabolism_bistability.simulation import (
    approximate_steady_state,
)


def test_total_pool_balance():
    rng = np.random.default_rng(21)
    p = default_parameters()
    for y in rng.uniform(0.1, 10, (20, 4)):
        expected = (
            p["k1"]
            + p["k2"]
            + p["k3"]
            - p["alpha"] * y[0]
            - p["beta"] * y[2]
            - p["gamma"] * y[3]
        )
        assert np.sum(rhs(0, y, p)) == pytest.approx(expected)


def test_no_false_steady_state_for_drift_or_short_horizon():
    with pytest.raises(RuntimeError, match="not reached"):
        approximate_steady_state(lambda t, y, p: np.ones_like(y), [0.0], {}, horizon=10)
    with pytest.raises(RuntimeError, match="not reached"):
        approximate_steady_state(lambda t, y, p: -y, [1.0], {}, horizon=0.1)
    np.testing.assert_allclose(
        approximate_steady_state(lambda t, y, p: -y, [1.0], {}), 0, atol=1e-8
    )


def test_default_state_residual_and_cross_solver():
    p = default_parameters()
    a = approximate_steady_state(rhs, np.ones(4), p)
    b = approximate_steady_state(rhs, np.ones(4), p, method="Radau")
    assert np.linalg.norm(rhs(0, a, p), np.inf) <= 1e-8
    np.testing.assert_allclose(a, b, atol=1e-7)


def test_jacobian_on_known_linear_system():
    matrix = np.array([[-1.0, 2.0], [0.0, -3.0]])
    np.testing.assert_allclose(
        finite_difference_jacobian(lambda t, y, p: matrix @ y, [1, 2], {}),
        matrix,
        atol=1e-8,
    )
    assert classify_stability(np.linalg.eigvals(matrix)) == "locally_stable"
    with pytest.raises(ValueError):
        classify_stability([np.nan])
