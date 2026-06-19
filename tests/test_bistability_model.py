import numpy as np
from fa_metabolism_bistability.model import default_parameters, rates, rhs
from fa_metabolism_bistability.rate_laws import michaelis_menten


def test_michaelis_menten_positive():
    assert michaelis_menten(2.0, 1.0, 1.0) > 0


def test_rhs_shape():
    params = default_parameters()
    y = np.array([1.0, 0.5, 0.2, 0.1])
    dy = rhs(0.0, y, params)
    assert dy.shape == (4,)


def test_rates_include_expected_keys():
    params = default_parameters()
    y = np.ones(4)
    r = rates(y, params)
    assert {"r1", "r2", "r3", "r4", "r5"}.issubset(r)
