import numpy as np
from fa_metabolism_bistability.model import default_parameters, rhs
from fa_metabolism_bistability.simulation import approximate_steady_state

base = default_parameters()
y0 = np.array([1.0, 0.0, 0.0, 0.0])
for value in np.linspace(0.1, 3.0, 10):
    params = dict(base)
    params["k1"] = float(value)
    ss = approximate_steady_state(rhs, y0, params, horizon=200)
    print(f"k1={value:.3f} steady_state={ss}")
