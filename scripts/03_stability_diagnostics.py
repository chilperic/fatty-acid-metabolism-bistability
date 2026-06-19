import numpy as np
from fa_metabolism_bistability.model import default_parameters, rhs
from fa_metabolism_bistability.simulation import approximate_steady_state
from fa_metabolism_bistability.analysis import eigenvalues_at_state, classify_stability

params = default_parameters()
y0 = np.array([1.0, 0.0, 0.0, 0.0])
ss = approximate_steady_state(rhs, y0, params, horizon=200)
eigs = eigenvalues_at_state(rhs, ss, params)
print("steady_state:", ss)
print("eigenvalues:", eigs)
print("classification:", classify_stability(eigs))
