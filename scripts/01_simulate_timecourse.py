import numpy as np
from fa_metabolism_bistability.model import default_parameters, rhs
from fa_metabolism_bistability.simulation import simulate_time_course
from fa_metabolism_bistability.plotting import plot_time_course

params = default_parameters()
y0 = np.array([1.0, 0.0, 0.0, 0.0])
result = simulate_time_course(rhs, y0, params, t_span=(0, 50), n_points=300)
print("success:", result.success)
print("final_state:", result.y[:, -1])
plot_time_course(result, output="timecourse_bistability_example.png")
