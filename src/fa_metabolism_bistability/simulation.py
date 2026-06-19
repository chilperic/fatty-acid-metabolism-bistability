"""Simulation utilities."""

from __future__ import annotations

from collections.abc import Callable, Mapping

import numpy as np
from scipy.integrate import solve_ivp


def simulate_time_course(
    rhs: Callable,
    y0: np.ndarray,
    params: Mapping[str, float],
    t_span: tuple[float, float] = (0.0, 50.0),
    n_points: int = 500,
    method: str = "LSODA",
    rtol: float = 1e-8,
    atol: float = 1e-10,
):
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    return solve_ivp(lambda t, y: rhs(t, y, params), t_span=t_span, y0=np.asarray(y0, dtype=float), t_eval=t_eval, method=method, rtol=rtol, atol=atol)


def approximate_steady_state(rhs: Callable, y0: np.ndarray, params: Mapping[str, float], horizon: float = 1_000.0) -> np.ndarray:
    result = simulate_time_course(rhs, y0, params, t_span=(0.0, horizon), n_points=2)
    if not result.success:
        raise RuntimeError(result.message)
    return result.y[:, -1]
