"""Local stability and parameter-scan utilities."""

from __future__ import annotations

from collections.abc import Callable, Mapping

import numpy as np


def finite_difference_jacobian(
    rhs: Callable, y: np.ndarray, params: Mapping[str, float], eps: float = 1e-6
) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if (
        y.ndim != 1
        or not y.size
        or not np.all(np.isfinite(y))
        or not np.isfinite(eps)
        or eps <= 0
    ):
        raise ValueError("finite vector and positive finite eps required")
    n = y.size
    jac = np.zeros((n, n), dtype=float)
    for i in range(n):
        step = np.zeros(n)
        step[i] = eps * max(1.0, abs(y[i]))
        f_plus = rhs(0.0, y + step, params)
        f_minus = rhs(0.0, y - step, params)
        jac[:, i] = (f_plus - f_minus) / (2.0 * step[i])
    return jac


def eigenvalues_at_state(
    rhs: Callable, y: np.ndarray, params: Mapping[str, float]
) -> np.ndarray:
    return np.linalg.eigvals(finite_difference_jacobian(rhs, y, params))


def classify_stability(eigenvalues: np.ndarray, tol: float = 1e-8) -> str:
    eigenvalues = np.asarray(eigenvalues)
    if (
        eigenvalues.ndim != 1
        or not eigenvalues.size
        or not np.all(np.isfinite(eigenvalues))
        or not np.isfinite(tol)
        or tol < 0
    ):
        raise ValueError("finite nonempty eigenvalues and nonnegative tol required")
    real_parts = np.real(eigenvalues)
    if np.all(real_parts < -tol):
        return "locally_stable"
    if np.any(real_parts > tol):
        return "unstable"
    return "non_hyperbolic_or_marginal"


def scan_parameter(
    base_params: Mapping[str, float], name: str, values, steady_state_fn: Callable
) -> list[dict[str, object]]:
    rows = []
    for value in values:
        params = dict(base_params)
        params[name] = float(value)
        state = steady_state_fn(params)
        rows.append({"parameter": name, "value": float(value), "steady_state": state})
    return rows
