"""Coarse-grained fatty-acid metabolism bi-stability model."""

from .model import COMPOUNDS, default_parameters, rates, rhs
from .simulation import simulate_time_course, approximate_steady_state
from .analysis import finite_difference_jacobian, eigenvalues_at_state, classify_stability

__all__ = [
    "COMPOUNDS",
    "default_parameters",
    "rates",
    "rhs",
    "simulate_time_course",
    "approximate_steady_state",
    "finite_difference_jacobian",
    "eigenvalues_at_state",
    "classify_stability",
]
