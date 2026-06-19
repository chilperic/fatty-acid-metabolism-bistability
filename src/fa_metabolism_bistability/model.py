"""Coarse-grained ODE model of fatty-acid metabolism.

State vector order:
    [acetyl_coa, malonyl_coa, fatty_acids, triglycerides]
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from .rate_laws import constant_influx, michaelis_menten, non_competitive_inhibition, proportional_outflux

COMPOUNDS = ("acetyl_coa", "malonyl_coa", "fatty_acids", "triglycerides")


def default_parameters() -> dict[str, float]:
    return {
        "k1": 1.0,      # acetyl-CoA influx
        "k2": 1.0,      # fatty-acid influx
        "k3": 1.0,      # triglyceride influx
        "alpha": 1.0,   # acetyl-CoA outflux coefficient
        "beta": 1.0,    # fatty-acid outflux coefficient
        "gamma": 1.0,   # triglyceride outflux coefficient
        "V1": 1.0,
        "V2": 1.0,
        "V3": 1.0,
        "V4": 1.0,
        "V5": 1.0,
        "Km1": 1.0,
        "Km2": 1.0,
        "Km3": 1.0,
        "Km4": 1.0,
        "Km5": 1.0,
        "q1": 1.0,
        "q4": 1.0,
    }


def _p(params: Mapping[str, float], key: str) -> float:
    return float(params[key])


def rates(y: np.ndarray, params: Mapping[str, float]) -> dict[str, float]:
    acetyl_coa, malonyl_coa, fatty_acids, triglycerides = np.asarray(y, dtype=float)
    return {
        "v_in_1": constant_influx(_p(params, "k1")),
        "v_in_2": constant_influx(_p(params, "k2")),
        "v_in_3": constant_influx(_p(params, "k3")),
        "r1": non_competitive_inhibition(_p(params, "V1"), _p(params, "Km1"), _p(params, "q1"), acetyl_coa, fatty_acids),
        "r2": michaelis_menten(_p(params, "V2"), _p(params, "Km2"), malonyl_coa),
        "r3": michaelis_menten(_p(params, "V3"), _p(params, "Km3"), fatty_acids),
        "r4": non_competitive_inhibition(_p(params, "V4"), _p(params, "Km4"), _p(params, "q4"), fatty_acids, malonyl_coa),
        "r5": michaelis_menten(_p(params, "V5"), _p(params, "Km5"), triglycerides),
        "v_out1": proportional_outflux(_p(params, "alpha"), acetyl_coa),
        "v_out2": proportional_outflux(_p(params, "beta"), fatty_acids),
        "v_out3": proportional_outflux(_p(params, "gamma"), triglycerides),
    }


def rhs(t: float, y: np.ndarray, params: Mapping[str, float]) -> np.ndarray:
    r = rates(y, params)
    d_acetyl = r["v_in_1"] - r["r1"] + r["r4"] - r["v_out1"]
    d_malonyl = r["r1"] - r["r2"]
    d_fatty_acids = r["v_in_2"] + r["r2"] - r["r3"] - r["r4"] + r["r5"] - r["v_out2"]
    d_triglycerides = r["v_in_3"] + r["r3"] - r["r5"] - r["v_out3"]
    return np.array([d_acetyl, d_malonyl, d_fatty_acids, d_triglycerides], dtype=float)
