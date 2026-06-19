"""Reusable kinetic rate laws used by the coarse-grained model."""

from __future__ import annotations


def constant_influx(k: float) -> float:
    return float(k)


def proportional_outflux(k: float, substrate: float) -> float:
    return float(k) * float(substrate)


def mass_action(k: float, *substrates: float) -> float:
    rate = float(k)
    for s in substrates:
        rate *= float(s)
    return rate


def michaelis_menten(vmax: float, km: float, substrate: float) -> float:
    substrate = float(substrate)
    return float(vmax) * substrate / (float(km) + substrate)


def non_competitive_inhibition(vmax: float, km: float, q: float, substrate: float, inhibitor: float) -> float:
    substrate = float(substrate)
    inhibitor = float(inhibitor)
    return float(vmax) * substrate / ((float(km) + substrate) * (1.0 + float(q) * inhibitor))
