"""Plotting helpers."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from .model import COMPOUNDS


def plot_time_course(result, output: str | Path | None = None):
    fig, ax = plt.subplots()
    for i, name in enumerate(COMPOUNDS):
        ax.plot(result.t, result.y[i], label=name)
    ax.set_xlabel("Time")
    ax.set_ylabel("Concentration")
    ax.legend(frameon=False)
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, bbox_inches="tight")
    return fig, ax
