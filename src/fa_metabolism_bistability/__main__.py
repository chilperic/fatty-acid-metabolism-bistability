"""Run an explicit JSON configuration and export results with provenance."""

import argparse
import json
from pathlib import Path

import numpy as np

from .provenance import write_manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if args.output.exists():
        parser.error("Output path already exists; choose a new run directory")
    from .model import COMPOUNDS, default_parameters, rhs
    from .simulation import simulate_time_course

    params = default_parameters()
    if set(config["parameters"]) - set(params):
        raise ValueError("Unknown model parameters")
    params.update(config["parameters"])
    config["parameters"] = params
    result = simulate_time_course(
        rhs, config["initial_state"], params, **config["solver"]
    )
    compounds = COMPOUNDS
    args.output.mkdir(parents=True)
    np.savetxt(
        args.output / "trajectory.csv",
        np.column_stack([result.t, result.y.T]),
        delimiter=",",
        header=",".join(["time", *compounds]),
        comments="",
    )
    (args.output / "solver.json").write_text(
        json.dumps(
            {
                "success": bool(result.success),
                "message": result.message,
                "nfev": int(result.nfev),
                "njev": int(result.njev),
            },
            indent=2,
        )
        + "\n"
    )
    write_manifest(
        args.output,
        config,
        args.config,
        "fa-metabolism-bistability",
        ["numpy", "scipy"],
    )


if __name__ == "__main__":
    main()
