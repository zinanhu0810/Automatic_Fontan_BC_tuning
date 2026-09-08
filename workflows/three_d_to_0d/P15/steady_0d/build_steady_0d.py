#!/usr/bin/env python3
"""Build and run the P15 constant-flow, resistance-outlet 0D model."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import fmean


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE.parent / "rom" / "P15_0d.json"
DEFAULT_MODEL = HERE / "P15_0d_resistance_constant_flow.json"
DEFAULT_RESULT = HERE / "P15_0d_resistance_constant_flow.csv"
DEFAULT_SUMMARY = HERE / "boundary_summary.csv"
MMHG_PER_DYN_CM2 = 1.0 / 1333.22


def trapz_mean(times: list[float], values: list[float]) -> float:
    if len(times) != len(values) or len(times) < 2:
        raise ValueError("FLOW boundary requires at least two paired t/Q values")
    duration = times[-1] - times[0]
    if duration <= 0.0:
        raise ValueError("FLOW boundary time values must increase")
    integral = sum(
        (t1 - t0) * (q0 + q1) / 2.0
        for t0, t1, q0, q1 in zip(times, times[1:], values, values[1:])
    )
    return integral / duration


def build_model(source: Path, output: Path) -> dict:
    model = json.loads(source.read_text())
    renamed_boundaries: dict[str, str] = {}

    for boundary in model["boundary_conditions"]:
        name = boundary["bc_name"]
        values = boundary["bc_values"]
        if boundary["bc_type"] == "FLOW":
            mean_flow = trapz_mean(values["t"], values["Q"])
            boundary["bc_values"] = {
                "t": [0.0, 1.0],
                "Q": [mean_flow, mean_flow],
            }
        elif boundary["bc_type"] == "RCR":
            new_name = name.replace("RCR_", "RESISTANCE_", 1)
            renamed_boundaries[name] = new_name
            boundary["bc_name"] = new_name
            boundary["bc_type"] = "RESISTANCE"
            boundary["bc_values"] = {
                "R": float(values["Rp"]) + float(values["Rd"]),
                "Pd": float(values.get("Pd", 0.0)),
            }
        else:
            raise ValueError(f"Unsupported boundary type: {boundary['bc_type']}")

    for vessel in model["vessels"]:
        for side, name in vessel.get("boundary_conditions", {}).items():
            vessel["boundary_conditions"][side] = renamed_boundaries.get(name, name)

    parameters = model["simulation_parameters"]
    parameters["model_name"] = "P15_resistance_constant_flow"
    # Retain the baseline duration and temporal resolution to let vessel RCL
    # transients decay before the solver exports its final cycle.
    parameters["number_of_cardiac_cycles"] = 18
    parameters["number_of_time_pts_per_cardiac_cycle"] = 1783

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(model, indent=2) + "\n")
    return model


def find_solver(explicit: Path | None) -> Path:
    candidates = [
        explicit,
        Path(shutil.which("svzerodsolver")) if shutil.which("svzerodsolver") else None,
        Path.home() / "anaconda3" / "bin" / "svzerodsolver",
    ]
    solver = next((path.expanduser() for path in candidates if path and path.is_file()), None)
    if solver is None:
        raise FileNotFoundError("svzerodsolver was not found; pass --solver PATH")
    return solver.resolve()


def boundary_attachments(model: dict) -> dict[str, tuple[str, str]]:
    attachments = {}
    for vessel in model["vessels"]:
        for side, boundary_name in vessel.get("boundary_conditions", {}).items():
            attachments[boundary_name] = (vessel["vessel_name"], side)
    return attachments


def summarize(model: dict, result: Path, summary: Path) -> dict:
    with result.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("svZeroDSolver produced an empty result")

    numeric = ("time", "flow_in", "flow_out", "pressure_in", "pressure_out")
    if not all(math.isfinite(float(row[key])) for row in rows for key in numeric):
        raise ValueError("svZeroDSolver output contains non-finite values")

    by_vessel: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_vessel[row["name"]].append(row)

    attachments = boundary_attachments(model)
    boundaries = {bc["bc_name"]: bc for bc in model["boundary_conditions"]}
    summary_rows = []
    root_flow = None
    terminal_sum = 0.0
    resistance_outflow = 0.0

    for name, boundary in boundaries.items():
        vessel_name, side = attachments[name]
        vessel_rows = by_vessel[vessel_name]
        flow_key = "flow_in" if side == "inlet" else "flow_out"
        pressure_key = "pressure_in" if side == "inlet" else "pressure_out"
        flows = [float(row[flow_key]) for row in vessel_rows]
        pressures = [float(row[pressure_key]) * MMHG_PER_DYN_CM2 for row in vessel_rows]
        mean_flow = fmean(flows)

        if side == "inlet":
            if root_flow is not None:
                raise ValueError("Expected exactly one rooted inlet boundary")
            root_flow = mean_flow
        else:
            terminal_sum += mean_flow
        if boundary["bc_type"] == "RESISTANCE":
            resistance_outflow += mean_flow

        prescribed = boundary["bc_values"].get("R")
        if boundary["bc_type"] == "FLOW":
            prescribed = boundary["bc_values"]["Q"][0]
        summary_rows.append(
            {
                "boundary": name,
                "type": boundary["bc_type"],
                "side": side,
                "prescribed_Q_or_R": prescribed,
                "mean_flow_ml_s": mean_flow,
                "flow_range_ml_s": max(flows) - min(flows),
                "mean_pressure_mmHg": fmean(pressures),
                "pressure_range_mmHg": max(pressures) - min(pressures),
            }
        )

    if root_flow is None:
        raise ValueError("No rooted inlet boundary was found")
    residual = root_flow - terminal_sum
    scale = abs(root_flow) + abs(terminal_sum)
    relative_residual = abs(residual) / scale if scale else 0.0

    with summary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)

    return {
        "number_of_rows": len(rows),
        "number_of_boundaries": len(boundaries),
        "number_of_flow_boundaries": sum(
            bc["bc_type"] == "FLOW" for bc in boundaries.values()
        ),
        "number_of_resistance_boundaries": sum(
            bc["bc_type"] == "RESISTANCE" for bc in boundaries.values()
        ),
        "mean_root_flow_ml_s": root_flow,
        "mean_terminal_flow_sum_ml_s": terminal_sum,
        "mean_resistance_outflow_ml_s": resistance_outflow,
        "mass_residual_ml_s": residual,
        "relative_mass_residual": relative_residual,
        "max_boundary_flow_range_ml_s": max(
            row["flow_range_ml_s"] for row in summary_rows
        ),
        "max_boundary_pressure_range_mmHg": max(
            row["pressure_range_mmHg"] for row in summary_rows
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--build-only", action="store_true")
    args = parser.parse_args()

    model = build_model(args.source.resolve(), args.model.resolve())
    print(f"Wrote model: {args.model.resolve()}")
    if args.build_only:
        return

    solver = find_solver(args.solver)
    completed = subprocess.run(
        [str(solver), str(args.model.resolve()), str(args.result.resolve())],
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"svZeroDSolver failed with exit code {completed.returncode}")

    metrics = summarize(model, args.result.resolve(), args.summary.resolve())
    metrics_path = args.summary.resolve().with_name("simulation_summary.json")
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"Wrote result: {args.result.resolve()}")
    print(f"Wrote boundary summary: {args.summary.resolve()}")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
