#!/usr/bin/env python3
"""Tune pulmonary terminal resistances in an svZeroD model to clinical targets."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import tempfile
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from statistics import fmean


MMHG_PER_DYN_CM2 = 1.0 / 1333.22
# The 0D BC tuner is inexpensive, so use a tighter target than the 3D
# acceptance criterion before passing the tuned BCs back to CFD.
TOLERANCE = 0.001
MIN_RESISTANCE_SCALE = 0.001
MAX_RESISTANCE_SCALE = 100.0


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--model", type=Path, required=True, help="Calibrated 0D JSON model.")
    parser.add_argument("--output-directory", default="tuned_bc")
    parser.add_argument("--svzerodsolver", type=Path, default=None)
    parser.add_argument("--max-evaluations", type=int, default=80)
    return parser.parse_args()


def solver_path(explicit: Path | None) -> Path:
    candidates = [explicit, Path(shutil.which("svzerodsolver")) if shutil.which("svzerodsolver") else None, Path.home() / "anaconda3/bin/svzerodsolver"]
    solver = next((path.resolve() for path in candidates if path and path.is_file()), None)
    if solver is None:
        raise FileNotFoundError("svzerodsolver was not found; pass --svzerodsolver PATH.")
    return solver


def attachments(model: dict) -> dict[str, tuple[str, str]]:
    return {name: (vessel["vessel_name"], side) for vessel in model["vessels"] for side, name in vessel.get("boundary_conditions", {}).items()}


def terminal_boundary_name(boundary_names: set[str], cap: str) -> str:
    for prefix in ("RCR_", "RESISTANCE_"):
        name = prefix + cap
        if name in boundary_names:
            return name
    raise KeyError(f"No RCR or RESISTANCE boundary is defined for outlet '{cap}'.")


def scale_terminal_resistances(original: dict, lpa: list[str], rpa: list[str], scales: tuple[float, float]) -> dict:
    common, bias = scales
    model = deepcopy(original)
    source_vessels = original["vessels"]
    names = {item["bc_name"] for item in model["boundary_conditions"]}
    lpa_names = {terminal_boundary_name(names, name) for name in lpa}
    rpa_names = {terminal_boundary_name(names, name) for name in rpa}
    seen: set[str] = set()
    for boundary in model["boundary_conditions"]:
        name = boundary["bc_name"]
        if name in lpa_names:
            factor = common * bias
        elif name in rpa_names:
            factor = common / bias
        else:
            continue
        if boundary["bc_type"] == "RCR":
            boundary["bc_values"]["Rp"] *= factor
            boundary["bc_values"]["Rd"] *= factor
        elif boundary["bc_type"] == "RESISTANCE":
            boundary["bc_values"]["R"] *= factor
        else:
            raise ValueError(f"Unsupported terminal boundary type: {boundary['bc_type']}")
        seen.add(name)
    if seen != lpa_names | rpa_names:
        raise ValueError("Every configured LPA/RPA outlet must have an RCR or RESISTANCE boundary in the 0D model.")
    if model["vessels"] != source_vessels:
        raise RuntimeError("Vessel parameters changed during RCR-only tuning.")
    return model


def run_solver(solver: Path, model: dict, folder: Path) -> Path:
    model_path, csv_path = folder / "model.json", folder / "result.csv"
    model_path.write_text(json.dumps(model, indent=2) + "\n")
    completed = subprocess.run([str(solver), str(model_path), str(csv_path)], capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip() or "svzerodsolver failed")
    return csv_path


def metrics(model: dict, result: Path, config: dict) -> dict[str, float]:
    by_vessel: dict[str, list[dict[str, str]]] = defaultdict(list)
    with result.open() as stream:
        for row in csv.DictReader(stream):
            by_vessel[row["name"]].append(row)
    count = int(config["simulation"]["steps_per_cycle"])
    mapping = attachments(model)

    def boundary_values(name: str) -> tuple[list[float], list[float]]:
        vessel, side = mapping[name]
        rows = by_vessel[vessel][-count:]
        flow = "flow_in" if side == "inlet" else "flow_out"
        pressure = "pressure_in" if side == "inlet" else "pressure_out"
        return [float(row[flow]) for row in rows], [float(row[pressure]) * MMHG_PER_DYN_CM2 for row in rows]

    names = set(mapping)
    lpa = [boundary_values(terminal_boundary_name(names, name))[0] for name in config["outlets"]["lpa"]]
    rpa = [boundary_values(terminal_boundary_name(names, name))[0] for name in config["outlets"]["rpa"]]
    _, fontan = boundary_values(f"FLOW_{config['metrics']['fontan_pressure_face']}")
    mean_lpa = fmean(sum(values[index] for values in lpa) for index in range(count))
    mean_rpa = fmean(sum(values[index] for values in rpa) for index in range(count))
    return {"mean_fontan_pressure_mmhg": fmean(fontan), "mean_lpa_flow_ml_s": mean_lpa, "mean_rpa_flow_ml_s": mean_rpa, "lpa_fraction": mean_lpa / (mean_lpa + mean_rpa)}


def objective(values: dict[str, float], target: dict[str, float]) -> float:
    errors = [abs(values["mean_fontan_pressure_mmhg"] / target["mean_fontan_pressure_mmhg"] - 1.0), abs(values["lpa_fraction"] / target["lpa_fraction"] - 1.0)]
    return max(errors) + 0.01 * sum(error**2 for error in errors)


def main() -> None:
    args = arguments()
    case = args.case_directory.resolve()
    config = json.loads((case / "step0_config.json").read_text())
    model_path = args.model if args.model.is_absolute() else case / args.model
    original = json.loads(model_path.read_text())
    solver = solver_path(args.svzerodsolver)
    output = Path(args.output_directory)
    output = output if output.is_absolute() else case / output
    output.mkdir(parents=True, exist_ok=True)
    target = {"mean_fontan_pressure_mmhg": float(config["targets"]["mean_fontan_pressure_mmhg"]), "lpa_fraction": float(config["targets"]["lpa_fraction_percent"]) / 100.0}
    cache: dict[tuple[float, float], tuple[dict[str, float], float]] = {}
    history: list[dict[str, float]] = []

    class EvaluationLimit(RuntimeError):
        pass

    def evaluate(position: tuple[float, float]) -> tuple[dict[str, float], float]:
        scales = tuple(math.exp(min(math.log(MAX_RESISTANCE_SCALE), max(math.log(MIN_RESISTANCE_SCALE), value))) for value in position)
        key = tuple(round(value, 10) for value in scales)
        if key not in cache:
            if len(history) >= args.max_evaluations:
                raise EvaluationLimit
            candidate = scale_terminal_resistances(original, config["outlets"]["lpa"], config["outlets"]["rpa"], scales)
            with tempfile.TemporaryDirectory(prefix="tune_rcr_") as folder:
                values = metrics(candidate, run_solver(solver, candidate, Path(folder)), config)
            score = objective(values, target)
            cache[key] = (values, score)
            history.append({"evaluation": len(history) + 1, "common_rcr_scale": scales[0], "lpa_rpa_bias": scales[1], "objective": score, **values})
        return cache[key]

    simplex = [(0.0, 0.0), (math.log(1.5), 0.0), (0.0, math.log(1.5))]
    for point in simplex:
        evaluate(point)
    try:
        while len(history) < args.max_evaluations:
            simplex.sort(key=lambda point: evaluate(point)[1])
            best, second, worst = simplex
            centroid = ((best[0] + second[0]) / 2.0, (best[1] + second[1]) / 2.0)
            reflected = (2.0 * centroid[0] - worst[0], 2.0 * centroid[1] - worst[1])
            reflected_value = evaluate(reflected)[1]
            if reflected_value < evaluate(best)[1]:
                expanded = (centroid[0] + 2.0 * (reflected[0] - centroid[0]), centroid[1] + 2.0 * (reflected[1] - centroid[1]))
                simplex[-1] = expanded if evaluate(expanded)[1] < reflected_value else reflected
            elif reflected_value < evaluate(second)[1]:
                simplex[-1] = reflected
            else:
                contracted = (centroid[0] + 0.5 * (worst[0] - centroid[0]), centroid[1] + 0.5 * (worst[1] - centroid[1]))
                simplex[-1] = contracted if evaluate(contracted)[1] < evaluate(worst)[1] else ((best[0] + worst[0]) / 2.0, (best[1] + worst[1]) / 2.0)
            best_values = evaluate(min(simplex, key=lambda point: evaluate(point)[1]))[0]
            if max(abs(best_values["mean_fontan_pressure_mmhg"] / target["mean_fontan_pressure_mmhg"] - 1.0), abs(best_values["lpa_fraction"] / target["lpa_fraction"] - 1.0)) <= TOLERANCE:
                break
    except EvaluationLimit:
        pass

    best = min(simplex, key=lambda point: evaluate(point)[1])
    scales = tuple(math.exp(min(math.log(MAX_RESISTANCE_SCALE), max(math.log(MIN_RESISTANCE_SCALE), value))) for value in best)
    tuned = scale_terminal_resistances(original, config["outlets"]["lpa"], config["outlets"]["rpa"], scales)
    with tempfile.TemporaryDirectory(prefix="tune_rcr_final_") as folder:
        final_csv = run_solver(solver, tuned, Path(folder))
        values = metrics(tuned, final_csv, config)
        (output / f"{config['case_name']}_0d_tuned_bc.csv").write_text(final_csv.read_text())
    (output / f"{config['case_name']}_0d_tuned_bc.json").write_text(json.dumps(tuned, indent=2) + "\n")
    with (output / "tuning_history.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(history[0])); writer.writeheader(); writer.writerows(history)
    lines = ["# Clinical Terminal-Resistance BC Tuning", "", "| Metric | Clinical target | Tuned 0D | Relative error | Status |", "| --- | ---: | ---: | ---: | --- |"]
    for label, key, target_value, unit, scale in (("Mean Fontan pressure", "mean_fontan_pressure_mmhg", target["mean_fontan_pressure_mmhg"], "mmHg", 1.0), ("LPA flow split", "lpa_fraction", target["lpa_fraction"], "%", 100.0)):
        error = abs(values[key] / target_value - 1.0)
        lines.append(f"| {label} | {target_value * scale:.4f} {unit} | {values[key] * scale:.4f} {unit} | {error * 100.0:.3f}% | {'PASS' if error <= TOLERANCE else 'FAIL'} |")
    lines.extend(["", "## Terminal Resistance Scales", "", f"- Common terminal-resistance scale: `{scales[0]:.8g}`", f"- LPA/RPA bias: `{scales[1]:.8g}`", f"- LPA terminal-resistance scale: `{scales[0] * scales[1]:.8g}`", f"- RPA terminal-resistance scale: `{scales[0] / scales[1]:.8g}`", "- Distal pressures, prescribed inflows, and geometric vessel parameters are unchanged.", f"- Evaluations: `{len(history)}`"])
    (output / "tuning_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
