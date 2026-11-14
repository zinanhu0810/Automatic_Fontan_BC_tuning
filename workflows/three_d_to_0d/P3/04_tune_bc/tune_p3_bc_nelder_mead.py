#!/usr/bin/env python3
"""Tune P3 pulmonary RCR BCs to clinical targets with bounded Nelder-Mead."""

from __future__ import annotations

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


HERE = Path(__file__).resolve().parent
PIPELINE = HERE.parent
SOURCE_MODEL = PIPELINE / "03_calibrate_0d" / "unsteady" / "P3_0d_calibrated.json"
TARGET_FILE = HERE / "clinical_targets.json"
OUTPUT = HERE / "unsteady"
MMHG_PER_DYN_CM2 = 1.0 / 1333.22
STEPS_PER_CYCLE = 1000
LPA_FACES = ("cap_LPA_1", "cap_LPA_2", "cap_LPA_3", "cap_LPA_4", "cap_LPA_5", "cap_LPA_6", "cap_LPA_7", "cap_LPA_8", "cap_PA_main_2")
RPA_FACES = ("cap_PA_main", "cap_RPA_1", "cap_RPA_2", "cap_RPA_3", "cap_RPA_4", "cap_RPA_5", "cap_RPA_6", "cap_RPA_7", "cap_RPA_8")
FONTAN_FACE = "cap_Fontan"
LOG_BOUNDS = ((math.log(0.1), math.log(10.0)), (math.log(0.2), math.log(5.0)))


def find_solver() -> Path:
    candidates = [
        Path(shutil.which("svzerodsolver")) if shutil.which("svzerodsolver") else None,
        Path.home() / "anaconda3" / "bin" / "svzerodsolver",
    ]
    solver = next((candidate for candidate in candidates if candidate and candidate.is_file()), None)
    if solver is None:
        raise FileNotFoundError("svzerodsolver was not found")
    return solver.resolve()


def attachments(model: dict) -> dict[str, tuple[str, str]]:
    return {
        boundary: (vessel["vessel_name"], side)
        for vessel in model["vessels"]
        for side, boundary in vessel.get("boundary_conditions", {}).items()
    }


def clamp(position: tuple[float, float]) -> tuple[float, float]:
    return tuple(
        min(max(value, bounds[0]), bounds[1])
        for value, bounds in zip(position, LOG_BOUNDS)
    )


def apply_rcr_scales(source: dict, position: tuple[float, float]) -> dict:
    """Scale paired Rp/Rd values only; preserve C, Pd, FLOW BCs, and vessels."""
    common_log, lpa_bias_log = position
    common = math.exp(common_log)
    lpa_scale = common * math.exp(lpa_bias_log)
    rpa_scale = common * math.exp(-lpa_bias_log)
    model = deepcopy(source)
    source_vessels = source["vessels"]
    source_flow = [bc for bc in source["boundary_conditions"] if bc["bc_type"] == "FLOW"]
    found = set()
    for bc in model["boundary_conditions"]:
        face = bc["bc_name"].removeprefix("RCR_")
        if face in LPA_FACES:
            scale = lpa_scale
        elif face in RPA_FACES:
            scale = rpa_scale
        else:
            continue
        bc["bc_values"]["Rp"] *= scale
        bc["bc_values"]["Rd"] *= scale
        found.add(bc["bc_name"])
    expected = {f"RCR_{face}" for face in LPA_FACES + RPA_FACES}
    if found != expected:
        raise ValueError(f"Unexpected RCR boundaries: {found}")
    if model["vessels"] != source_vessels:
        raise RuntimeError("Geometric vessel parameters changed during BC tuning")
    if [bc for bc in model["boundary_conditions"] if bc["bc_type"] == "FLOW"] != source_flow:
        raise RuntimeError("Prescribed inflow changed during BC tuning")
    model["simulation_parameters"]["model_name"] = "P3_clinically_tuned_RCR_BC"
    return model


def run_solver(solver: Path, model: dict, model_path: Path, result_path: Path) -> None:
    model_path.write_text(json.dumps(model, indent=2) + "\n")
    completed = subprocess.run([str(solver), str(model_path), str(result_path)], capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip() or "svzerodsolver failed")


def metrics_from_result(model: dict, result_path: Path) -> dict[str, float]:
    rows_by_vessel: dict[str, list[dict[str, str]]] = defaultdict(list)
    with result_path.open() as stream:
        for row in csv.DictReader(stream):
            rows_by_vessel[row["name"]].append(row)
    mapping = attachments(model)

    def values(boundary_name: str) -> tuple[list[float], list[float]]:
        vessel, side = mapping[boundary_name]
        rows = rows_by_vessel[vessel][-STEPS_PER_CYCLE:]
        flow_key = "flow_in" if side == "inlet" else "flow_out"
        pressure_key = "pressure_in" if side == "inlet" else "pressure_out"
        return [float(row[flow_key]) for row in rows], [float(row[pressure_key]) * MMHG_PER_DYN_CM2 for row in rows]

    lpa_flow = [sum(flow[index] for flow, _ in (values(f"RCR_{face}") for face in LPA_FACES)) for index in range(STEPS_PER_CYCLE)]
    rpa_flow = [sum(flow[index] for flow, _ in (values(f"RCR_{face}") for face in RPA_FACES)) for index in range(STEPS_PER_CYCLE)]
    _, fontan_pressure = values(f"FLOW_{FONTAN_FACE}")
    lpa, rpa = fmean(lpa_flow), fmean(rpa_flow)
    return {
        "mean_fontan_pressure_mmHg": fmean(fontan_pressure),
        "mean_lpa_flow_ml_s": lpa,
        "mean_rpa_flow_ml_s": rpa,
        "lpa_flow_fraction": lpa / (lpa + rpa),
        "rpa_flow_fraction": rpa / (lpa + rpa),
    }


def objective(metrics: dict[str, float], targets: dict[str, float]) -> float:
    pressure = metrics["mean_fontan_pressure_mmHg"] / targets["mean_fontan_pressure_mmHg"] - 1.0
    split = metrics["lpa_flow_fraction"] / targets["lpa_flow_fraction"] - 1.0
    return pressure * pressure + split * split


def nelder_mead(evaluate, initial: tuple[float, float], max_iterations: int = 10, tolerance: float = 1.0e-6) -> tuple[tuple[float, float], float]:
    """Small bounded, dependency-free Nelder-Mead implementation for two controls."""
    simplex = [clamp(initial), clamp((initial[0] + 0.30, initial[1])), clamp((initial[0], initial[1] + 0.30))]
    values = [evaluate(point) for point in simplex]
    for _ in range(max_iterations):
        order = sorted(range(3), key=lambda index: values[index])
        simplex, values = [simplex[index] for index in order], [values[index] for index in order]
        if max(values) - min(values) < tolerance:
            break
        best, second, worst = simplex
        centroid = ((best[0] + second[0]) / 2.0, (best[1] + second[1]) / 2.0)
        reflected = clamp((centroid[0] + (centroid[0] - worst[0]), centroid[1] + (centroid[1] - worst[1])))
        reflected_value = evaluate(reflected)
        if values[0] <= reflected_value < values[1]:
            simplex[2], values[2] = reflected, reflected_value
            continue
        if reflected_value < values[0]:
            expanded = clamp((centroid[0] + 2.0 * (reflected[0] - centroid[0]), centroid[1] + 2.0 * (reflected[1] - centroid[1])))
            expanded_value = evaluate(expanded)
            simplex[2], values[2] = (expanded, expanded_value) if expanded_value < reflected_value else (reflected, reflected_value)
            continue
        contracted = clamp((centroid[0] + 0.5 * (worst[0] - centroid[0]), centroid[1] + 0.5 * (worst[1] - centroid[1])))
        contracted_value = evaluate(contracted)
        if contracted_value < values[2]:
            simplex[2], values[2] = contracted, contracted_value
            continue
        simplex[1] = clamp(((best[0] + second[0]) / 2.0, (best[1] + second[1]) / 2.0))
        simplex[2] = clamp(((best[0] + worst[0]) / 2.0, (best[1] + worst[1]) / 2.0))
        values[1], values[2] = evaluate(simplex[1]), evaluate(simplex[2])
    best_index = min(range(3), key=lambda index: values[index])
    return simplex[best_index], values[best_index]


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    source = json.loads(SOURCE_MODEL.read_text())
    targets = json.loads(TARGET_FILE.read_text())
    solver = find_solver()
    history, cache = [], {}

    def evaluate(position: tuple[float, float]) -> float:
        position = clamp(position)
        key = tuple(round(value, 12) for value in position)
        if key in cache:
            return cache[key][1]
        model = apply_rcr_scales(source, position)
        with tempfile.TemporaryDirectory(dir=OUTPUT, prefix="candidate_") as directory:
            directory = Path(directory)
            result = directory / "result.csv"
            run_solver(solver, model, directory / "model.json", result)
            metrics = metrics_from_result(model, result)
        value = objective(metrics, targets)
        scales = {"common_rcr_resistance_scale": math.exp(position[0]), "lpa_rcr_resistance_scale": math.exp(position[0] + position[1]), "rpa_rcr_resistance_scale": math.exp(position[0] - position[1])}
        history.append({"evaluation": len(history) + 1, "objective": value, **scales, **metrics})
        cache[key] = (metrics, value)
        return value

    baseline_position = (0.0, 0.0)
    baseline_value = evaluate(baseline_position)
    best_position, best_value = nelder_mead(evaluate, baseline_position)
    best_model = apply_rcr_scales(source, best_position)
    tuned_json, tuned_csv = OUTPUT / "P3_0d_tuned_bc.json", OUTPUT / "P3_0d_tuned_bc.csv"
    run_solver(solver, best_model, tuned_json, tuned_csv)
    baseline_metrics = cache[tuple(round(value, 12) for value in baseline_position)][0]
    tuned_metrics = metrics_from_result(best_model, tuned_csv)
    write_csv(OUTPUT / "tuning_history.csv", history)
    errors = {
        "mean_fontan_pressure": abs(tuned_metrics["mean_fontan_pressure_mmHg"] / targets["mean_fontan_pressure_mmHg"] - 1.0),
        "mean_lpa_flow": abs(tuned_metrics["mean_lpa_flow_ml_s"] / (targets["lpa_flow_fraction"] * (tuned_metrics["mean_lpa_flow_ml_s"] + tuned_metrics["mean_rpa_flow_ml_s"])) - 1.0),
        "mean_rpa_flow": abs(tuned_metrics["mean_rpa_flow_ml_s"] / (targets["rpa_flow_fraction"] * (tuned_metrics["mean_lpa_flow_ml_s"] + tuned_metrics["mean_rpa_flow_ml_s"])) - 1.0),
    }
    summary = {
        "source_model": str(SOURCE_MODEL), "method": "bounded two-parameter Nelder-Mead RCR tuning", "targets": targets,
        "baseline": baseline_metrics, "tuned": tuned_metrics, "best_objective": best_value,
        "scales": {"common_rcr_resistance_scale": math.exp(best_position[0]), "lpa_rcr_resistance_scale": math.exp(best_position[0] + best_position[1]), "rpa_rcr_resistance_scale": math.exp(best_position[0] - best_position[1])},
        "relative_errors": errors, "number_of_solver_evaluations": len(history),
        "acceptance": {"fontan_pressure": errors["mean_fontan_pressure"] <= targets["mean_fontan_pressure_relative_tolerance"], "lpa_flow": errors["mean_lpa_flow"] <= targets["pulmonary_flow_relative_tolerance"], "rpa_flow": errors["mean_rpa_flow"] <= targets["pulmonary_flow_relative_tolerance"]},
    }
    (OUTPUT / "tuning_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = ["# P3 Clinical RCR BC Tuning", "", "| Metric | Clinical target | Calibrated-geometry baseline | Tuned 0D | Relative error | Status |", "| --- | ---: | ---: | ---: | ---: | --- |"]
    metrics = (("Mean Fontan pressure", "mean_fontan_pressure_mmHg", "mmHg", targets["mean_fontan_pressure_mmHg"], "fontan_pressure"), ("LPA flow", "mean_lpa_flow_ml_s", "mL/s", targets["lpa_flow_fraction"] * (tuned_metrics["mean_lpa_flow_ml_s"] + tuned_metrics["mean_rpa_flow_ml_s"]), "lpa_flow"), ("RPA flow", "mean_rpa_flow_ml_s", "mL/s", targets["rpa_flow_fraction"] * (tuned_metrics["mean_lpa_flow_ml_s"] + tuned_metrics["mean_rpa_flow_ml_s"]), "rpa_flow"))
    for label, key, unit, target, status_key in metrics:
        error = errors["mean_fontan_pressure"] if status_key == "fontan_pressure" else errors[f"mean_{status_key}"]
        lines.append(f"| {label} | {target:.4f} {unit} | {baseline_metrics[key]:.4f} {unit} | {tuned_metrics[key]:.4f} {unit} | {error * 100.0:.3f}% | {'PASS' if summary['acceptance'][status_key] else 'FAIL'} |")
    lines.extend(["", "## Tuned RCR resistance scales", "", f"- Common scale: `{summary['scales']['common_rcr_resistance_scale']:.6g}`", f"- LPA RCR scale: `{summary['scales']['lpa_rcr_resistance_scale']:.6g}`", f"- RPA RCR scale: `{summary['scales']['rpa_rcr_resistance_scale']:.6g}`", "- RCR capacitance, distal pressure, prescribed inflows, and all geometric vessel parameters are unchanged."])
    (OUTPUT / "tuning_summary.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
