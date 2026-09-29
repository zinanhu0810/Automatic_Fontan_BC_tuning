#!/usr/bin/env python3
"""Calibrate a rigid geometric 0D model to the final 3D CFD cycle.

Fits three vessel-only scale factors while preserving every solver.xml boundary
condition: global resistance for mean Fontan pressure, LPA/RPA resistance balance
for pulmonary split, and global inductance for Fontan-pressure amplitude.
"""

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
TOLERANCE = 0.001


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--model", type=Path, default=None, help="Default: CASE/0d/CASE_0d.json")
    parser.add_argument("--results", type=Path, default=Path("."), help="Directory with B_NS_* files; default: case directory")
    parser.add_argument("--output-directory", default="calibrated_0d")
    parser.add_argument("--svzerodsolver", type=Path, default=None)
    parser.add_argument("--max-evaluations", type=int, default=100)
    return parser.parse_args()


def cfd_rows(path: Path) -> list[dict[str, float]]:
    lines = path.read_text().splitlines()
    start = next(index for index, line in enumerate(lines) if line.startswith("step"))
    header = lines[start].split()
    start += 1
    return [dict(zip(header, map(float, line.split()))) for line in lines[start:] if line.split()]


def final_cfd_cycle(results: Path, config: dict) -> list[dict[str, float]]:
    pressure = {row["step"]: row for row in cfd_rows(results / "B_NS_Pressure_average.txt")}
    flow = {row["step"]: row for row in cfd_rows(results / "B_NS_Velocity_flux.txt")}
    steps = sorted(set(pressure) & set(flow))[-int(config["simulation"]["steps_per_cycle"]):]
    if len(steps) != int(config["simulation"]["steps_per_cycle"]): raise ValueError("No complete final CFD cycle is available.")
    fontan = config["metrics"]["fontan_pressure_face"]
    lpa, rpa = config["outlets"]["lpa"], config["outlets"]["rpa"]
    t0 = pressure[steps[0]]["time"]
    return [{"time": pressure[step]["time"] - t0, "pressure": pressure[step][fontan] * MMHG_PER_DYN_CM2, "lpa": sum(flow[step][name] for name in lpa), "rpa": sum(flow[step][name] for name in rpa)} for step in steps]


def cfd_metrics(rows: list[dict[str, float]]) -> dict[str, float]:
    lpa, rpa = fmean(row["lpa"] for row in rows), fmean(row["rpa"] for row in rows)
    values = [row["pressure"] for row in rows]
    return {"mean_pressure": fmean(values), "pressure_amplitude": max(values) - min(values), "lpa_fraction": lpa / (lpa + rpa)}


def solver_path(explicit: Path | None) -> Path:
    choices = [explicit, Path(shutil.which("svzerodsolver")) if shutil.which("svzerodsolver") else None, Path.home() / "anaconda3/bin/svzerodsolver"]
    solver = next((item.resolve() for item in choices if item and item.is_file()), None)
    if solver is None: raise FileNotFoundError("svzerodsolver was not found; pass --svzerodsolver PATH.")
    return solver


def attachments(model: dict) -> dict[str, tuple[str, str]]:
    return {name: (vessel["vessel_name"], side) for vessel in model["vessels"] for side, name in vessel.get("boundary_conditions", {}).items()}


def terminal_boundary_name(mapping: dict[str, tuple[str, str]], cap: str) -> str:
    for prefix in ("RCR_", "RESISTANCE_"):
        name = prefix + cap
        if name in mapping:
            return name
    raise KeyError(f"No RCR or RESISTANCE boundary is attached to outlet '{cap}'.")


def scaled_model(original: dict, lpa_caps: list[str], rpa_caps: list[str], scales: tuple[float, float, float]) -> dict:
    resistance, balance, inductance = scales; model = deepcopy(original); original_bc = original["boundary_conditions"]
    attached = attachments(model)
    lpa_vessels = {attached[terminal_boundary_name(attached, cap)][0] for cap in lpa_caps}
    rpa_vessels = {attached[terminal_boundary_name(attached, cap)][0] for cap in rpa_caps}
    for vessel in model["vessels"]:
        values = vessel["zero_d_element_values"]; factor = resistance
        if vessel["vessel_name"] in lpa_vessels: factor *= balance
        if vessel["vessel_name"] in rpa_vessels: factor /= balance
        values["R_poiseuille"] *= factor; values["stenosis_coefficient"] *= factor
        values["L"] *= inductance; values["C"] = 0.0
    model["simulation_parameters"]["number_of_cardiac_cycles"] = max(12, int(model["simulation_parameters"].get("number_of_cardiac_cycles", 8)))
    if model["boundary_conditions"] != original_bc: raise RuntimeError("Boundary conditions changed during vessel-only calibration.")
    return model


def run(solver: Path, model: dict, folder: Path) -> Path:
    json_path, csv_path = folder / "model.json", folder / "result.csv"
    json_path.write_text(json.dumps(model, indent=2) + "\n")
    result = subprocess.run([str(solver), str(json_path), str(csv_path)], capture_output=True, text=True)
    if result.returncode: raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    return csv_path


def zero_d_cycle(model: dict, result: Path, config: dict) -> list[dict[str, float]]:
    by_name: dict[str, list[dict[str, str]]] = defaultdict(list)
    with result.open() as stream:
        for row in csv.DictReader(stream): by_name[row["name"]].append(row)
    count = int(config["simulation"]["steps_per_cycle"]); mapping = attachments(model)
    def bc(name: str) -> tuple[list[float], list[float], list[float]]:
        vessel, side = mapping[terminal_boundary_name(mapping, name)]; rows = by_name[vessel][-count:]
        flow = "flow_in" if side == "inlet" else "flow_out"; pressure = "pressure_in" if side == "inlet" else "pressure_out"
        return [float(row["time"]) for row in rows], [float(row[flow]) for row in rows], [float(row[pressure]) * MMHG_PER_DYN_CM2 for row in rows]
    root = config["metrics"]["fontan_pressure_face"]; vessel, side = mapping[f"FLOW_{root}"]
    rows = by_name[vessel][-count:]; pressure = "pressure_in" if side == "inlet" else "pressure_out"
    lpa = [bc(cap)[1] for cap in config["outlets"]["lpa"]]; rpa = [bc(cap)[1] for cap in config["outlets"]["rpa"]]
    time = [float(row["time"]) for row in rows]; return [{"time": value - time[0], "pressure": float(row[pressure]) * MMHG_PER_DYN_CM2, "lpa": sum(values[index] for values in lpa), "rpa": sum(values[index] for values in rpa)} for index, (value, row) in enumerate(zip(time, rows))]


def objective(metrics: dict[str, float], target: dict[str, float]) -> float:
    errors = [abs(metrics[key] / target[key] - 1.0) for key in target]
    # Match the acceptance rule: no metric may hide behind excellent agreement
    # in another metric. The small secondary term breaks ties smoothly.
    return max(errors) + 0.01 * sum(error**2 for error in errors)


def plot(path: Path, cfd: list[dict[str, float]], zero_d: list[dict[str, float]]) -> None:
    width, height, left, right, top, bottom = 900, 480, 90, 35, 55, 75
    duration = max(cfd[-1]["time"], zero_d[-1]["time"]); values = [item["pressure"] for item in cfd + zero_d]
    low, high = min(values), max(values); pad = max(0.1, 0.08 * (high - low)); low -= pad; high += pad
    x = lambda value: left + (width - left - right) * value / duration
    y = lambda value: top + (height - top - bottom) * (high - value) / (high - low)
    points = lambda rows: " ".join(f"{x(row['time']):.2f},{y(row['pressure']):.2f}" for row in rows)
    grid = "".join(f'<line x1="{left}" y1="{y(low + i * (high-low) / 4):.2f}" x2="{width-right}" y2="{y(low + i * (high-low) / 4):.2f}" class="grid"/>' for i in range(5))
    path.write_text(f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"><style>.grid{{stroke:#d1d5db}}.axis{{stroke:#374151}}text{{font-family:Arial;fill:#111827}}.tick{{font-size:12px}}</style><rect width="100%" height="100%" fill="white"/><text x="{left}" y="28" font-size="20" font-weight="600">Final-cycle Fontan pressure: 3D CFD vs calibrated 0D</text>{grid}<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" class="axis"/><line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" class="axis"/><polyline fill="none" stroke="#2563eb" stroke-width="2" points="{points(cfd)}"/><polyline fill="none" stroke="#dc2626" stroke-width="2" points="{points(zero_d)}"/><line x1="{width-270}" y1="52" x2="{width-240}" y2="52" stroke="#2563eb" stroke-width="2"/><text x="{width-232}" y="57" class="tick">3D CFD</text><line x1="{width-150}" y1="52" x2="{width-120}" y2="52" stroke="#dc2626" stroke-width="2"/><text x="{width-112}" y="57" class="tick">Calibrated 0D</text><text x="{left+(width-left-right)/2}" y="{height-25}" class="tick" text-anchor="middle">Time within cardiac cycle [s]</text><text x="20" y="{top+(height-top-bottom)/2}" class="tick" text-anchor="middle" transform="rotate(-90 20 {top+(height-top-bottom)/2})">Fontan pressure [mmHg]</text></svg>''')


def main() -> None:
    args = arguments(); case = args.case_directory.resolve(); config = json.loads((case / "step0_config.json").read_text())
    results = args.results if args.results.is_absolute() else case / args.results; output = Path(args.output_directory); output = output if output.is_absolute() else case / output; output.mkdir(parents=True, exist_ok=True)
    model_path = case / "0d" / f"{config['case_name']}_0d.json" if args.model is None else (args.model if args.model.is_absolute() else case / args.model)
    original = json.loads(model_path.read_text()); solver = solver_path(args.svzerodsolver)
    cfd = final_cfd_cycle(results, config); target = cfd_metrics(cfd)
    # A steady, constant-inflow CFD run has no physiological pressure pulse.
    # Its small residual numerical ripple cannot identify vessel inductance.
    is_steady = config.get("outlet_boundary_type") == "resistance" and int(config["simulation"].get("number_of_cycles", 0)) == 1
    if is_steady:
        target.pop("pressure_amplitude")
    cache: dict[tuple[float, float, float], tuple[dict[str, float], list[dict[str, float]], float]] = {}; evaluations = 0
    if args.max_evaluations < 4:
        raise ValueError("Nelder-Mead calibration requires at least four evaluations.")

    class EvaluationLimit(RuntimeError):
        pass

    def evaluate(scales: tuple[float, float, float]):
        nonlocal evaluations
        key = tuple(round(value, 9) for value in scales)
        if key not in cache:
            if evaluations >= args.max_evaluations:
                raise EvaluationLimit
            candidate = scaled_model(original, config["outlets"]["lpa"], config["outlets"]["rpa"], scales)
            with tempfile.TemporaryDirectory(prefix="calibrate_0d_") as folder:
                waveform = zero_d_cycle(candidate, run(solver, candidate, Path(folder)), config)
            metrics = cfd_metrics(waveform); cache[key] = (metrics, waveform, objective(metrics, target)); evaluations += 1
        return cache[key]
    bounds = ((0.05, 50.0), (0.1, 10.0), (0.05, 20.0))
    log_bounds = tuple((math.log(low), math.log(high)) for low, high in bounds)
    def scales_from_log(values: tuple[float, float, float]) -> tuple[float, float, float]:
        return tuple(math.exp(min(high, max(low, value))) for value, (low, high) in zip(values, log_bounds))
    def evaluation(values: tuple[float, float, float]):
        return evaluate(scales_from_log(values))
    def meets_tolerance(metrics: dict[str, float]) -> bool:
        return all(abs(metrics[key] / target[key] - 1.0) <= TOLERANCE for key in target)

    # Nelder-Mead can make coupled R/balance/L moves, unlike the former
    # coordinate search that could only perturb one calibration factor at once.
    origin = (0.0, 0.0, 0.0)
    initial_step = math.log(2.0)
    simplex = [(origin, evaluation(origin))]
    for index in range(3):
        point = list(origin); point[index] = initial_step
        simplex.append((tuple(point), evaluation(tuple(point))))
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5
    while evaluations < args.max_evaluations:
        simplex.sort(key=lambda item: item[1][2])
        if meets_tolerance(simplex[0][1][0]):
            break
        best, second_worst, worst = simplex[0], simplex[-2], simplex[-1]
        centroid = tuple(sum(point[index] for point, _ in simplex[:-1]) / 3.0 for index in range(3))
        try:
            reflected = tuple(centroid[index] + alpha * (centroid[index] - worst[0][index]) for index in range(3))
            reflected_value = evaluation(reflected)
            if reflected_value[2] < best[1][2]:
                expanded = tuple(centroid[index] + gamma * (reflected[index] - centroid[index]) for index in range(3))
                expanded_value = evaluation(expanded)
                simplex[-1] = (expanded, expanded_value) if expanded_value[2] < reflected_value[2] else (reflected, reflected_value)
            elif reflected_value[2] < second_worst[1][2]:
                simplex[-1] = (reflected, reflected_value)
            else:
                contracted = tuple(centroid[index] + rho * ((worst[0] if reflected_value[2] >= worst[1][2] else reflected)[index] - centroid[index]) for index in range(3))
                contracted_value = evaluation(contracted)
                if contracted_value[2] < min(worst[1][2], reflected_value[2]):
                    simplex[-1] = (contracted, contracted_value)
                else:
                    for index in range(1, len(simplex)):
                        point = tuple(best[0][axis] + sigma * (simplex[index][0][axis] - best[0][axis]) for axis in range(3))
                        simplex[index] = (point, evaluation(point))
        except EvaluationLimit:
            break
    simplex.sort(key=lambda item: item[1][2])
    current = scales_from_log(simplex[0][0])
    calibrated = scaled_model(original, config["outlets"]["lpa"], config["outlets"]["rpa"], current); metrics, waveform, cost = evaluate(current)
    output_model = output / f"{config['case_name']}_0d_calibrated.json"; output_model.write_text(json.dumps(calibrated, indent=2) + "\n")
    labels = [("Mean Fontan pressure", "mean_pressure", "mmHg", 1.0), ("LPA flow split", "lpa_fraction", "%", 100.0)]
    if "pressure_amplitude" in target:
        labels.append(("Fontan pressure amplitude", "pressure_amplitude", "mmHg", 1.0))
    lines = ["# 0D Calibration to Final 3D Cycle", "", "| Metric | 3D target | Calibrated 0D | Relative error | Status |", "| --- | ---: | ---: | ---: | --- |"]
    for label, key, unit, scale in labels:
        error = abs(metrics[key] / target[key] - 1.0); lines.append(f"| {label} | {target[key]*scale:.4f} {unit} | {metrics[key]*scale:.4f} {unit} | {error*100:.3f}% | {'PASS' if error <= TOLERANCE else 'FAIL'} |")
    lines.extend(["", "## Vessel Scales", "", f"- Global resistance: `{current[0]:.8g}`", f"- LPA/RPA resistance balance: `{current[1]:.8g}`", f"- Global inductance: `{current[2]:.8g}`", "- Geometric vessel capacitance: `0.0`"])
    if is_steady:
        lines.append("- Fontan pressure amplitude is omitted for this steady constant-inflow case; inductance is not identifiable.")
    lines.extend(["- Optimizer: `bounded Nelder-Mead (log-scaled parameters)`", f"- Evaluations: `{evaluations}`"])
    (output / "calibration_summary.md").write_text("\n".join(lines) + "\n")
    plot(output / "fontan_pressure_3d_vs_0d.svg", cfd, waveform)
    print("\n".join(lines))


if __name__ == "__main__": main()
