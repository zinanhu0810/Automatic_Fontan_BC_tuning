#!/usr/bin/env python3
"""Calibrate P3 geometric 0D vessels to final-cycle unsteady 3D metrics."""

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
SOURCE_MODEL = PIPELINE / "02_generate_0d" / "rom" / "P3_0d.json"
THREE_D_RESULTS = PIPELINE / "01_initial_3d" / "unsteady" / "48-procs"
THREE_D_PRESSURE = THREE_D_RESULTS / "B_NS_Pressure_average.txt"
THREE_D_FLOW = THREE_D_RESULTS / "B_NS_Velocity_flux.txt"
OUTPUT = HERE / "unsteady"
MMHG_PER_DYN_CM2 = 1.0 / 1333.22
STEPS_PER_CYCLE = 1000
NUMBER_OF_CYCLES = 9
MATCH_RELATIVE_TOLERANCE = 0.001
LPA_FACES = ("cap_LPA_1", "cap_LPA_2", "cap_LPA_3", "cap_LPA_4", "cap_LPA_5", "cap_LPA_6", "cap_LPA_7", "cap_LPA_8", "cap_PA_main_2")
RPA_FACES = ("cap_PA_main", "cap_RPA_1", "cap_RPA_2", "cap_RPA_3", "cap_RPA_4", "cap_RPA_5", "cap_RPA_6", "cap_RPA_7", "cap_RPA_8")
FONTAN_FACE = "cap_Fontan"


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


def read_integrals(path: Path) -> list[dict[str, float]]:
    lines = path.read_text().splitlines()
    header = next(index for index, line in enumerate(lines) if line.startswith("step"))
    names = lines[header].split()
    return [dict(zip(names, map(float, line.split()))) for line in lines[header + 1 :] if line.split()]


def three_d_final_cycle() -> list[dict[str, float]]:
    pressure = {row["step"]: row for row in read_integrals(THREE_D_PRESSURE)}
    flow = {row["step"]: row for row in read_integrals(THREE_D_FLOW)}
    steps = sorted(set(pressure) & set(flow))[-STEPS_PER_CYCLE:]
    start_time = pressure[steps[0]]["time"]
    return [{
        "step": step,
        "time_in_window_s": pressure[step]["time"] - start_time,
        "fontan_pressure_mmHg": pressure[step][FONTAN_FACE] * MMHG_PER_DYN_CM2,
        "lpa_flow_ml_s": sum(flow[step][face] for face in LPA_FACES),
        "rpa_flow_ml_s": sum(flow[step][face] for face in RPA_FACES),
    } for step in steps]


def target_metrics() -> dict[str, float]:
    final = three_d_final_cycle()
    lpa, rpa = fmean(row["lpa_flow_ml_s"] for row in final), fmean(row["rpa_flow_ml_s"] for row in final)
    pressure = [row["fontan_pressure_mmHg"] for row in final]
    return {
        "mean_fontan_pressure_mmHg": fmean(pressure),
        "fontan_pressure_amplitude_mmHg": max(pressure) - min(pressure),
        "mean_lpa_flow_ml_s": lpa,
        "mean_rpa_flow_ml_s": rpa,
        "lpa_split_fraction": lpa / (lpa + rpa),
        "rpa_split_fraction": rpa / (lpa + rpa),
    }


def prepare_model(source: Path, scales: tuple[float, float, float]) -> dict:
    """Scale only vessel R/L; preserve every FLOW/RCR boundary unchanged."""
    global_r, lpa_to_rpa_r, vessel_l = scales
    original = json.loads(source.read_text())
    model = deepcopy(original)
    initial_boundaries = original["boundary_conditions"]
    model_attachments = attachments(model)
    lpa_vessels = {model_attachments[f"RCR_{face}"][0] for face in LPA_FACES}
    rpa_vessels = {model_attachments[f"RCR_{face}"][0] for face in RPA_FACES}

    for vessel in model["vessels"]:
        values = vessel["zero_d_element_values"]
        resistance_scale = global_r
        if vessel["vessel_name"] in lpa_vessels:
            resistance_scale *= lpa_to_rpa_r
        elif vessel["vessel_name"] in rpa_vessels:
            resistance_scale /= lpa_to_rpa_r
        values["R_poiseuille"] *= resistance_scale
        values["stenosis_coefficient"] *= resistance_scale
        values["L"] *= vessel_l
        values["C"] = 0.0

    parameters = model["simulation_parameters"]
    parameters["model_name"] = "P3_3D_calibrated_rigid_geometry"
    parameters["number_of_cardiac_cycles"] = 12
    parameters["number_of_time_pts_per_cardiac_cycle"] = STEPS_PER_CYCLE
    if model["boundary_conditions"] != initial_boundaries:
        raise RuntimeError("Boundary conditions changed during vessel-only calibration")
    return model


def run_solver(solver: Path, model: dict, model_path: Path, result_path: Path) -> None:
    model_path.write_text(json.dumps(model, indent=2) + "\n")
    completed = subprocess.run([str(solver), str(model_path), str(result_path)], capture_output=True, text=True, timeout=20)
    if completed.returncode:
        raise RuntimeError(f"svzerodsolver failed: {completed.stderr.strip() or completed.stdout.strip()}")


def result_metrics(model: dict, result_path: Path) -> dict[str, float]:
    rows_by_vessel: dict[str, list[dict[str, str]]] = defaultdict(list)
    with result_path.open() as stream:
        for row in csv.DictReader(stream):
            rows_by_vessel[row["name"]].append(row)
    mapping = attachments(model)

    def boundary_values(boundary_name: str) -> tuple[list[float], list[float]]:
        vessel, side = mapping[boundary_name]
        rows = rows_by_vessel[vessel][-STEPS_PER_CYCLE:]
        flow_key = "flow_in" if side == "inlet" else "flow_out"
        pressure_key = "pressure_in" if side == "inlet" else "pressure_out"
        return [float(row[flow_key]) for row in rows], [float(row[pressure_key]) * MMHG_PER_DYN_CM2 for row in rows]

    lpa_flow = [sum(values[index] for values, _ in (boundary_values(f"RCR_{face}") for face in LPA_FACES)) for index in range(STEPS_PER_CYCLE)]
    rpa_flow = [sum(values[index] for values, _ in (boundary_values(f"RCR_{face}") for face in RPA_FACES)) for index in range(STEPS_PER_CYCLE)]
    _, fontan_pressure = boundary_values(f"FLOW_{FONTAN_FACE}")
    lpa, rpa = fmean(lpa_flow), fmean(rpa_flow)
    return {
        "mean_fontan_pressure_mmHg": fmean(fontan_pressure),
        "fontan_pressure_amplitude_mmHg": max(fontan_pressure) - min(fontan_pressure),
        "min_fontan_pressure_mmHg": min(fontan_pressure),
        "max_fontan_pressure_mmHg": max(fontan_pressure),
        "mean_lpa_flow_ml_s": lpa,
        "mean_rpa_flow_ml_s": rpa,
        "lpa_split_fraction": lpa / (lpa + rpa),
        "rpa_split_fraction": rpa / (lpa + rpa),
    }


def residuals(metrics: dict[str, float], target: dict[str, float]) -> list[float]:
    return [
        metrics["mean_fontan_pressure_mmHg"] / target["mean_fontan_pressure_mmHg"] - 1.0,
        metrics["lpa_split_fraction"] / target["lpa_split_fraction"] - 1.0,
        metrics["fontan_pressure_amplitude_mmHg"] / target["fontan_pressure_amplitude_mmHg"] - 1.0,
    ]


def objective(metrics: dict[str, float], target: dict[str, float]) -> float:
    return sum(value * value for value in residuals(metrics, target))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_markdown_summary(path: Path, target: dict[str, float], baseline: dict[str, float], calibrated: dict[str, float], scales: tuple[float, float, float]) -> None:
    rows = (
        ("Mean Fontan pressure", target["mean_fontan_pressure_mmHg"], baseline["mean_fontan_pressure_mmHg"], calibrated["mean_fontan_pressure_mmHg"], "mmHg"),
        ("LPA flow split", target["lpa_split_fraction"] * 100.0, baseline["lpa_split_fraction"] * 100.0, calibrated["lpa_split_fraction"] * 100.0, "%"),
        ("Fontan pressure amplitude", target["fontan_pressure_amplitude_mmHg"], baseline["fontan_pressure_amplitude_mmHg"], calibrated["fontan_pressure_amplitude_mmHg"], "mmHg"),
    )
    lines = [
        "# P3 Unsteady 0D Calibration",
        "",
        "The 0D geometric vessels were fit to the final 3D cycle. All prescribed FLOW and RCR boundary conditions are unchanged; vessel capacitances are zero for the rigid-wall comparison.",
        "",
        "| Metric | 3D target | Rigid 0D baseline | Calibrated 0D | Relative error | Status |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for label, target_value, baseline_value, calibrated_value, unit in rows:
        error = abs(calibrated_value / target_value - 1.0) * 100.0
        status = "PASS" if error <= MATCH_RELATIVE_TOLERANCE * 100.0 else "FAIL"
        lines.append(f"| {label} | {target_value:.4f} {unit} | {baseline_value:.4f} {unit} | {calibrated_value:.4f} {unit} | {error:.3f}% | {status} |")
    lines.extend([
        "",
        "## Vessel-only scales",
        "",
        f"- Global vessel resistance: `{scales[0]:.6g}`",
        f"- LPA/RPA terminal-vessel resistance balance: `{scales[1]:.6g}`",
        f"- Global vessel inductance: `{scales[2]:.6g}`",
        "- Geometric vessel capacitance: `0.0`",
        "- Acceptance criterion: relative error `<= 0.1%` for each fitted metric.",
    ])
    path.write_text("\n".join(lines) + "\n")


def write_fontan_waveform_svg(path: Path, model: dict, result_path: Path) -> None:
    """Write a dependency-free overlay of the final 3D and 0D Fontan cycles."""
    three_d_rows = three_d_final_cycle()
    three_d_time = [float(row["time_in_window_s"]) for row in three_d_rows]
    three_d_pressure = [float(row["fontan_pressure_mmHg"]) for row in three_d_rows]

    rows_by_vessel: dict[str, list[dict[str, str]]] = defaultdict(list)
    with result_path.open() as stream:
        for row in csv.DictReader(stream):
            rows_by_vessel[row["name"]].append(row)
    vessel, side = attachments(model)[f"FLOW_{FONTAN_FACE}"]
    rows = rows_by_vessel[vessel][-STEPS_PER_CYCLE:]
    pressure_key = "pressure_in" if side == "inlet" else "pressure_out"
    zero_d_time = [float(row["time"]) for row in rows]
    zero_d_time = [value - zero_d_time[0] for value in zero_d_time]
    zero_d_pressure = [float(row[pressure_key]) * MMHG_PER_DYN_CM2 for row in rows]

    width, height = 980, 500
    left, right, top, bottom = 100, 40, 35, 85
    plot_width, plot_height = width - left - right, height - top - bottom
    y_min = min(three_d_pressure + zero_d_pressure)
    y_max = max(three_d_pressure + zero_d_pressure)
    y_pad = max(0.15, 0.08 * (y_max - y_min))
    y_min, y_max = y_min - y_pad, y_max + y_pad
    duration = max(three_d_time[-1], zero_d_time[-1])

    def points(time: list[float], pressure: list[float]) -> str:
        return " ".join(
            f"{left + plot_width * value / duration:.2f},{top + plot_height * (1.0 - (press - y_min) / (y_max - y_min)):.2f}"
            for value, press in zip(time, pressure)
        )

    ticks = 5
    grid = []
    labels = []
    for index in range(ticks):
        value = y_min + index * (y_max - y_min) / (ticks - 1)
        y = top + plot_height * (1.0 - (value - y_min) / (y_max - y_min))
        grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}" class="grid"/>')
        labels.append(f'<text x="{left - 12}" y="{y + 5:.2f}" class="tick" text-anchor="end">{value:.2f}</text>')
    x_ticks = []
    for index in range(5):
        value = index * duration / 4.0
        x = left + plot_width * index / 4.0
        x_ticks.append(f'<line x1="{x:.2f}" y1="{height - bottom}" x2="{x:.2f}" y2="{height - bottom + 6}" class="axis"/><text x="{x:.2f}" y="{height - bottom + 28}" class="tick" text-anchor="middle">{value:.3f}</text>')

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>.title{{font:600 20px Arial,sans-serif;fill:#111827}}.label{{font:14px Arial,sans-serif;fill:#374151}}.tick{{font:12px Arial,sans-serif;fill:#4b5563}}.axis{{stroke:#374151;stroke-width:1.4}}.grid{{stroke:#d1d5db;stroke-width:1}}.legend{{font:14px Arial,sans-serif;fill:#111827}}</style>
<rect width="100%" height="100%" fill="white"/>
<text x="{left}" y="25" class="title">P3 final-cycle Fontan pressure: 3D CFD vs calibrated 0D</text>
{''.join(grid)}
<line x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}" class="axis"/><line x1="{left}" y1="{height - bottom}" x2="{width - right}" y2="{height - bottom}" class="axis"/>
{''.join(labels)}{''.join(x_ticks)}
<polyline fill="none" stroke="#2563eb" stroke-width="2.2" points="{points(three_d_time, three_d_pressure)}"/>
<polyline fill="none" stroke="#dc2626" stroke-width="2.2" points="{points(zero_d_time, zero_d_pressure)}"/>
<line x1="{width - 310}" y1="58" x2="{width - 280}" y2="58" stroke="#2563eb" stroke-width="2.2"/><text x="{width - 270}" y="63" class="legend">3D CFD</text>
<line x1="{width - 170}" y1="58" x2="{width - 140}" y2="58" stroke="#dc2626" stroke-width="2.2"/><text x="{width - 130}" y="63" class="legend">Calibrated 0D</text>
<text x="{left - 62}" y="{top + plot_height / 2}" class="label" text-anchor="middle" transform="rotate(-90 {left - 62} {top + plot_height / 2})">Fontan pressure [mmHg]</text>
<text x="{left + plot_width / 2}" y="{height - 22}" class="label" text-anchor="middle">Time within cardiac cycle [s]</text>
</svg>'''
    path.write_text(svg)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    solver, target = find_solver(), target_metrics()
    history: list[dict] = []
    cache: dict[tuple[float, float, float], tuple[dict[str, float], float]] = {}

    def evaluate(scales: tuple[float, float, float]) -> tuple[dict[str, float], float]:
        key = tuple(round(value, 10) for value in scales)
        if key in cache:
            return cache[key]
        model = prepare_model(SOURCE_MODEL, scales)
        with tempfile.TemporaryDirectory(dir=OUTPUT, prefix="candidate_") as temporary:
            result = Path(temporary) / "result.csv"
            run_solver(solver, model, Path(temporary) / "model.json", result)
            metrics = result_metrics(model, result)
        cost = objective(metrics, target)
        history.append({
            "evaluation": len(history) + 1,
            "global_vessel_resistance_scale": scales[0],
            "lpa_to_rpa_terminal_resistance_scale": scales[1],
            "global_vessel_inductance_scale": scales[2],
            "objective": cost,
            **metrics,
        })
        cache[key] = (metrics, cost)
        return cache[key]

    # P3's geometry has a higher vessel-resistance response than P1. Start in
    # the verified local region; the search remains multiplicative and keeps
    # all R/L values positive.
    current = (16.37, 0.533, 2.879)
    _, current_cost = evaluate(current)
    # These ranges cover the local P3 pressure/split response while avoiding
    # nonphysical trials that can make the 0D solve unnecessarily stiff.
    bounds = ((0.3, 30.0), (0.5, 2.0), (0.2, 5.0))
    steps = [math.log(1.1), math.log(1.05), math.log(1.1)]
    for _ in range(3):
        improved = False
        for index in range(3):
            candidates = [current]
            for direction in (-1.0, 1.0):
                candidate = list(current)
                candidate[index] = min(bounds[index][1], max(bounds[index][0], current[index] * math.exp(direction * steps[index])))
                candidates.append(tuple(candidate))
            best = min(candidates, key=lambda candidate: evaluate(candidate)[1])
            _, best_cost = evaluate(best)
            if best_cost + 1.0e-12 < current_cost:
                current, current_cost, improved = best, best_cost, True
        steps = [value * 0.55 for value in steps]
        if not improved:
            break

    baseline_model, calibrated_model = prepare_model(SOURCE_MODEL, (1.0, 1.0, 1.0)), prepare_model(SOURCE_MODEL, current)
    baseline_json, baseline_csv = OUTPUT / "P3_0d_rigid_baseline.json", OUTPUT / "P3_0d_rigid_baseline.csv"
    calibrated_json, calibrated_csv = OUTPUT / "P3_0d_calibrated.json", OUTPUT / "P3_0d_calibrated.csv"
    run_solver(solver, baseline_model, baseline_json, baseline_csv)
    run_solver(solver, calibrated_model, calibrated_json, calibrated_csv)
    baseline, calibrated = result_metrics(baseline_model, baseline_csv), result_metrics(calibrated_model, calibrated_csv)
    rows = []
    for label, key, unit, scale in (
        ("Mean Fontan pressure", "mean_fontan_pressure_mmHg", "mmHg", 1.0),
        ("LPA flow split", "lpa_split_fraction", "%", 100.0),
        ("Fontan pressure amplitude", "fontan_pressure_amplitude_mmHg", "mmHg", 1.0),
    ):
        rows.append({"metric": label, "unit": unit, "3D_target": target[key] * scale, "0D_baseline": baseline[key] * scale, "0D_calibrated": calibrated[key] * scale, "calibrated_relative_error": abs(calibrated[key] / target[key] - 1.0)})
    write_csv(OUTPUT / "calibration_metrics.csv", rows)
    write_csv(OUTPUT / "calibration_history.csv", history)
    write_markdown_summary(OUTPUT / "calibration_summary.md", target, baseline, calibrated, current)
    write_fontan_waveform_svg(OUTPUT / "final_cycle_fontan_pressure_3d_vs_0d.svg", calibrated_model, calibrated_csv)
    (OUTPUT / "calibration_targets.json").write_text(json.dumps(target, indent=2) + "\n")
    summary = {
        "target_source": str(THREE_D_RESULTS), "calibration_method": "vessel-only coordinate search against final 3D cycle",
        "parameters": {"global_vessel_resistance_scale": current[0], "lpa_to_rpa_terminal_resistance_scale": current[1], "global_vessel_inductance_scale": current[2], "geometric_vessel_capacitance": 0.0, "boundary_conditions_preserved": True},
        "baseline": baseline, "calibrated": calibrated, "objective": current_cost, "evaluations": len(history),
    }
    (OUTPUT / "calibration_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"target": target, "calibrated": calibrated, "parameters": summary["parameters"]}, indent=2))


if __name__ == "__main__":
    main()
