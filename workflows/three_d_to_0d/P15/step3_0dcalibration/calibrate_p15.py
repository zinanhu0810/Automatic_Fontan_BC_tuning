#!/usr/bin/env python3
"""Calibrate the P15 steady 0D vessel resistances to the steady 3D result."""

from __future__ import annotations

import argparse
import csv
import copy
import json
import math
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import fmean

import numpy as np
import vtk
from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
P15 = HERE.parent
DEFAULT_MODEL = P15 / "steady_0d" / "P15_0d_resistance_constant_flow.json"
DEFAULT_BASELINE_RESULT = P15 / "steady_0d" / "P15_0d_resistance_constant_flow.csv"
DEFAULT_CENTERLINE = P15 / "rom" / "centerlines.vtp"
DEFAULT_3D_RESULT = P15 / "steady" / "48-procs" / "result_4000.vtu"
DEFAULT_3D_INTEGRALS = P15 / "steady" / "48-procs"
MMHG_PER_DYN_CM2 = 1.0 / 1333.22

LEFT_PULMONARY_FACES = {
    "cap_l_pa_1_1_x",
    "cap_l_pa_1_x",
    "cap_l_pa_2_x",
    "cap_l_pa_3_1_x",
    "cap_l_pa_3_x",
    "cap_l_pa_4_x",
    "cap_l_pa_5_x",
    "cap_l_pa_7_1_x",
    "cap_l_pa_7_x",
    "cap_pa_x_2",
}
RIGHT_PULMONARY_FACES = {
    "cap_r_pa_1_2_x",
    "cap_r_pa_1_3_x",
    "cap_r_pa_1_x",
    "cap_r_pa_2_1_x",
    "cap_r_pa_2_x",
    "cap_r_pa_3_x",
    "cap_r_pa_4_x",
    "cap_r_pa_5_x",
    "cap_pa_x",
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def read_boundary_integral(path: Path) -> list[dict[str, float]]:
    lines = path.read_text().splitlines()
    header_index = next(
        index for index, line in enumerate(lines) if line.strip().startswith("step")
    )
    names = lines[header_index].split()
    rows = []
    for line in lines[header_index + 1 :]:
        fields = line.split()
        if not fields:
            continue
        if len(fields) != len(names):
            raise ValueError(f"Malformed boundary-integral row in {path}: {line}")
        rows.append({name: float(value) for name, value in zip(names, fields)})
    if not rows:
        raise ValueError(f"No boundary-integral data found in {path}")
    return rows


def final_window_boundary_means(directory: Path, number: int) -> tuple[dict, dict]:
    pressure_rows = read_boundary_integral(directory / "B_NS_Pressure_average.txt")
    flow_rows = read_boundary_integral(directory / "B_NS_Velocity_flux.txt")
    if len(pressure_rows) != len(flow_rows):
        raise ValueError("3D pressure and flow boundary histories have different lengths")
    number = min(number, len(flow_rows))
    pressure_rows = pressure_rows[-number:]
    flow_rows = flow_rows[-number:]
    for pressure, flow in zip(pressure_rows, flow_rows):
        if pressure["step"] != flow["step"]:
            raise ValueError("3D pressure and flow samples are not aligned")
    names = [name for name in flow_rows[0] if name not in {"step", "time", "walls_combined"}]
    pressures = {name: fmean(row[name] for row in pressure_rows) for name in names}
    flows = {name: fmean(row[name] for row in flow_rows) for name in names}
    return pressures, flows


def boundary_face_name(boundary_name: str) -> str:
    for prefix in ("FLOW_", "RESISTANCE_", "RCR_"):
        if boundary_name.startswith(prefix):
            return boundary_name[len(prefix) :]
    raise ValueError(f"Cannot infer 3D face from boundary name {boundary_name}")


def boundary_attachments(model: dict) -> dict[str, tuple[int, str]]:
    attachments = {}
    for vessel in model["vessels"]:
        vessel_id = int(vessel["vessel_id"])
        for side, boundary_name in vessel.get("boundary_conditions", {}).items():
            if boundary_name in attachments:
                raise ValueError(f"Boundary {boundary_name} is attached more than once")
            attachments[boundary_name] = (vessel_id, side)
    return attachments


def solve_steady_vessel_flows(model: dict, face_flows: dict[str, float]) -> dict[int, float]:
    """Recover every directed vessel flow from terminal 3D fluxes and continuity."""
    vessels = model["vessels"]
    vessel_ids = [int(vessel["vessel_id"]) for vessel in vessels]
    columns = {vessel_id: index for index, vessel_id in enumerate(vessel_ids)}
    rows = []
    rhs = []

    for boundary_name, (vessel_id, side) in boundary_attachments(model).items():
        face = boundary_face_name(boundary_name)
        row = np.zeros(len(vessel_ids))
        row[columns[vessel_id]] = 1.0
        # svFSI flux is positive outward. A vessel inlet uses the opposite sign;
        # a vessel outlet uses the same sign as the directed vessel flow.
        target = -face_flows[face] if side == "inlet" else face_flows[face]
        rows.append(row)
        rhs.append(target)

    for junction in model["junctions"]:
        row = np.zeros(len(vessel_ids))
        for vessel_id in junction.get("inlet_vessels", []):
            row[columns[int(vessel_id)]] += 1.0
        for vessel_id in junction.get("outlet_vessels", []):
            row[columns[int(vessel_id)]] -= 1.0
        rows.append(row)
        rhs.append(0.0)

    matrix = np.vstack(rows)
    targets = np.asarray(rhs)
    solution, _, rank, _ = np.linalg.lstsq(matrix, targets, rcond=None)
    if rank != len(vessel_ids):
        raise ValueError(f"Vessel-flow system rank is {rank}, expected {len(vessel_ids)}")
    residual = matrix @ solution - targets
    if np.max(np.abs(residual)) > 1e-6:
        raise ValueError(f"Vessel-flow continuity residual is {np.max(np.abs(residual)):.3e} mL/s")
    return {vessel_id: float(solution[columns[vessel_id]]) for vessel_id in vessel_ids}


def read_centerline(path: Path) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()
    if data.GetNumberOfPoints() == 0:
        raise ValueError(f"Centerline is empty: {path}")
    return data


def read_volume(path: Path) -> vtk.vtkUnstructuredGrid:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()
    if data.GetNumberOfPoints() == 0:
        raise ValueError(f"3D result is empty: {path}")
    return data


def probe_centerline_pressure(
    centerline: vtk.vtkPolyData, volume: vtk.vtkUnstructuredGrid
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    probe = vtk.vtkProbeFilter()
    probe.SetInputData(centerline)
    probe.SetSourceData(volume)
    probe.Update()
    output = probe.GetOutput()
    pressure = output.GetPointData().GetArray("Pressure")
    valid = output.GetPointData().GetArray("vtkValidPointMask")
    branch = centerline.GetPointData().GetArray("BranchId")
    path = centerline.GetPointData().GetArray("Path")
    if any(array is None for array in (pressure, valid, branch, path)):
        raise ValueError("Required Pressure, BranchId, Path, or validity array is missing")

    samples: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for point_id in range(centerline.GetNumberOfPoints()):
        if valid.GetTuple1(point_id) < 0.5:
            continue
        branch_id = int(branch.GetTuple1(point_id))
        if branch_id < 0:
            continue
        samples[branch_id].append(
            (float(path.GetTuple1(point_id)), float(pressure.GetTuple1(point_id)))
        )

    result = {}
    for branch_id, values in samples.items():
        values.sort()
        paths = np.asarray([item[0] for item in values])
        pressures = np.asarray([item[1] for item in values])
        unique_paths, unique_indices = np.unique(paths, return_index=True)
        result[branch_id] = (unique_paths, pressures[unique_indices])
    return result


def vessel_path_bounds(model: dict) -> dict[int, tuple[int, float, float]]:
    pattern = re.compile(r"branch(-?\d+)_seg(\d+)$")
    branches: dict[int, list[tuple[int, int, float]]] = defaultdict(list)
    for vessel in model["vessels"]:
        match = pattern.fullmatch(vessel["vessel_name"])
        if not match:
            raise ValueError(f"Unexpected vessel name {vessel['vessel_name']}")
        branches[int(match.group(1))].append(
            (int(match.group(2)), int(vessel["vessel_id"]), float(vessel["vessel_length"]))
        )

    bounds = {}
    for branch_id, segments in branches.items():
        cumulative = 0.0
        for _, vessel_id, length in sorted(segments):
            bounds[vessel_id] = (branch_id, cumulative, cumulative + length)
            cumulative += length
    return bounds


def vessel_endpoint_pressures(
    model: dict,
    centerline_pressures: dict[int, tuple[np.ndarray, np.ndarray]],
    boundary_pressures: dict[str, float],
) -> dict[int, tuple[float, float]]:
    bounds = vessel_path_bounds(model)
    attachments_by_endpoint = {
        (vessel_id, side): boundary_name
        for boundary_name, (vessel_id, side) in boundary_attachments(model).items()
    }
    result = {}
    for vessel in model["vessels"]:
        vessel_id = int(vessel["vessel_id"])
        branch_id, start, end = bounds[vessel_id]
        paths, pressures = centerline_pressures[branch_id]
        p_in = float(np.interp(start, paths, pressures))
        p_out = float(np.interp(end, paths, pressures))
        inlet_name = attachments_by_endpoint.get((vessel_id, "inlet"))
        outlet_name = attachments_by_endpoint.get((vessel_id, "outlet"))
        if inlet_name:
            p_in = boundary_pressures[boundary_face_name(inlet_name)]
        if outlet_name:
            p_out = boundary_pressures[boundary_face_name(outlet_name)]
        result[vessel_id] = (p_in, p_out)
    return result


def harmonize_normal_junction_pressures(
    model: dict,
    pressures: dict[int, tuple[float, float]],
    vessel_flows: dict[int, float],
) -> tuple[dict[int, tuple[float, float]], dict[str, float]]:
    """Project each 3D junction region onto one NORMAL_JUNCTION pressure."""
    adjusted = {vessel_id: list(values) for vessel_id, values in pressures.items()}
    junction_pressures = {}
    for junction in model["junctions"]:
        endpoint_values = []
        endpoint_weights = []
        for vessel_id in junction.get("inlet_vessels", []):
            vessel_id = int(vessel_id)
            endpoint_values.append(adjusted[vessel_id][1])
            endpoint_weights.append(abs(vessel_flows[vessel_id]))
        for vessel_id in junction.get("outlet_vessels", []):
            vessel_id = int(vessel_id)
            endpoint_values.append(adjusted[vessel_id][0])
            endpoint_weights.append(abs(vessel_flows[vessel_id]))
        if not endpoint_values:
            continue
        if sum(endpoint_weights) > 0.0:
            pressure = float(np.average(endpoint_values, weights=endpoint_weights))
        else:
            pressure = float(np.mean(endpoint_values))
        junction_pressures[junction["junction_name"]] = pressure
        for vessel_id in junction.get("inlet_vessels", []):
            adjusted[int(vessel_id)][1] = pressure
        for vessel_id in junction.get("outlet_vessels", []):
            adjusted[int(vessel_id)][0] = pressure
    return {
        vessel_id: (values[0], values[1]) for vessel_id, values in adjusted.items()
    }, junction_pressures


def connection_keys(model: dict) -> dict[int, tuple[tuple[str, str], tuple[str, str]]]:
    inlet_junction = {}
    outlet_junction = {}
    for junction in model["junctions"]:
        name = junction["junction_name"]
        for vessel_id in junction.get("outlet_vessels", []):
            inlet_junction[int(vessel_id)] = name
        for vessel_id in junction.get("inlet_vessels", []):
            outlet_junction[int(vessel_id)] = name

    keys = {}
    for vessel in model["vessels"]:
        vessel_id = int(vessel["vessel_id"])
        name = vessel["vessel_name"]
        boundary = vessel.get("boundary_conditions", {})
        if "inlet" in boundary:
            inlet_pair = (boundary["inlet"], name)
        else:
            inlet_pair = (inlet_junction[vessel_id], name)
        if "outlet" in boundary:
            outlet_pair = (name, boundary["outlet"])
        else:
            outlet_pair = (name, outlet_junction[vessel_id])
        keys[vessel_id] = (inlet_pair, outlet_pair)
    return keys


def make_observations(
    model: dict,
    vessel_flows: dict[int, float],
    pressures: dict[int, tuple[float, float]],
    repeats: int = 2,
) -> tuple[dict, list[dict]]:
    observations = {"y": {}, "dy": {}}
    rows = []
    for vessel in model["vessels"]:
        vessel_id = int(vessel["vessel_id"])
        vessel_name = vessel["vessel_name"]
        q = vessel_flows[vessel_id]
        p_in, p_out = pressures[vessel_id]
        inlet_pair, outlet_pair = connection_keys(model)[vessel_id]
        endpoint_data = (
            ("inlet", inlet_pair, p_in),
            ("outlet", outlet_pair, p_out),
        )
        for endpoint, pair, pressure in endpoint_data:
            pressure_key = f"pressure:{pair[0]}:{pair[1]}"
            flow_key = f"flow:{pair[0]}:{pair[1]}"
            observations["y"][pressure_key] = [pressure] * repeats
            observations["y"][flow_key] = [q] * repeats
            observations["dy"][pressure_key] = [0.0] * repeats
            observations["dy"][flow_key] = [0.0] * repeats
            rows.append(
                {
                    "vessel_id": vessel_id,
                    "vessel_name": vessel_name,
                    "endpoint": endpoint,
                    "connection": f"{pair[0]}:{pair[1]}",
                    "pressure_dyn_cm2": pressure,
                    "pressure_mmHg": pressure * MMHG_PER_DYN_CM2,
                    "flow_ml_s": q,
                    "pressure_source": "3D boundary average"
                    if pair[0].startswith(("FLOW_", "RESISTANCE_", "RCR_"))
                    or pair[1].startswith(("FLOW_", "RESISTANCE_", "RCR_"))
                    else "flow-weighted normal-junction pressure from 3D endpoint probes",
                }
            )
    return observations, rows


def find_executable(explicit: Path | None, name: str) -> Path:
    candidates = [
        explicit,
        Path(shutil.which(name)) if shutil.which(name) else None,
        Path.home() / "anaconda3" / "bin" / name,
    ]
    executable = next((item.expanduser() for item in candidates if item and item.is_file()), None)
    if executable is None:
        raise FileNotFoundError(f"{name} was not found")
    return executable.resolve()


def run_command(arguments: list[str], log_path: Path | None = None) -> None:
    print("Running:", " ".join(arguments))
    if log_path is None:
        completed = subprocess.run(arguments, check=False)
    else:
        with log_path.open("w") as stream:
            completed = subprocess.run(
                arguments,
                check=False,
                stdout=stream,
                stderr=subprocess.STDOUT,
                text=True,
            )
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {arguments[0]}")


def read_solver_result(path: Path) -> dict[str, dict[str, float]]:
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Solver result is empty: {path}")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["name"]].append(row)
    averages = {}
    for name, vessel_rows in grouped.items():
        averages[name] = {
            key: fmean(float(row[key]) for row in vessel_rows)
            for key in ("flow_in", "flow_out", "pressure_in", "pressure_out")
        }
    return averages


def metrics_from_solver(model: dict, result: Path) -> dict[str, float]:
    values = read_solver_result(result)
    attachments = boundary_attachments(model)
    vessels = {int(vessel["vessel_id"]): vessel for vessel in model["vessels"]}
    face_flow = {}
    face_pressure = {}
    for boundary_name, (vessel_id, side) in attachments.items():
        vessel_name = vessels[vessel_id]["vessel_name"]
        endpoint = values[vessel_name]
        flow = endpoint["flow_in" if side == "inlet" else "flow_out"]
        pressure = endpoint["pressure_in" if side == "inlet" else "pressure_out"]
        face = boundary_face_name(boundary_name)
        face_flow[face] = -flow if side == "inlet" else flow
        face_pressure[face] = pressure
    left = sum(face_flow[name] for name in LEFT_PULMONARY_FACES)
    right = sum(face_flow[name] for name in RIGHT_PULMONARY_FACES)
    total = left + right
    return {
        "fontan_pressure_mmHg": face_pressure["cap_fontan_x"] * MMHG_PER_DYN_CM2,
        "lpa_flow_ml_s": left,
        "rpa_flow_ml_s": right,
        "pulmonary_flow_ml_s": total,
        "lpa_fraction": left / total,
        "rpa_fraction": right / total,
    }


def metrics_from_3d(
    pressures: dict[str, float], flows: dict[str, float]
) -> dict[str, float]:
    left = sum(flows[name] for name in LEFT_PULMONARY_FACES)
    right = sum(flows[name] for name in RIGHT_PULMONARY_FACES)
    total = left + right
    return {
        "fontan_pressure_mmHg": pressures["cap_fontan_x"] * MMHG_PER_DYN_CM2,
        "lpa_flow_ml_s": left,
        "rpa_flow_ml_s": right,
        "pulmonary_flow_ml_s": total,
        "lpa_fraction": left / total,
        "rpa_fraction": right / total,
    }


def apply_fontan_pressure_correction(
    model: dict,
    target_pressure_mmHg: float,
    current_pressure_mmHg: float,
    face_flows: dict[str, float],
) -> dict[str, float | str]:
    """Correct the prescribed-flow root resistance without changing flow split."""
    attachments = boundary_attachments(model)
    boundary_name = next(
        name for name in attachments if boundary_face_name(name) == "cap_fontan_x"
    )
    vessel_id, side = attachments[boundary_name]
    vessels = {int(vessel["vessel_id"]): vessel for vessel in model["vessels"]}
    vessel = vessels[vessel_id]
    directed_flow = -face_flows["cap_fontan_x"] if side == "inlet" else face_flows["cap_fontan_x"]
    if abs(directed_flow) < 1.0e-12:
        raise ValueError("Cannot correct Fontan pressure with zero root flow")
    pressure_difference = (target_pressure_mmHg - current_pressure_mmHg) / MMHG_PER_DYN_CM2
    resistance_correction = pressure_difference / directed_flow
    values = vessel["zero_d_element_values"]
    initial_resistance = float(values["R_poiseuille"])
    corrected_resistance = initial_resistance + resistance_correction
    if corrected_resistance <= 0.0:
        raise ValueError("Fontan pressure correction would create nonpositive resistance")
    values["R_poiseuille"] = corrected_resistance
    return {
        "vessel_name": vessel["vessel_name"],
        "resistance_before_dyn_s_cm5": initial_resistance,
        "resistance_correction_dyn_s_cm5": resistance_correction,
        "resistance_after_dyn_s_cm5": corrected_resistance,
    }


def clean_calibrated_model(raw: dict, baseline: dict) -> tuple[dict, list[str]]:
    """Make calibrator output simulation-ready and restore steady-unidentifiable terms."""
    calibrated = copy.deepcopy(raw)
    calibrated.pop("y", None)
    calibrated.pop("dy", None)
    calibrated.pop("calibration_parameters", None)
    calibrated["boundary_conditions"] = copy.deepcopy(baseline["boundary_conditions"])
    calibrated["simulation_parameters"] = copy.deepcopy(baseline["simulation_parameters"])
    calibrated["simulation_parameters"]["model_name"] = "P15_steady_3D_calibrated"

    baseline_vessels = {int(vessel["vessel_id"]): vessel for vessel in baseline["vessels"]}
    restored_resistances = []
    for vessel in calibrated["vessels"]:
        original = baseline_vessels[int(vessel["vessel_id"])]["zero_d_element_values"]
        values = vessel["zero_d_element_values"]
        if not math.isfinite(float(values["R_poiseuille"])) or float(values["R_poiseuille"]) <= 0.0:
            values["R_poiseuille"] = original["R_poiseuille"]
            restored_resistances.append(vessel["vessel_name"])
        for parameter in ("C", "L", "stenosis_coefficient"):
            if parameter in original:
                values[parameter] = original[parameter]
    return calibrated, restored_resistances


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parameter_rows(baseline: dict, raw: dict, calibrated: dict) -> list[dict]:
    raw_by_id = {int(vessel["vessel_id"]): vessel for vessel in raw["vessels"]}
    calibrated_by_id = {int(vessel["vessel_id"]): vessel for vessel in calibrated["vessels"]}
    rows = []
    for vessel in baseline["vessels"]:
        vessel_id = int(vessel["vessel_id"])
        initial = vessel["zero_d_element_values"]
        raw_values = raw_by_id[vessel_id]["zero_d_element_values"]
        final = calibrated_by_id[vessel_id]["zero_d_element_values"]
        r_initial = float(initial["R_poiseuille"])
        r_raw = float(raw_values["R_poiseuille"])
        r_final = float(final["R_poiseuille"])
        rows.append(
            {
                "vessel_id": vessel_id,
                "vessel_name": vessel["vessel_name"],
                "R_initial_dyn_s_cm5": r_initial,
                "R_calibrator_raw_dyn_s_cm5": r_raw,
                "R_applied_dyn_s_cm5": r_final,
                "R_ratio": r_final / r_initial if r_initial else math.nan,
                "raw_resistance_restored": r_raw <= 0.0 or not math.isfinite(r_raw),
                "C_preserved_cm5_dyn": final.get("C", math.nan),
                "L_preserved_dyn_s2_cm5": final.get("L", math.nan),
                "stenosis_preserved_dyn_s2_cm8": final.get("stenosis_coefficient", math.nan),
            }
        )
    return rows


def comparison_rows(target: dict, baseline: dict, calibrated: dict) -> list[dict]:
    definitions = (
        ("Mean Fontan pressure", "fontan_pressure_mmHg", "mmHg"),
        ("LPA flow", "lpa_flow_ml_s", "mL/s"),
        ("RPA flow", "rpa_flow_ml_s", "mL/s"),
        ("Total pulmonary flow", "pulmonary_flow_ml_s", "mL/s"),
        ("LPA split", "lpa_fraction", "%"),
        ("RPA split", "rpa_fraction", "%"),
    )
    rows = []
    for metric, key, unit in definitions:
        scale = 100.0 if unit == "%" else 1.0
        target_value = target[key] * scale
        baseline_value = baseline[key] * scale
        calibrated_value = calibrated[key] * scale
        rows.append(
            {
                "metric": metric,
                "unit": unit,
                "3D_target": target_value,
                "0D_baseline": baseline_value,
                "baseline_relative_error": abs(baseline_value - target_value) / abs(target_value),
                "0D_calibrated": calibrated_value,
                "calibrated_relative_error": abs(calibrated_value - target_value) / abs(target_value),
            }
        )
    return rows


def plot_summary(comparison: list[dict], parameters: list[dict], path: Path) -> None:
    width, height = 1800, 620
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    regular_path = "/System/Library/Fonts/Supplemental/Arial.ttf"
    bold_path = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
    regular = ImageFont.truetype(regular_path, 20)
    small = ImageFont.truetype(regular_path, 17)
    bold = ImageFont.truetype(bold_path, 25)
    colors = ("#343A40", "#3B82C4", "#D1495B")

    def categorical_bars(box, title, ylabel, labels, values):
        left, top, right, bottom = box
        plot_left, plot_top = left + 75, top + 56
        plot_right, plot_bottom = right - 25, bottom - 78
        maximum = max(values) * 1.14
        draw.text((left + 18, top + 10), title, fill="#1D2329", font=bold)
        draw.line((plot_left, plot_top, plot_left, plot_bottom), fill="#555555", width=2)
        draw.line((plot_left, plot_bottom, plot_right, plot_bottom), fill="#555555", width=2)
        draw.text((left + 5, plot_top), ylabel, fill="#333333", font=small)
        spacing = (plot_right - plot_left) / len(values)
        bar_width = spacing * 0.52
        for index, (label, value, color) in enumerate(zip(labels, values, colors)):
            center = plot_left + spacing * (index + 0.5)
            y = plot_bottom - (value / maximum) * (plot_bottom - plot_top)
            draw.rectangle((center - bar_width / 2, y, center + bar_width / 2, plot_bottom), fill=color)
            value_text = f"{value:.3f}"
            value_box = draw.textbbox((0, 0), value_text, font=regular)
            draw.text((center - (value_box[2] - value_box[0]) / 2, y - 29), value_text, fill="#222222", font=regular)
            label_box = draw.textbbox((0, 0), label, font=small)
            draw.text((center - (label_box[2] - label_box[0]) / 2, plot_bottom + 14), label, fill="#333333", font=small)

    pressure = comparison[0]
    categorical_bars(
        (20, 20, 590, 600),
        "Mean Fontan pressure",
        "mmHg",
        ("3D target", "0D baseline", "0D calibrated"),
        (pressure["3D_target"], pressure["0D_baseline"], pressure["0D_calibrated"]),
    )
    lpa = comparison[4]
    categorical_bars(
        (610, 20, 1180, 600),
        "Pulmonary flow split",
        "LPA fraction [%]",
        ("3D target", "0D baseline", "0D calibrated"),
        (lpa["3D_target"], lpa["0D_baseline"], lpa["0D_calibrated"]),
    )

    left, top, right, bottom = 1200, 20, 1780, 600
    plot_left, plot_top = left + 75, top + 56
    plot_right, plot_bottom = right - 25, bottom - 70
    ratios = [float(row["R_ratio"]) for row in parameters]
    low = min(0.0, min(ratios))
    high = max(1.0, max(ratios))
    span = max(high - low, 1e-12)
    zero_y = plot_bottom - (0.0 - low) / span * (plot_bottom - plot_top)
    one_y = plot_bottom - (1.0 - low) / span * (plot_bottom - plot_top)
    draw.text((left + 18, top + 10), "Branch resistance update", fill="#1D2329", font=bold)
    draw.text((left + 5, plot_top), "R calibrated / initial", fill="#333333", font=small)
    draw.line((plot_left, plot_top, plot_left, plot_bottom), fill="#555555", width=2)
    draw.line((plot_left, zero_y, plot_right, zero_y), fill="#555555", width=2)
    draw.line((plot_left, one_y, plot_right, one_y), fill="#222222", width=2)
    spacing = (plot_right - plot_left) / len(ratios)
    for index, ratio in enumerate(ratios):
        center = plot_left + spacing * (index + 0.5)
        ratio_y = plot_bottom - (ratio - low) / span * (plot_bottom - plot_top)
        draw.rectangle(
            (center - spacing * 0.36, min(zero_y, ratio_y), center + spacing * 0.36, max(zero_y, ratio_y)),
            fill="#4C956C" if ratio >= 0.0 else "#D1495B",
        )
        if index % 4 == 0:
            draw.text((center - 7, plot_bottom + 13), str(index), fill="#333333", font=small)
    draw.text((plot_left, plot_bottom + 39), "Vessel ID", fill="#333333", font=small)
    image.save(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--baseline-result", type=Path, default=DEFAULT_BASELINE_RESULT)
    parser.add_argument("--centerline", type=Path, default=DEFAULT_CENTERLINE)
    parser.add_argument("--three-d-result", type=Path, default=DEFAULT_3D_RESULT)
    parser.add_argument("--three-d-integrals", type=Path, default=DEFAULT_3D_INTEGRALS)
    parser.add_argument("--last", type=int, default=100)
    parser.add_argument("--calibrator", type=Path)
    parser.add_argument("--solver", type=Path)
    args = parser.parse_args()
    if args.last <= 0:
        raise ValueError("--last must be positive")

    model = read_json(args.model.resolve())
    boundary_pressures, boundary_flows = final_window_boundary_means(
        args.three_d_integrals.resolve(), args.last
    )
    vessel_flows = solve_steady_vessel_flows(model, boundary_flows)
    centerline = read_centerline(args.centerline.resolve())
    volume = read_volume(args.three_d_result.resolve())
    centerline_pressures = probe_centerline_pressure(centerline, volume)
    endpoint_pressures = vessel_endpoint_pressures(
        model, centerline_pressures, boundary_pressures
    )
    endpoint_pressures, junction_pressures = harmonize_normal_junction_pressures(
        model, endpoint_pressures, vessel_flows
    )
    observations, observation_rows = make_observations(
        model, vessel_flows, endpoint_pressures
    )

    calibration_input = copy.deepcopy(model)
    calibration_input["simulation_parameters"]["number_of_cardiac_cycles"] = 1
    calibration_input["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"] = 2
    calibration_input["simulation_parameters"]["output_all_cycles"] = True
    calibration_input["calibration_parameters"] = {
        "tolerance_gradient": 1.0e-6,
        "tolerance_increment": 1.0e-2,
        "maximum_iterations": 2000,
        "calibrate_stenosis_coefficient": False,
        "set_capacitance_to_zero": False,
        "initial_damping_factor": 1.0,
    }
    calibration_input.update(observations)

    input_path = HERE / "P15_steady_calibration_input.json"
    raw_path = HERE / "P15_0d_calibrated_raw.json"
    calibrated_path = HERE / "P15_0d_calibrated.json"
    calibrated_result_path = HERE / "P15_0d_calibrated.csv"
    observations_path = HERE / "calibration_observations.csv"
    parameters_path = HERE / "vessel_parameter_comparison.csv"
    comparison_path = HERE / "calibration_metrics.csv"
    summary_path = HERE / "calibration_summary.json"
    figure_path = HERE / "calibration_summary.png"
    calibrator_log_path = HERE / "svzerodcalibrator.log"

    write_json(input_path, calibration_input)
    write_csv(observations_path, observation_rows)

    calibrator = find_executable(args.calibrator, "svzerodcalibrator")
    solver = find_executable(args.solver, "svzerodsolver")
    run_command(
        [str(calibrator), str(input_path), str(raw_path)],
        log_path=calibrator_log_path,
    )

    raw_calibrated = read_json(raw_path)
    calibrated, restored_resistances = clean_calibrated_model(raw_calibrated, model)
    write_json(calibrated_path, calibrated)
    run_command([str(solver), str(calibrated_path), str(calibrated_result_path)])

    target_metrics = metrics_from_3d(boundary_pressures, boundary_flows)
    baseline_metrics = metrics_from_solver(model, args.baseline_result.resolve())
    provisional_metrics = metrics_from_solver(calibrated, calibrated_result_path)
    fontan_correction = apply_fontan_pressure_correction(
        calibrated,
        target_metrics["fontan_pressure_mmHg"],
        provisional_metrics["fontan_pressure_mmHg"],
        boundary_flows,
    )
    write_json(calibrated_path, calibrated)
    run_command([str(solver), str(calibrated_path), str(calibrated_result_path)])
    calibrated_metrics = metrics_from_solver(calibrated, calibrated_result_path)
    parameters = parameter_rows(model, raw_calibrated, calibrated)
    comparison = comparison_rows(target_metrics, baseline_metrics, calibrated_metrics)
    write_csv(parameters_path, parameters)
    write_csv(comparison_path, comparison)
    plot_summary(comparison, parameters, figure_path)

    summary = {
        "method": "svZeroDCalibrator steady branch-resistance calibration",
        "calibrator": str(calibrator),
        "solver": str(solver),
        "calibrator_termination": calibrator_log_path.read_text().splitlines()[-1],
        "three_d_window": {
            "number_of_samples": min(args.last, len(read_boundary_integral(args.three_d_integrals.resolve() / "B_NS_Velocity_flux.txt"))),
            "source": str(args.three_d_integrals.resolve()),
        },
        "identifiability": {
            "calibrated": ["R_poiseuille"],
            "preserved": ["C", "L", "stenosis_coefficient"],
            "reason": "Steady data have zero temporal derivatives, so C and L are not identifiable.",
        },
        "normal_junction_pressures_mmHg": {
            name: value * MMHG_PER_DYN_CM2
            for name, value in junction_pressures.items()
        },
        "fontan_root_pressure_correction": fontan_correction,
        "three_d_target": target_metrics,
        "zero_d_baseline": baseline_metrics,
        "zero_d_calibrated": calibrated_metrics,
        "maximum_calibrated_target_relative_error": max(
            row["calibrated_relative_error"] for row in comparison
        ),
        "raw_negative_calibrated_resistance_count": len(restored_resistances),
        "raw_negative_resistances_restored_to_geometric_values": restored_resistances,
        "applied_negative_resistance_count": sum(
            row["R_applied_dyn_s_cm5"] < 0.0 for row in parameters
        ),
        "outputs": {
            "calibration_input": str(input_path),
            "calibrated_model": str(calibrated_path),
            "calibrated_result": str(calibrated_result_path),
            "metric_table": str(comparison_path),
            "parameter_table": str(parameters_path),
            "observation_table": str(observations_path),
            "summary_figure": str(figure_path),
            "calibrator_log": str(calibrator_log_path),
        },
    }
    write_json(summary_path, summary)

    print("\nCalibration metrics")
    for row in comparison:
        print(
            f"  {row['metric']:<25} target={row['3D_target']:.6g} "
            f"baseline={row['0D_baseline']:.6g} calibrated={row['0D_calibrated']:.6g} "
            f"error={100.0 * row['calibrated_relative_error']:.3f}%"
        )
    print(f"\nWrote calibrated model: {calibrated_path}")
    print(f"Wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
