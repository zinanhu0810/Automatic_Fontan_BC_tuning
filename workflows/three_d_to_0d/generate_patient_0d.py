#!/usr/bin/env python3
"""Generate a multi-inlet Fontan svZeroD model from a 3D CFD case.

The geometry-to-ROM step delegates centerline extraction and vessel parameter
generation to SimVascular's ``sv_rom_simulation`` package. This driver supplies
the GUI parameters from files, then adapts the package's single-inlet tree to
the multiple prescribed-flow inlets used by a Fontan model.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import json
import math
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from statistics import fmean
from typing import Any


DEFAULT_SIMVASCULAR = Path(
    "/Applications/SimVascular.app/Contents/Resources/simvascular"
)
FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a patient-specific 0D ROM from a SimVascular 3D case."
    )
    parser.add_argument("patient_dir", type=Path, help="Patient case directory")
    parser.add_argument(
        "--config",
        type=Path,
        help="Patient config JSON (default: PATIENT_DIR/patient_0d_config.json)",
    )
    parser.add_argument(
        "--simvascular",
        type=Path,
        default=DEFAULT_SIMVASCULAR,
        help="SimVascular launcher script",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare ROM inputs without launching SimVascular",
    )
    parser.add_argument(
        "--reuse-centerlines",
        action="store_true",
        help="Reuse OUTPUT_DIR/centerlines.vtp instead of recomputing centerlines",
    )
    parser.add_argument(
        "--run-baseline",
        action="store_true",
        help="Run the generated JSON with svZeroDSolver and check mass conservation",
    )
    parser.add_argument(
        "--svzerodsolver",
        type=Path,
        help="svZeroDSolver executable used with --run-baseline",
    )
    return parser.parse_args()


def as_float(text: str | None, field: str) -> float:
    if text is None:
        raise ValueError(f"Missing numeric value for {field}")
    match = FLOAT_PATTERN.search(text.strip())
    if not match:
        raise ValueError(f"Invalid numeric value for {field}: {text!r}")
    return float(match.group(0))


def text_at(parent: ET.Element, path: str, field: str) -> str:
    node = parent.find(path)
    if node is None or node.text is None:
        raise ValueError(f"Missing XML field {field} at {path}")
    return node.text.strip()


def resolve(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def read_config(patient_dir: Path, config_path: Path) -> dict[str, Any]:
    with config_path.open() as stream:
        config = json.load(stream)
    required = [
        "patient_id",
        "solver_xml",
        "surface_model",
        "boundary_surfaces_directory",
        "centerline_root",
        "output_directory",
    ]
    missing = [key for key in required if not config.get(key)]
    if missing:
        raise ValueError(f"Missing config fields: {', '.join(missing)}")
    config["_config_path"] = str(config_path)
    config["_patient_dir"] = str(patient_dir)
    return config


def find_fluid_equation(root: ET.Element) -> ET.Element:
    for equation in root.findall("Add_equation"):
        if equation.attrib.get("type", "").lower() == "fluid":
            return equation
    raise ValueError("No fluid Add_equation was found in solver XML")


def parse_solver_tree(solver_xml: Path) -> ET.Element:
    xml_text = solver_xml.read_text(errors="replace")
    # Some legacy svMultiPhysics files contain this exact opening/closing tag
    # mismatch. Normalize it in memory so the source simulation file is untouched.
    xml_text = xml_text.replace(
        "</Number_of_time_isteps>", "</Number_of_time_steps>"
    )
    try:
        return ET.fromstring(xml_text)
    except ET.ParseError as error:
        raise ValueError(f"Cannot parse solver XML {solver_xml}: {error}") from error


def parse_solver_xml(solver_xml: Path, patient_dir: Path) -> dict[str, Any]:
    root = parse_solver_tree(solver_xml)
    general = root.find("GeneralSimulationParameters")
    if general is None:
        raise ValueError("No GeneralSimulationParameters found in solver XML")
    fluid = find_fluid_equation(root)

    number_steps_node = general.find("Number_of_time_steps")
    if number_steps_node is None:
        number_steps_node = general.find("Number_of_time_isteps")
    if number_steps_node is None or number_steps_node.text is None:
        raise ValueError("No number of time steps found in solver XML")

    data: dict[str, Any] = {
        "time_step": as_float(
            text_at(general, "Time_step_size", "time step"), "time step"
        ),
        "number_of_time_steps": int(as_float(number_steps_node.text, "time steps")),
        "density": as_float(text_at(fluid, "Density", "density"), "density"),
        "viscosity": as_float(
            text_at(fluid, "Viscosity/Value", "viscosity"), "viscosity"
        ),
        "flow_boundaries": {},
        "rcr_boundaries": {},
        "wall_boundaries": [],
    }

    for bc in fluid.findall("Add_BC"):
        name = bc.attrib.get("name")
        if not name:
            continue
        bc_type = (bc.findtext("Type") or "").strip().lower()
        time_dependence = (bc.findtext("Time_dependence") or "").strip().lower()
        if bc_type == "dir" and time_dependence == "unsteady":
            waveform = text_at(
                bc, "Temporal_values_file_path", f"waveform for {name}"
            )
            data["flow_boundaries"][name] = {
                "waveform": str(resolve(patient_dir, waveform))
            }
        elif bc_type == "neu" and time_dependence == "rcr":
            values = bc.find("RCR_values")
            if values is None:
                raise ValueError(f"RCR boundary {name} has no RCR_values")
            data["rcr_boundaries"][name] = {
                "Rp": as_float(values.findtext("Proximal_resistance"), f"{name} Rp"),
                "C": as_float(values.findtext("Capacitance"), f"{name} C"),
                "Rd": as_float(values.findtext("Distal_resistance"), f"{name} Rd"),
                "Pd": as_float(values.findtext("Distal_pressure"), f"{name} Pd"),
            }
        elif bc_type == "dir" and time_dependence == "steady":
            data["wall_boundaries"].append(name)

    if not data["flow_boundaries"]:
        raise ValueError("No prescribed-flow boundaries found in fluid equation")
    if not data["rcr_boundaries"]:
        raise ValueError("No RCR boundaries found in fluid equation")
    return data


def read_flow_file(path: Path) -> tuple[list[float], list[float]]:
    rows: list[tuple[float, float]] = []
    with path.open() as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            fields = re.split(r"[\s,]+", stripped)
            if len(fields) < 2:
                continue
            try:
                time = float(fields[0])
                flow = float(fields[1])
            except ValueError:
                continue
            if not rows and time > 2.0 and flow > 2.0:
                # SimVascular .flow files may start with a point/mode count header.
                continue
            rows.append((time, flow))
    if len(rows) < 2:
        raise ValueError(f"Flow file has fewer than two data rows: {path}")
    times, flows = zip(*rows)
    if any(t1 <= t0 for t0, t1 in zip(times, times[1:])):
        raise ValueError(f"Flow times must be strictly increasing: {path}")
    return list(times), list(flows)


def write_root_inflow(path: Path, times: list[float], native_flows: list[float]) -> None:
    with path.open("w") as stream:
        for time, flow in zip(times, native_flows):
            # The root vessel points into the network; SimVascular's surface flux is
            # negative for inflow, so reverse it at the root inlet.
            stream.write(f"{time:.10g} {-flow:.12g}\n")


def write_outlet_names(path: Path, names: list[str]) -> None:
    with path.open("w") as stream:
        for name in names:
            stream.write(f"{name}\n")


def write_rcrt(
    path: Path,
    outlet_names: list[str],
    rcr_boundaries: dict[str, dict[str, float]],
) -> None:
    # The ROM package requires one supported terminal BC per target centerline.
    # Flow terminals receive a harmless placeholder that is replaced afterward.
    with path.open("w") as stream:
        stream.write("2\n")
        for name in outlet_names:
            values = rcr_boundaries.get(
                name, {"Rp": 1.0, "C": 1.0e-8, "Rd": 1.0, "Pd": 0.0}
            )
            stream.write("2\n")
            stream.write(f"{name}\n")
            stream.write(f"{values['Rp']:.16g}\n")
            stream.write(f"{values['C']:.16g}\n")
            stream.write(f"{values['Rd']:.16g}\n")
            stream.write(f"0.0 {values['Pd']:.16g}\n")
            stream.write(f"1.0 {values['Pd']:.16g}\n")


def cap_surface_names(boundary_dir: Path) -> set[str]:
    return {
        path.stem
        for path in boundary_dir.glob("*.vtp")
        if not path.name.lower().startswith("wall")
    }


def prepare_inputs(
    patient_dir: Path,
    config: dict[str, Any],
    cfd: dict[str, Any],
) -> dict[str, Any]:
    output_dir = resolve(patient_dir, config["output_directory"])
    input_dir = output_dir / "rom_inputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)

    root_name = config["centerline_root"]
    flow_boundaries = cfd["flow_boundaries"]
    rcr_boundaries = cfd["rcr_boundaries"]
    if root_name not in flow_boundaries:
        raise ValueError(f"Centerline root {root_name!r} is not a FLOW boundary")

    boundary_dir = resolve(patient_dir, config["boundary_surfaces_directory"])
    surface_caps = cap_surface_names(boundary_dir)
    xml_caps = set(flow_boundaries) | set(rcr_boundaries)
    if surface_caps != xml_caps:
        missing_xml = sorted(surface_caps - xml_caps)
        missing_surface = sorted(xml_caps - surface_caps)
        raise ValueError(
            "Cap/BC mismatch. Caps without fluid BC: "
            f"{missing_xml}; fluid BCs without cap surface: {missing_surface}"
        )

    outlet_names = sorted(xml_caps - {root_name})
    root_times, root_native_flows = read_flow_file(
        Path(flow_boundaries[root_name]["waveform"])
    )
    root_inflow = input_dir / "root_inflow.dat"
    outlet_file = input_dir / "centerlines_outlets.dat"
    rcrt_file = input_dir / "rcrt.dat"
    write_root_inflow(root_inflow, root_times, root_native_flows)
    write_outlet_names(outlet_file, outlet_names)
    write_rcrt(rcrt_file, outlet_names, rcr_boundaries)

    raw_json = output_dir / f"{config['patient_id']}_0d_raw.json"
    centerlines = output_dir / "centerlines.vtp"
    runtime_config = output_dir / "rom_runtime_config.json"
    rom_options = config.get("rom_options", {})
    rom_arguments: dict[str, Any] = {
        "model_name": config["patient_id"],
        "model_order": 0,
        "boundary_surfaces_directory": str(boundary_dir),
        "centerlines_output_file": str(centerlines),
        "compute_centerlines": True,
        "density": cfd["density"],
        "inlet_face_input_file": f"{root_name}.vtp",
        "inflow_input_file": str(root_inflow),
        "material_model": rom_options.get("material_model", "OLUFSEN"),
        "num_time_steps": cfd["number_of_time_steps"],
        "outflow_bc_input_file": str(input_dir),
        "outflow_bc_type": "rcrt.dat",
        "outlet_face_names_input_file": str(outlet_file),
        "output_directory": str(output_dir),
        "solver_output_file": raw_json.name,
        "surface_model": str(resolve(patient_dir, config["surface_model"])),
        "time_step": cfd["time_step"],
        "uniform_bc": "false",
        "units": config.get("units", "cm"),
        "viscosity": cfd["viscosity"],
    }
    option_map = {
        "element_size": "element_size",
        "min_num_elements": "min_num_elements",
        "seg_min_num": "seg_min_num",
        "seg_size": "seg_size",
        "linear_material_ehr": "linear_material_ehr",
        "linear_material_pressure": "linear_material_pressure",
        "olufsen_material_k1": "olufsen_material_k1",
        "olufsen_material_k2": "olufsen_material_k2",
        "olufsen_material_k3": "olufsen_material_k3",
        "olufsen_material_exp": "olufsen_material_exp",
        "olufsen_material_pressure": "olufsen_material_pressure",
    }
    for config_key, argument_key in option_map.items():
        if config_key in rom_options:
            rom_arguments[argument_key] = rom_options[config_key]

    runtime = {
        "source": "SimVascular/Python/site-packages/sv_rom_simulation",
        "rom_arguments": rom_arguments,
        "expected_output_file": str(raw_json),
    }
    with runtime_config.open("w") as stream:
        json.dump(runtime, stream, indent=2)
        stream.write("\n")

    return {
        "output_dir": output_dir,
        "input_dir": input_dir,
        "raw_json": raw_json,
        "centerlines": centerlines,
        "runtime_config": runtime_config,
        "outlet_names": outlet_names,
        "root_name": root_name,
        "root_times": root_times,
        "root_native_flows": root_native_flows,
        "surface_caps": surface_caps,
    }


def run_simvascular(
    launcher: Path,
    runtime_config: Path,
    reuse_centerlines: bool,
) -> None:
    if not launcher.is_file():
        raise FileNotFoundError(f"SimVascular launcher not found: {launcher}")
    entry_point = Path(__file__).with_name("sv_rom_generate.py").resolve()
    if reuse_centerlines:
        with runtime_config.open() as stream:
            runtime = json.load(stream)
        centerlines = Path(runtime["rom_arguments"]["centerlines_output_file"])
        if not centerlines.is_file():
            raise FileNotFoundError(f"Cannot reuse missing centerlines: {centerlines}")
        runtime["rom_arguments"]["compute_centerlines"] = False
        runtime["rom_arguments"]["centerlines_input_file"] = str(centerlines)
        with runtime_config.open("w") as stream:
            json.dump(runtime, stream, indent=2)
            stream.write("\n")

    command = [str(launcher), "--python", "--", str(entry_point), str(runtime_config)]
    print("Running:", " ".join(command))
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"SimVascular ROM generation failed with exit code {completed.returncode}"
        )


def replace_boundary_name(model: dict[str, Any], old: str, new: str) -> None:
    replacements = 0
    for vessel in model.get("vessels", []):
        for location, value in vessel.get("boundary_conditions", {}).items():
            if value == old:
                vessel["boundary_conditions"][location] = new
                replacements += 1
    if replacements != 1:
        raise ValueError(
            f"Expected one vessel attachment for {old}, found {replacements}"
        )


def adapt_multi_inlet_model(
    patient_dir: Path,
    config: dict[str, Any],
    cfd: dict[str, Any],
    prepared: dict[str, Any],
) -> tuple[Path, dict[str, Any]]:
    raw_json = prepared["raw_json"]
    with raw_json.open() as stream:
        model = json.load(stream)

    # SimVascular 2023 used the legacy name ``internal_junction`` for serial
    # connections. Current svZeroDSolver and the upstream ROM writer use
    # NORMAL_JUNCTION for both serial and branching vessel connections.
    for junction in model.get("junctions", []):
        if junction.get("junction_type", "").lower() == "internal_junction":
            junction["junction_type"] = "NORMAL_JUNCTION"

    generated_outlet_file = prepared["output_dir"] / "centerlines_outlets.dat"
    if not generated_outlet_file.is_file():
        raise FileNotFoundError(
            "SimVascular did not write centerlines_outlets.dat; terminal mapping is unavailable"
        )
    generated_outlet_names = [
        line.strip()
        for line in generated_outlet_file.read_text().splitlines()
        if line.strip()
    ]

    boundary_by_name = {
        bc["bc_name"]: bc for bc in model.get("boundary_conditions", [])
    }
    root_name = prepared["root_name"]
    root_bc_name = f"FLOW_{root_name}"
    root_bc = boundary_by_name.get("INFLOW")
    if root_bc is None:
        raise ValueError("Raw ROM has no INFLOW boundary condition")
    root_bc["bc_name"] = root_bc_name
    root_bc["bc_values"] = {
        "t": prepared["root_times"],
        "Q": [-value for value in prepared["root_native_flows"]],
    }
    replace_boundary_name(model, "INFLOW", root_bc_name)

    for index, cap_name in enumerate(generated_outlet_names):
        old_name = f"RCR_{index}"
        bc = boundary_by_name.get(old_name)
        if bc is None:
            raise ValueError(f"Raw ROM has no boundary condition {old_name} for {cap_name}")
        if cap_name in cfd["flow_boundaries"]:
            times, native_flows = read_flow_file(
                Path(cfd["flow_boundaries"][cap_name]["waveform"])
            )
            new_name = f"FLOW_{cap_name}"
            bc.clear()
            bc.update(
                {
                    "bc_name": new_name,
                    "bc_type": "FLOW",
                    # This BC is attached at a terminal vessel outlet. Native CFD
                    # flux is negative for flow entering the 3D domain, which is
                    # also the correct signed flow relative to the rooted tree.
                    "bc_values": {"t": times, "Q": native_flows},
                }
            )
        elif cap_name in cfd["rcr_boundaries"]:
            values = cfd["rcr_boundaries"][cap_name]
            new_name = f"RCR_{cap_name}"
            bc.clear()
            bc.update(
                {
                    "bc_name": new_name,
                    "bc_type": "RCR",
                    "bc_values": dict(values),
                }
            )
        else:
            raise ValueError(f"Generated terminal {cap_name} has no CFD boundary data")
        replace_boundary_name(model, old_name, new_name)

    output_json = prepared["output_dir"] / f"{config['patient_id']}_0d.json"
    with output_json.open("w") as stream:
        json.dump(model, stream, indent=2, sort_keys=False)
        stream.write("\n")
    return output_json, model


def trapz_mean(times: list[float], values: list[float]) -> float:
    duration = times[-1] - times[0]
    if duration <= 0.0:
        raise ValueError("Waveform duration must be positive")
    integral = sum(
        0.5 * (q0 + q1) * (t1 - t0)
        for t0, t1, q0, q1 in zip(times, times[1:], values, values[1:])
    )
    return integral / duration


def validate_model(
    model: dict[str, Any],
    cfd: dict[str, Any],
    prepared: dict[str, Any],
) -> dict[str, Any]:
    vessels = model.get("vessels", [])
    junctions = model.get("junctions", [])
    boundaries = model.get("boundary_conditions", [])
    boundary_names = {bc["bc_name"] for bc in boundaries}
    attached: list[str] = []
    for vessel in vessels:
        attached.extend(vessel.get("boundary_conditions", {}).values())
    if set(attached) != boundary_names:
        raise ValueError(
            "Boundary attachment mismatch: unattached or undefined names are present"
        )
    if len(attached) != len(boundary_names):
        raise ValueError("A boundary condition is attached to more than one vessel")
    unsupported_junctions = sorted(
        {
            junction.get("junction_type")
            for junction in junctions
            if junction.get("junction_type") != "NORMAL_JUNCTION"
        }
    )
    if unsupported_junctions:
        raise ValueError(f"Unsupported junction types: {unsupported_junctions}")

    expected_flow_names = {f"FLOW_{name}" for name in cfd["flow_boundaries"]}
    expected_rcr_names = {f"RCR_{name}" for name in cfd["rcr_boundaries"]}
    if boundary_names != expected_flow_names | expected_rcr_names:
        raise ValueError("Generated boundary names do not match CFD cap names")
    boundary_map = {bc["bc_name"]: bc for bc in boundaries}
    for cap_name, expected_values in cfd["rcr_boundaries"].items():
        actual = boundary_map[f"RCR_{cap_name}"]
        if actual["bc_type"] != "RCR" or actual["bc_values"] != expected_values:
            raise ValueError(f"RCR values changed during generation for {cap_name}")
    for cap_name, waveform in cfd["flow_boundaries"].items():
        times, native = read_flow_file(Path(waveform["waveform"]))
        expected_flow = (
            [-value for value in native]
            if cap_name == prepared["root_name"]
            else native
        )
        actual = boundary_map[f"FLOW_{cap_name}"]
        if (
            actual["bc_type"] != "FLOW"
            or actual["bc_values"]["t"] != times
            or actual["bc_values"]["Q"] != expected_flow
        ):
            raise ValueError(f"FLOW waveform changed during generation for {cap_name}")

    flow_means: dict[str, float] = {}
    for name, values in cfd["flow_boundaries"].items():
        times, native = read_flow_file(Path(values["waveform"]))
        flow_means[name] = trapz_mean(times, native)

    vessel_values = [vessel.get("zero_d_element_values", {}) for vessel in vessels]
    for index, values in enumerate(vessel_values):
        for key in ("R_poiseuille", "L", "C", "stenosis_coefficient"):
            value = values.get(key)
            if value is None or not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"Invalid {key} for vessel index {index}: {value}")

    return {
        "status": "PASS",
        "number_of_vessels": len(vessels),
        "number_of_junctions": len(junctions),
        "number_of_boundaries": len(boundaries),
        "number_of_flow_inlets": len(cfd["flow_boundaries"]),
        "number_of_rcr_outlets": len(cfd["rcr_boundaries"]),
        "exact_cfd_boundary_values": True,
        "centerline_root": prepared["root_name"],
        "prescribed_inlet_mean_native_flow_ml_s": flow_means,
        "total_mean_inward_prescribed_flow_ml_s": -sum(flow_means.values()),
        "mean_vessel_R_poiseuille": fmean(
            float(values["R_poiseuille"]) for values in vessel_values
        ),
        "mean_vessel_L": fmean(float(values["L"]) for values in vessel_values),
        "mean_vessel_C": fmean(float(values["C"]) for values in vessel_values),
    }


def find_svzerodsolver(explicit: Path | None) -> Path:
    if explicit is not None:
        solver = explicit.expanduser().resolve()
    else:
        discovered = shutil.which("svzerodsolver")
        candidates = [
            Path(discovered) if discovered else None,
            Path.home() / "anaconda3" / "bin" / "svzerodsolver",
        ]
        solver = next((path for path in candidates if path and path.is_file()), None)
        if solver is None:
            raise FileNotFoundError(
                "svzerodsolver was not found; pass --svzerodsolver PATH"
            )
    if not solver.is_file():
        raise FileNotFoundError(f"svzerodsolver not found: {solver}")
    return solver


def run_baseline_solver(
    model_path: Path,
    model: dict[str, Any],
    solver_path: Path,
) -> dict[str, Any]:
    output_csv = model_path.with_name(f"{model_path.stem}_baseline.csv")
    completed = subprocess.run(
        [str(solver_path), str(model_path), str(output_csv)], check=False
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"svZeroDSolver failed with exit code {completed.returncode}"
        )

    with output_csv.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("svZeroDSolver produced an empty CSV")
    numeric_fields = ("time", "flow_in", "flow_out", "pressure_in", "pressure_out")
    if not all(
        math.isfinite(float(row[field])) for row in rows for field in numeric_fields
    ):
        raise ValueError("svZeroDSolver output contains non-finite values")

    by_vessel: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_vessel[row["name"]].append(row)
    attachment_to_vessel = {
        next(iter(vessel["boundary_conditions"].values())): vessel["vessel_name"]
        for vessel in model["vessels"]
        if vessel.get("boundary_conditions")
    }
    root_pairs = [
        (vessel["boundary_conditions"]["inlet"], vessel["vessel_name"])
        for vessel in model["vessels"]
        if "inlet" in vessel.get("boundary_conditions", {})
    ]
    if len(root_pairs) != 1:
        raise ValueError(f"Expected one rooted inlet vessel, found {len(root_pairs)}")
    root_name, root_vessel = root_pairs[0]
    mean = lambda values: sum(values) / len(values)
    root_inflow = mean(
        [float(row["flow_in"]) for row in by_vessel[root_vessel]]
    )
    terminal_sum = 0.0
    rcr_outflow = 0.0
    terminal_flow = 0.0
    for boundary_name, vessel_name in attachment_to_vessel.items():
        if boundary_name == root_name:
            continue
        value = mean(
            [float(row["flow_out"]) for row in by_vessel[vessel_name]]
        )
        terminal_sum += value
        if boundary_name.startswith("RCR_"):
            rcr_outflow += value
        elif boundary_name.startswith("FLOW_"):
            terminal_flow += value
    residual = root_inflow - terminal_sum
    scale = abs(root_inflow) + abs(terminal_sum)
    return {
        "status": "PASS",
        "solver": str(solver_path),
        "output_csv": str(output_csv),
        "number_of_rows": len(rows),
        "number_of_time_points_per_vessel": len(by_vessel[root_vessel]),
        "mean_root_inflow_ml_s": root_inflow,
        "mean_other_terminal_flow_ml_s": terminal_flow,
        "mean_rcr_outflow_ml_s": rcr_outflow,
        "mean_external_flow_residual_ml_s": residual,
        "relative_external_flow_residual": abs(residual) / scale if scale else 0.0,
    }


def write_report(
    patient_dir: Path,
    config: dict[str, Any],
    cfd: dict[str, Any],
    prepared: dict[str, Any],
    output_json: Path,
    validation: dict[str, Any],
) -> None:
    manifest = {
        "patient_id": config["patient_id"],
        "source_solver_xml": str(resolve(patient_dir, config["solver_xml"])),
        "source_surface_model": str(resolve(patient_dir, config["surface_model"])),
        "centerlines": str(prepared["centerlines"]),
        "raw_0d_json": str(prepared["raw_json"]),
        "multi_inlet_0d_json": str(output_json),
        "fluid": {
            "density": cfd["density"],
            "viscosity": cfd["viscosity"],
        },
        "simulation": {
            "time_step": cfd["time_step"],
            "number_of_time_steps": cfd["number_of_time_steps"],
        },
        "flow_boundaries": cfd["flow_boundaries"],
        "rcr_boundaries": cfd["rcr_boundaries"],
        "validation": validation,
    }
    manifest_path = prepared["output_dir"] / "generation_manifest.json"
    with manifest_path.open("w") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")

    report_path = prepared["output_dir"] / "generation_report.md"
    with report_path.open("w") as stream:
        stream.write(f"# {config['patient_id']} 3D-to-0D Generation Report\n\n")
        stream.write(f"- Status: **{validation['status']}**\n")
        stream.write(f"- Centerline root: `{prepared['root_name']}`\n")
        stream.write(f"- Geometry-derived vessels: {validation['number_of_vessels']}\n")
        stream.write(f"- Junctions: {validation['number_of_junctions']}\n")
        stream.write(f"- Prescribed-flow inlets: {validation['number_of_flow_inlets']}\n")
        stream.write(f"- RCR outlets: {validation['number_of_rcr_outlets']}\n")
        stream.write(
            "- Mean prescribed inflow: "
            f"{validation['total_mean_inward_prescribed_flow_ml_s']:.6g} mL/s\n"
        )
        stream.write(f"- Final svZeroD file: `{output_json.name}`\n\n")
        if "baseline_solver" in validation:
            baseline = validation["baseline_solver"]
            stream.write("## Baseline solver check\n\n")
            stream.write(f"- Status: **{baseline['status']}**\n")
            stream.write(
                "- Mean RCR outflow: "
                f"{baseline['mean_rcr_outflow_ml_s']:.6g} mL/s\n"
            )
            stream.write(
                "- Mean external-flow residual: "
                f"{baseline['mean_external_flow_residual_ml_s']:.6g} mL/s\n"
            )
            stream.write(
                "- Relative external-flow residual: "
                f"{baseline['relative_external_flow_residual']:.6g}\n\n"
            )
        stream.write("## Multi-inlet convention\n\n")
        stream.write(
            "The root cap is attached to the root vessel inlet and uses positive "
            "flow into the network. Other physical inlets terminate branches of "
            "the rooted centerline tree and therefore retain the negative native "
            "CFD surface-flux sign.\n"
        )


def main() -> int:
    args = parse_args()
    patient_dir = args.patient_dir.expanduser().resolve()
    config_path = (
        args.config.expanduser().resolve()
        if args.config
        else patient_dir / "patient_0d_config.json"
    )
    config = read_config(patient_dir, config_path)
    cfd = parse_solver_xml(resolve(patient_dir, config["solver_xml"]), patient_dir)
    prepared = prepare_inputs(patient_dir, config, cfd)
    print(
        f"Prepared {len(cfd['flow_boundaries'])} FLOW and "
        f"{len(cfd['rcr_boundaries'])} RCR boundaries in {prepared['output_dir']}"
    )
    if args.prepare_only:
        print("Preparation complete; SimVascular was not launched.")
        return 0

    run_simvascular(
        args.simvascular.expanduser().resolve(),
        prepared["runtime_config"],
        args.reuse_centerlines,
    )
    output_json, model = adapt_multi_inlet_model(patient_dir, config, cfd, prepared)
    validation = validate_model(model, cfd, prepared)
    if args.run_baseline:
        validation["baseline_solver"] = run_baseline_solver(
            output_json,
            model,
            find_svzerodsolver(args.svzerodsolver),
        )
    write_report(patient_dir, config, cfd, prepared, output_json, validation)
    print(json.dumps(validation, indent=2))
    print(f"Generated patient 0D model: {output_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
