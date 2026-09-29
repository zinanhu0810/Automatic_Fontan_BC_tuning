#!/usr/bin/env python3
"""Generate an initialized svMultiPhysics solver.xml from JSON input.

The script intentionally uses only the Python standard library.  It sizes RCR
outlet resistances from cap area and the clinical mean Fontan pressure/split,
rescales pulsatile inflow files to the requested cardiac period, and computes
phase-aware periodic RCR capacitor initial pressures before writing solver.xml.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from statistics import fmean
from xml.dom import minidom


MMHG_TO_DYN_CM2 = 1333.22


def read_config(path: Path) -> dict:
    with path.open() as stream:
        config = json.load(stream)
    for key in ("case_name", "mesh", "inlets", "outlets", "targets", "simulation"):
        if key not in config:
            raise ValueError(f"Configuration is missing required key '{key}'.")
    if not config["inlets"]:
        raise ValueError("At least one inlet must be supplied.")
    return config


def parse_cap_info(path: Path) -> dict[str, float]:
    caps: dict[str, float] = {}
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields or line.lstrip().startswith("#"):
            continue
        if len(fields) != 2:
            raise ValueError(f"{path}:{number}: expected 'cap_name area_cm2'.")
        caps[fields[0]] = float(fields[1])
    if not caps:
        raise ValueError(f"No cap areas were read from {path}.")
    return caps


def read_waveform(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        fields = line.split()
        if len(fields) < 2:
            continue
        try:
            time, flow = float(fields[0]), float(fields[1])
        except ValueError:
            continue
        # SimVascular flow files commonly begin with a numeric metadata line
        # such as "31 16" (number of samples and format/version). It is not a
        # time/flow sample and otherwise breaks the monotonic time check.
        if index == 0 and len(fields) == 2 and time > 1.0 and flow.is_integer():
            continue
        rows.append((time, flow))
    if len(rows) < 2:
        raise ValueError(f"{path} must contain at least two numeric time/flow rows.")
    if any(right[0] <= left[0] for left, right in zip(rows, rows[1:])):
        raise ValueError(f"Waveform times must be strictly increasing: {path}")
    return rows


def rescale_waveform(rows: list[tuple[float, float]], period: float) -> list[tuple[float, float]]:
    start, end = rows[0][0], rows[-1][0]
    if end <= start:
        raise ValueError("Waveform duration must be positive.")
    return [((time - start) * period / (end - start), flow) for time, flow in rows]


def mean_waveform(rows: list[tuple[float, float]]) -> float:
    integral = sum(0.5 * (t1 - t0) * (q0 + q1) for (t0, q0), (t1, q1) in zip(rows, rows[1:]))
    return integral / (rows[-1][0] - rows[0][0])


def interpolate(rows: list[tuple[float, float]], time: float) -> float:
    if time <= rows[0][0]:
        return rows[0][1]
    if time >= rows[-1][0]:
        return rows[-1][1]
    for (t0, q0), (t1, q1) in zip(rows, rows[1:]):
        if t0 <= time <= t1:
            return q0 + (q1 - q0) * (time - t0) / (t1 - t0)
    raise RuntimeError("Interpolation failed.")


def periodic_rcr_pressure(times: list[float], flows: list[float], resistance: float, capacitance: float, distal_pressure: float) -> float:
    """Periodic capacitor pressure using exact updates with midpoint flow."""
    if capacitance <= 0:
        return distal_pressure
    rd = resistance
    a_total, b_total = 1.0, 0.0
    for t0, t1, q0, q1 in zip(times, times[1:], flows, flows[1:]):
        decay = math.exp(-(t1 - t0) / (rd * capacitance))
        source = (1.0 - decay) * rd * 0.5 * (q0 + q1)
        b_total = decay * b_total + source
        a_total *= decay
    return distal_pressure + b_total / max(1.0 - a_total, 1e-15)


def indent_xml(root: ET.Element) -> str:
    raw = ET.tostring(root, encoding="utf-8")
    return "\n".join(line for line in minidom.parseString(raw).toprettyxml(indent="  ").splitlines() if line.strip()) + "\n"


def add_text(parent: ET.Element, tag: str, value: object) -> ET.Element:
    child = ET.SubElement(parent, tag)
    child.text = str(value)
    return child


def wall_paths(mesh_dir: Path, min_polygons: int) -> list[Path]:
    result = []
    for path in sorted(mesh_dir.glob("walls*.vtp")):
        header = path.read_bytes()[:4096].decode("ascii", errors="ignore")
        match = re.search(r'NumberOfPolys="(\d+)"', header)
        if match and int(match.group(1)) >= min_polygons:
            result.append(path)
    if not result:
        raise ValueError("No usable wall VTP was found. Set mesh.wall_surfaces explicitly if needed.")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--config", type=Path, default=Path("step0_config.json"))
    parser.add_argument("--output", type=Path, default=Path("."))
    args = parser.parse_args()

    case_dir = args.case_directory.resolve()
    config_path = args.config if args.config.is_absolute() else case_dir / args.config
    config = read_config(config_path)
    output = args.output if args.output.is_absolute() else case_dir / args.output
    output.mkdir(parents=True, exist_ok=True)
    mesh = config["mesh"]
    mesh_dir = case_dir / mesh.get("directory", "mesh")
    surface_dir = mesh_dir / mesh.get("surface_directory", "mesh-surfaces")
    volume = mesh_dir / mesh.get("volume_file", "mesh-complete.mesh.vtu")
    cap_areas = parse_cap_info(case_dir / mesh.get("cap_info_file", "cap_info.txt"))
    surface_paths = {path.stem: path for path in surface_dir.glob("*.vtp")}
    if not volume.is_file() or not surface_dir.is_dir():
        raise FileNotFoundError("Expected mesh volume file and mesh-surfaces directory were not found.")

    period = float(config["simulation"]["cardiac_period_s"])
    steps_per_cycle = int(config["simulation"].get("steps_per_cycle", 1000))
    cycles = int(config["simulation"].get("number_of_cycles", 8))
    time_step = float(config["simulation"].get("time_step_s", period / steps_per_cycle))
    if not math.isclose(time_step * steps_per_cycle, period, rel_tol=2e-3, abs_tol=1e-8):
        raise ValueError("time_step_s * steps_per_cycle must equal cardiac_period_s.")

    inflow_rows: dict[str, list[tuple[float, float]]] = {}
    inlet_means: dict[str, float] = {}
    for inlet in config["inlets"]:
        name = inlet["name"]
        if name not in surface_paths:
            raise ValueError(f"Inlet '{name}' has no matching cap VTP.")
        if "flow_file" in inlet:
            rows = rescale_waveform(read_waveform(case_dir / inlet["flow_file"]), period)
            out_name = f"{name}_physical.flow"
            with (output / out_name).open("w") as stream:
                # svMultiPhysics BCT files begin with numeric metadata rather
                # than a text header: number of samples, then Fourier modes.
                stream.write(f"{len(rows)} {max(1, len(rows) // 2)}\n")
                for time, flow in rows:
                    stream.write(f"{time:.10g} {flow:.12g}\n")
            inflow_rows[name] = rows
            inlet_means[name] = mean_waveform(rows)
            inlet["generated_flow_file"] = out_name
        else:
            inlet_means[name] = float(inlet["steady_flow_ml_s"])

    lpa = list(config["outlets"]["lpa"])
    rpa = list(config["outlets"]["rpa"])
    outlet_names = lpa + rpa
    declared = {item["name"] for item in config["inlets"]} | set(outlet_names)
    missing = declared - set(surface_paths)
    if missing:
        raise ValueError(f"Configured caps missing from mesh-surfaces: {sorted(missing)}")
    cap_paths = {name: surface_paths[name] for name in declared}
    if set(lpa) & set(rpa) or not lpa or not rpa:
        raise ValueError("LPA and RPA outlet lists must be non-empty and disjoint.")

    total_inflow = abs(sum(inlet_means.values()))
    if total_inflow <= 0:
        raise ValueError("Total mean inflow is zero.")
    target = config["targets"]
    lpa_fraction = float(target["lpa_fraction_percent"]) / 100.0
    if not 0.0 < lpa_fraction < 1.0:
        raise ValueError("targets.lpa_fraction_percent must be between 0 and 100.")
    pressure = float(target["mean_fontan_pressure_mmhg"]) * MMHG_TO_DYN_CM2
    group_flows = {"lpa": total_inflow * lpa_fraction, "rpa": total_inflow * (1.0 - lpa_fraction)}
    group_caps = {"lpa": lpa, "rpa": rpa}
    outlet_type = str(config.get("outlet_boundary_type", "rcr")).lower()
    if outlet_type not in {"rcr", "resistance"}:
        raise ValueError("outlet_boundary_type must be 'rcr' or 'resistance'.")
    rcr_config = config.get("rcr", {})
    capacitance = float(rcr_config.get("capacitance_s", 1.5e-4))
    rp_fraction = float(rcr_config.get("proximal_resistance_fraction", 1.0 / 11.0))
    if outlet_type == "rcr" and not 0 < rp_fraction < 1:
        raise ValueError("rcr.proximal_resistance_fraction must be between 0 and 1.")

    rcr: dict[str, dict[str, float]] = {}
    for side, names in group_caps.items():
        area_total = sum(cap_areas[name] for name in names)
        if area_total <= 0:
            raise ValueError(f"{side.upper()} cap area is non-positive.")
        equivalent = pressure / group_flows[side]
        for name in names:
            total_r = equivalent / (cap_areas[name] / area_total)
            if outlet_type == "rcr":
                rcr[name] = {"rp": total_r * rp_fraction, "rd": total_r * (1 - rp_fraction), "c": capacitance}
            else:
                rcr[name] = {"resistance": total_r}

    if outlet_type == "rcr" and inflow_rows:
        time_grid = sorted({time for rows in inflow_rows.values() for time, _ in rows})
        stored_total = [sum(interpolate(rows, time) for rows in inflow_rows.values()) for time in time_grid]
        physical_total = [-value for value in stored_total]
        conductance = {name: 1.0 / (values["rp"] + values["rd"]) for name, values in rcr.items()}
        denom = sum(conductance.values())
        for name, values in rcr.items():
            outlet_flow = [value * conductance[name] / denom for value in physical_total]
            values["initial_pressure"] = periodic_rcr_pressure(time_grid, outlet_flow, values["rd"], values["c"], 0.0)
    elif outlet_type == "rcr":
        initial = float(config["rcr"].get("initial_pressure_dyn_cm2", pressure))
        for values in rcr.values():
            values["initial_pressure"] = initial

    root = ET.Element("svMultiPhysicsFile", {"version": "0.1"})
    general = ET.SubElement(root, "GeneralSimulationParameters")
    for tag, value in (("Continue_previous_simulation", "false"), ("Number_of_spatial_dimensions", 3), ("Number_of_time_steps", cycles * steps_per_cycle), ("Time_step_size", f"{time_step:.10g}"), ("Spectral_radius_of_infinite_time_step", "0.50"), ("Searched_file_name_to_trigger_stop", "STOP_SIM"), ("Save_results_to_VTK_format", 1), ("Name_prefix_of_saved_VTK_files", "result"), ("Increment_in_saving_VTK_files", config["simulation"].get("vtk_interval", 50)), ("Start_saving_after_time_step", 1), ("Increment_in_saving_restart_files", steps_per_cycle), ("Convert_BIN_to_VTK_format", 0), ("Verbose", 1), ("Warning", 0), ("Debug", 0)):
        add_text(general, tag, value)
    add_mesh = ET.SubElement(root, "Add_mesh", {"name": "msh"})
    add_text(add_mesh, "Mesh_file_path", os.path.relpath(volume, output).replace(os.sep, "/"))
    if "scale_factor" in mesh:
        add_text(add_mesh, "Mesh_scale_factor", mesh["scale_factor"])
    for name, path in sorted(cap_paths.items()):
        face = ET.SubElement(add_mesh, "Add_face", {"name": name})
        add_text(face, "Face_file_path", os.path.relpath(path, output).replace(os.sep, "/"))
    walls = [mesh_dir / value for value in mesh["wall_surfaces"]] if "wall_surfaces" in mesh else wall_paths(mesh_dir, int(mesh.get("minimum_wall_polygons", 10)))
    for path in walls:
        face = ET.SubElement(add_mesh, "Add_face", {"name": path.stem})
        add_text(face, "Face_file_path", os.path.relpath(path, output).replace(os.sep, "/"))

    fluid = ET.SubElement(root, "Add_equation", {"type": "fluid"})
    for tag, value in (("Coupled", "true"), ("Min_iterations", 2), ("Max_iterations", 7), ("Tolerance", "1e-4"), ("Backflow_stabilization_coefficient", "0.2"), ("Density", "1.06")):
        add_text(fluid, tag, value)
    viscosity = ET.SubElement(fluid, "Viscosity", {"model": "Constant"}); add_text(viscosity, "Value", "0.04")
    spatial = ET.SubElement(fluid, "Output", {"type": "Spatial"})
    for tag in ("Velocity", "Pressure", "Traction", "Vorticity", "Divergence", "WSS"): add_text(spatial, tag, "true")
    bint = ET.SubElement(fluid, "Output", {"type": "B_INT"}); add_text(bint, "Pressure", "true"); add_text(bint, "Velocity", "true")
    vint = ET.SubElement(fluid, "Output", {"type": "V_INT"}); add_text(vint, "Pressure", "true")
    ls = ET.SubElement(fluid, "LS", {"type": "NS"}); la = ET.SubElement(ls, "Linear_algebra", {"type": "fsils"}); add_text(la, "Preconditioner", "fsils")
    for tag, value in (("Max_iterations", 10), ("NS_GM_max_iterations", 200), ("NS_CG_max_iterations", 500), ("Tolerance", "0.4"), ("NS_GM_tolerance", "0.01"), ("NS_CG_tolerance", "0.2"), ("Absolute_tolerance", "1e-17"), ("Krylov_space_dimension", 300)): add_text(ls, tag, value)
    for name, values in rcr.items():
        bc = ET.SubElement(fluid, "Add_BC", {"name": name}); add_text(bc, "Type", "Neu")
        if outlet_type == "rcr":
            add_text(bc, "Time_dependence", "RCR")
            values_xml = ET.SubElement(bc, "RCR_values")
            for tag, value in (("Capacitance", values["c"]), ("Distal_resistance", values["rd"]), ("Proximal_resistance", values["rp"]), ("Distal_pressure", 0), ("Initial_pressure", values["initial_pressure"])): add_text(values_xml, tag, f"{value:.12g}" if isinstance(value, float) else value)
        else:
            add_text(bc, "Time_dependence", "Resistance")
            add_text(bc, "Value", f"{values['resistance']:.12g}")
    for inlet in config["inlets"]:
        bc = ET.SubElement(fluid, "Add_BC", {"name": inlet["name"]}); add_text(bc, "Type", "Dir")
        if "flow_file" in inlet:
            add_text(bc, "Time_dependence", "Unsteady"); add_text(bc, "Temporal_values_file_path", inlet["generated_flow_file"])
        else:
            add_text(bc, "Time_dependence", "Steady"); add_text(bc, "Value", inlet["steady_flow_ml_s"])
        add_text(bc, "Profile", inlet.get("profile", "Parabolic")); add_text(bc, "Impose_flux", "true")
    for wall in walls:
        bc = ET.SubElement(fluid, "Add_BC", {"name": wall.stem}); add_text(bc, "Type", "Dir"); add_text(bc, "Time_dependence", "Steady"); add_text(bc, "Value", 0.0)

    (output / "solver.xml").write_text(indent_xml(root))
    print(f"Wrote {output / 'solver.xml'}")


if __name__ == "__main__":
    main()
