#!/usr/bin/env python3
"""Generate a rigid geometric svZeroD model with VMTK, without SimVascular.

Run this script with the VMTK Conda Python. It extracts centerlines from the
surface mesh using cap-centroid point seeds, converts the resulting rooted tree
to one R/L/C vessel per VMTK branch, and restores exact FLOW/RCR/Resistance
conditions from solver.xml.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import defaultdict, deque
from pathlib import Path
from statistics import mean
from typing import Any, Optional, Sequence

try:
    import vtk
    from vmtk import vmtkscripts
except ImportError as error:
    raise SystemExit(
        "Run with the VMTK Python environment"
    ) from error


FLOAT_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def as_float(text: Optional[str], field: str) -> float:
    match = FLOAT_PATTERN.search((text or "").strip())
    if match is None:
        raise ValueError(f"Invalid numeric value for {field}: {text!r}")
    value = float(match.group(0))
    if not math.isfinite(value): raise ValueError(f"Non-finite numeric value for {field}")
    return value


def parse_bool(text: Optional[str]) -> bool:
    return (text or "").strip().lower() in {"1", "true", "yes", "on"}


def read_flow_file(path: Path) -> tuple[list[float], list[float]]:
    rows: list[tuple[float, float]] = []
    for line in path.read_text().splitlines():
        fields = re.split(r"[\s,]+", line.strip())
        if len(fields) < 2: continue
        try: time, flow = float(fields[0]), float(fields[1])
        except ValueError: continue
        if not rows and time > 2.0 and flow > 2.0: continue
        rows.append((time, flow))
    if len(rows) < 2 or any(next_time <= time for (time, _), (next_time, _) in zip(rows, rows[1:])):
        raise ValueError(f"Invalid flow waveform: {path}")
    return [time for time, _ in rows], [flow for _, flow in rows]


def parse_solver_xml(path: Path) -> dict[str, Any]:
    root = __import__("xml.etree.ElementTree", fromlist=["ElementTree"]).fromstring(path.read_text(errors="replace"))
    general = root.find("GeneralSimulationParameters")
    fluid = next((item for item in root.findall("Add_equation") if item.attrib.get("type", "").lower() == "fluid"), None)
    if general is None or fluid is None: raise ValueError("solver.xml requires GeneralSimulationParameters and fluid equation.")
    dt = as_float(general.findtext("Time_step_size"), "time step")
    steps = int(as_float(general.findtext("Number_of_time_steps"), "number of time steps"))
    parsed: dict[str, Any] = {"time_step": dt, "number_of_time_steps": steps, "density": as_float(fluid.findtext("Density"), "density"), "viscosity": as_float(fluid.findtext("Viscosity/Value"), "viscosity"), "flow_boundaries": {}, "rcr_boundaries": {}, "resistance_boundaries": {}}
    for bc in fluid.findall("Add_BC"):
        name = (bc.attrib.get("name") or "").strip(); bc_type = (bc.findtext("Type") or "").lower(); dependence = (bc.findtext("Time_dependence") or "").lower()
        if bc_type in {"dir", "dirichlet"} and parse_bool(bc.findtext("Impose_flux")):
            if dependence == "unsteady":
                waveform = path.parent / (bc.findtext("Temporal_values_file_path") or "")
                times, flows = read_flow_file(waveform)
            elif dependence == "steady":
                value = as_float(bc.findtext("Value"), f"{name} flow"); times, flows = [0.0, max(dt * steps, dt)], [value, value]
            else: raise ValueError(f"Unsupported FLOW dependence for {name}: {dependence}")
            parsed["flow_boundaries"][name] = {"t": times, "Q_native": flows}
        elif bc_type in {"neu", "neumann"} and dependence == "rcr":
            values = bc.find("RCR_values")
            if values is None: raise ValueError(f"RCR BC {name} lacks RCR_values")
            parsed["rcr_boundaries"][name] = {"Rp": as_float(values.findtext("Proximal_resistance"), f"{name} Rp"), "C": as_float(values.findtext("Capacitance"), f"{name} C"), "Rd": as_float(values.findtext("Distal_resistance"), f"{name} Rd"), "Pd": as_float(values.findtext("Distal_pressure"), f"{name} Pd")}
        elif bc_type in {"neu", "neumann"} and dependence == "resistance":
            parsed["resistance_boundaries"][name] = {"R": as_float(bc.findtext("Value"), f"{name} R"), "Pd": as_float(bc.findtext("Distal_pressure") or "0", f"{name} Pd")}
    if not parsed["flow_boundaries"] or not (parsed["rcr_boundaries"] or parsed["resistance_boundaries"]): raise ValueError("solver.xml needs at least one FLOW and terminal outlet BC.")
    return parsed


def validate_model(model: dict[str, Any], cfd: dict[str, Any], root: str) -> None:
    expected = {f"FLOW_{name}" for name in cfd["flow_boundaries"]} | {f"RCR_{name}" for name in cfd["rcr_boundaries"]} | {f"RESISTANCE_{name}" for name in cfd["resistance_boundaries"]}
    boundaries = {item["bc_name"] for item in model["boundary_conditions"]}
    attached = [name for vessel in model["vessels"] for name in vessel["boundary_conditions"].values()]
    if boundaries != expected or set(attached) != boundaries or len(attached) != len(boundaries): raise ValueError("Final model does not attach every solver.xml BC exactly once.")
    for vessel in model["vessels"]:
        for key in ("R_poiseuille", "L", "C", "stenosis_coefficient"):
            value = vessel["zero_d_element_values"].get(key)
            if value is None or not math.isfinite(float(value)) or float(value) < 0: raise ValueError(f"Invalid {key} on {vessel['vessel_name']}")


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--step0-config", type=Path, default=Path("step0_config.json"))
    parser.add_argument("--root-surface", help="Override metrics.fontan_pressure_face")
    parser.add_argument("--output-directory", default="0d")
    parser.add_argument("--vessel-capacitance", type=float, default=0.0, help="Rigid-wall geometric vessel capacitance")
    parser.add_argument("--overwrite-centerlines", action="store_true")
    return parser.parse_args()


def read_polydata(path: Path) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path)); reader.Update()
    data = reader.GetOutput()
    if data is None or data.GetNumberOfPoints() == 0:
        raise ValueError(f"Could not read nonempty VTP: {path}")
    return data


def centroid(polydata: vtk.vtkPolyData) -> list[float]:
    points = polydata.GetPoints()
    result = [0.0, 0.0, 0.0]
    for index in range(points.GetNumberOfPoints()):
        point = points.GetPoint(index)
        for coordinate in range(3): result[coordinate] += point[coordinate]
    return [value / points.GetNumberOfPoints() for value in result]


def write_polydata(path: Path, data: vtk.vtkPolyData) -> None:
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path)); writer.SetInputData(data)
    if not writer.Write(): raise RuntimeError(f"Could not write {path}")


def centerlines(surface: vtk.vtkPolyData, root: list[float], targets: list[list[float]]) -> vtk.vtkPolyData:
    calculator = vmtkscripts.vmtkCenterlines()
    calculator.Surface = surface
    calculator.SeedSelectorName = "pointlist"
    calculator.SourcePoints = root
    calculator.TargetPoints = [value for point in targets for value in point]
    calculator.AppendEndPoints = 1
    calculator.Execute()
    branches = vmtkscripts.vmtkBranchExtractor(); branches.Centerlines = calculator.Centerlines; branches.Execute()
    merger = vmtkscripts.vmtkCenterlineMerge(); merger.Centerlines = branches.Centerlines; merger.Execute()
    return merger.Centerlines


def line_records(data: vtk.vtkPolyData, density: float, viscosity: float, capacitance: float) -> list[dict]:
    radius = data.GetPointData().GetArray("MaximumInscribedSphereRadius")
    if radius is None: raise ValueError("VMTK centerlines lack MaximumInscribedSphereRadius.")
    blanking = data.GetCellData().GetArray("Blanking")
    ids = vtk.vtkIdList(); records = []
    for cell in range(data.GetNumberOfLines()):
        if blanking is not None and blanking.GetTuple1(cell) != 0: continue
        data.GetCellPoints(cell, ids)
        point_ids = [ids.GetId(i) for i in range(ids.GetNumberOfIds())]
        if len(point_ids) < 2: continue
        length = sum(math.dist(data.GetPoint(a), data.GetPoint(b)) for a, b in zip(point_ids, point_ids[1:]))
        radii = [radius.GetTuple1(index) for index in point_ids]
        effective_radius = sum(radii) / len(radii)
        if length <= 0 or effective_radius <= 0: raise ValueError("Nonpositive VMTK branch geometry.")
        area = math.pi * effective_radius ** 2
        records.append({
            "cell": cell, "start": point_ids[0], "end": point_ids[-1], "length": length,
            "radius": effective_radius, "R": 8.0 * viscosity * length / (math.pi * effective_radius ** 4),
            "L": density * length / area, "C": capacitance,
        })
    if not records: raise ValueError("VMTK produced no usable centerline branches.")
    return records


def endpoint_graph(records: list[dict]) -> tuple[dict[int, list[int]], dict[int, list[int]]]:
    edges: dict[int, list[int]] = defaultdict(list)
    nodes: dict[int, list[int]] = defaultdict(list)
    for index, record in enumerate(records):
        edges[record["start"]].append(index); edges[record["end"]].append(index)
        nodes[record["start"]].append(index); nodes[record["end"]].append(index)
    return edges, nodes


def nearest_leaf(data: vtk.vtkPolyData, leaf_nodes: list[int], point: list[float], label: str) -> int:
    node = min(leaf_nodes, key=lambda item: math.dist(data.GetPoint(item), point))
    if math.dist(data.GetPoint(node), point) > 2.0:
        raise ValueError(f"Cap '{label}' did not map near a VMTK centerline leaf; inspect cap topology.")
    return node


def oriented_tree(records: list[dict], root_node: int) -> tuple[dict[int, tuple[int, int]], dict[int, int]]:
    adjacency, _ = endpoint_graph(records)
    orientation: dict[int, tuple[int, int]] = {}; parents: dict[int, int] = {root_node: -1}; queue = deque([root_node])
    while queue:
        node = queue.popleft()
        for edge in adjacency[node]:
            record = records[edge]; other = record["end"] if record["start"] == node else record["start"]
            if other in parents: continue
            parents[other] = edge; orientation[edge] = (node, other); queue.append(other)
    if len(orientation) != len(records):
        raise ValueError("VMTK centerline network is not a connected rooted tree.")
    return orientation, parents


def boundary(cfd: dict, cap: str, is_root: bool) -> dict:
    if cap in cfd["flow_boundaries"]:
        source = cfd["flow_boundaries"][cap]
        values = source["Q_native"]
        return {"bc_name": f"FLOW_{cap}", "bc_type": "FLOW", "bc_values": {"t": source["t"], "Q": [-q for q in values] if is_root else values}}
    if cap in cfd["rcr_boundaries"]:
        return {"bc_name": f"RCR_{cap}", "bc_type": "RCR", "bc_values": dict(cfd["rcr_boundaries"][cap])}
    if cap in cfd["resistance_boundaries"]:
        return {"bc_name": f"RESISTANCE_{cap}", "bc_type": "RESISTANCE", "bc_values": dict(cfd["resistance_boundaries"][cap])}
    raise ValueError(f"No supported solver.xml boundary for cap '{cap}'.")


def main() -> None:
    options = args(); case = options.case_directory.resolve()
    config_path = options.step0_config if options.step0_config.is_absolute() else case / options.step0_config
    config = json.loads(config_path.read_text()); mesh = config["mesh"]
    root_cap = options.root_surface or config["metrics"]["fontan_pressure_face"]
    if root_cap not in {item["name"] for item in config["inlets"]}: raise ValueError("Root must be a Step 0 inlet cap.")
    output = Path(options.output_directory); output = output if output.is_absolute() else case / output; output.mkdir(parents=True, exist_ok=True)
    mesh_dir = case / mesh.get("directory", "mesh"); cap_dir = mesh_dir / mesh.get("surface_directory", "mesh-surfaces")
    cfd = parse_solver_xml(case / "solver.xml")
    cap_names = sorted(set(cfd["flow_boundaries"]) | set(cfd["rcr_boundaries"]) | set(cfd["resistance_boundaries"]))
    cap_data = {name: read_polydata(cap_dir / f"{name}.vtp") for name in cap_names}
    surface = read_polydata(mesh_dir / "mesh-complete.exterior.vtp")
    centerline_file = output / "centerlines.vtp"
    if centerline_file.is_file() and not options.overwrite_centerlines:
        branches = read_polydata(centerline_file)
    else:
        branches = centerlines(surface, centroid(cap_data[root_cap]), [centroid(cap_data[name]) for name in cap_names if name != root_cap])
        write_polydata(centerline_file, branches)
    records = line_records(branches, cfd["density"], cfd["viscosity"], options.vessel_capacitance)
    graph, nodes = endpoint_graph(records); leaves = [node for node, edge_ids in nodes.items() if len(edge_ids) == 1]
    cap_nodes = {name: nearest_leaf(branches, leaves, centroid(data), name) for name, data in cap_data.items()}
    if len(set(cap_nodes.values())) != len(cap_nodes): raise ValueError("Multiple caps mapped to one centerline leaf.")
    orientation, _ = oriented_tree(records, cap_nodes[root_cap])
    vessels = []; vessel_for_edge = {}
    for index, record in enumerate(records):
        vessel_for_edge[index] = index
        vessels.append({"vessel_id": index, "vessel_name": f"branch{index}_seg0", "vessel_length": record["length"], "zero_d_element_type": "BloodVessel", "zero_d_element_values": {"R_poiseuille": record["R"], "L": record["L"], "C": record["C"], "stenosis_coefficient": 0.0}, "boundary_conditions": {}})
    for cap, node in cap_nodes.items():
        edge = nodes[node][0]; parent, child = orientation[edge]
        vessels[vessel_for_edge[edge]]["boundary_conditions"]["inlet" if node == parent else "outlet"] = boundary(cfd, cap, cap == root_cap)["bc_name"]
    junctions = []
    for node, edges in nodes.items():
        if len(edges) <= 2: continue
        incoming = [vessel_for_edge[edge] for edge in edges if orientation[edge][1] == node]
        outgoing = [vessel_for_edge[edge] for edge in edges if orientation[edge][0] == node]
        if len(incoming) != 1 or not outgoing: raise ValueError("Invalid rooted VMTK junction orientation.")
        junctions.append({"junction_name": f"J{len(junctions)}", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": incoming, "outlet_vessels": outgoing})
    model = {"boundary_conditions": [boundary(cfd, cap, cap == root_cap) for cap in cap_names], "junctions": junctions, "simulation_parameters": {"model_name": f"{config['case_name']}_vmtk", "density": cfd["density"], "viscosity": cfd["viscosity"], "number_of_cardiac_cycles": int(config["simulation"].get("number_of_cycles", 8)), "number_of_time_pts_per_cardiac_cycle": int(config["simulation"].get("steps_per_cycle", 1000))}, "vessels": vessels}
    validate_model(model, cfd, root_cap)
    model_file = output / f"{config['case_name']}_0d.json"; model_file.write_text(json.dumps(model, indent=2) + "\n")
    print(f"Wrote {centerline_file}")
    print(f"Wrote {model_file} ({len(vessels)} vessels, {len(junctions)} junctions)")


if __name__ == "__main__":
    main()
