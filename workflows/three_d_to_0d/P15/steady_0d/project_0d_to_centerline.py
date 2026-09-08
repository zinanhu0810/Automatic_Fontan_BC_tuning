#!/usr/bin/env python3
"""Project steady svZeroD pressure and flow results onto a ROM centerline."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import vtk
from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
DEFAULT_MODEL = HERE / "P15_0d_resistance_constant_flow.json"
DEFAULT_RESULT = HERE / "P15_0d_resistance_constant_flow.csv"
DEFAULT_CENTERLINE = HERE.parent / "rom" / "centerlines.vtp"
DEFAULT_OUTPUT = HERE / "P15_0d_pressure_centerlines.vtp"
DEFAULT_NODES = HERE / "P15_0d_pressure_centerline_nodes.csv"
DEFAULT_FIGURE = HERE / "P15_0d_pressure_centerlines.png"
MMHG_PER_DYN_CM2 = 1.0 / 1333.22


def read_polydata(path: Path) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()
    if data.GetNumberOfPoints() == 0:
        raise ValueError(f"Centerline has no points: {path}")
    return data


def read_vessel_results(path: Path) -> dict[str, dict[str, float]]:
    totals: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "count": 0.0,
            "flow_in": 0.0,
            "flow_out": 0.0,
            "pressure_in": 0.0,
            "pressure_out": 0.0,
        }
    )
    with path.open() as stream:
        for row in csv.DictReader(stream):
            values = totals[row["name"]]
            values["count"] += 1.0
            for key in ("flow_in", "flow_out", "pressure_in", "pressure_out"):
                values[key] += float(row[key])
    if not totals:
        raise ValueError(f"Result CSV is empty: {path}")
    return {
        name: {
            key: value / values["count"]
            for key, value in values.items()
            if key != "count"
        }
        for name, values in totals.items()
    }


def model_segments(model: dict) -> dict[int, list[dict]]:
    branches: dict[int, list[dict]] = defaultdict(list)
    pattern = re.compile(r"branch(-?\d+)_seg(\d+)$")
    for vessel in model["vessels"]:
        match = pattern.fullmatch(vessel["vessel_name"])
        if not match:
            raise ValueError(f"Unexpected vessel name: {vessel['vessel_name']}")
        branches[int(match.group(1))].append(
            {
                "segment_id": int(match.group(2)),
                "vessel_id": int(vessel["vessel_id"]),
                "vessel_name": vessel["vessel_name"],
                "length": float(vessel["vessel_length"]),
            }
        )
    for segments in branches.values():
        segments.sort(key=lambda item: item["segment_id"])
        cumulative = 0.0
        for segment in segments:
            segment["path_start"] = cumulative
            cumulative += segment["length"]
            segment["path_end"] = cumulative
    return dict(branches)


def find_segment(segments: list[dict], path_value: float) -> dict:
    tolerance = 1e-8 * max(1.0, segments[-1]["path_end"])
    for segment in segments:
        if path_value <= segment["path_end"] + tolerance:
            return segment
    raise ValueError(
        f"Centerline path {path_value} exceeds model branch length "
        f"{segments[-1]['path_end']}"
    )


def graph_adjacency(data: vtk.vtkPolyData) -> list[set[int]]:
    adjacency = [set() for _ in range(data.GetNumberOfPoints())]
    for cell_id in range(data.GetNumberOfCells()):
        ids = data.GetCell(cell_id).GetPointIds()
        for index in range(ids.GetNumberOfIds() - 1):
            first = ids.GetId(index)
            second = ids.GetId(index + 1)
            adjacency[first].add(second)
            adjacency[second].add(first)
    return adjacency


def fill_junction_pressures(
    values: list[float], known: list[bool], adjacency: list[set[int]]
) -> None:
    unknown = [index for index, is_known in enumerate(known) if not is_known]
    if not unknown:
        return
    known_values = [value for value, is_known in zip(values, known) if is_known]
    initial = sum(known_values) / len(known_values)
    for index in unknown:
        values[index] = initial

    for _ in range(50000):
        maximum_change = 0.0
        for index in unknown:
            neighbors = adjacency[index]
            if not neighbors:
                continue
            updated = sum(values[neighbor] for neighbor in neighbors) / len(neighbors)
            maximum_change = max(maximum_change, abs(updated - values[index]))
            values[index] = updated
        if maximum_change < 1e-10:
            return
    raise RuntimeError("Junction pressure interpolation did not converge")


def add_array(
    data: vtk.vtkPolyData, name: str, values: list[float], integer: bool = False
) -> None:
    array = vtk.vtkIntArray() if integer else vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples(len(values))
    for index, value in enumerate(values):
        if integer:
            array.SetValue(index, int(value))
        else:
            array.SetValue(index, float(value))
    data.GetPointData().AddArray(array)


def project(model: dict, results: dict[str, dict[str, float]], data: vtk.vtkPolyData):
    branch_array = data.GetPointData().GetArray("BranchId")
    path_array = data.GetPointData().GetArray("Path")
    if branch_array is None or path_array is None:
        raise ValueError("Centerline requires BranchId and Path point arrays")

    branches = model_segments(model)
    number_of_points = data.GetNumberOfPoints()
    pressure_dyn = [math.nan] * number_of_points
    flow = [math.nan] * number_of_points
    vessel_ids = [-1] * number_of_points
    segment_ids = [-1] * number_of_points
    known_pressure = [False] * number_of_points

    for point_id in range(number_of_points):
        branch_id = int(branch_array.GetTuple1(point_id))
        if branch_id < 0:
            continue
        if branch_id not in branches:
            raise ValueError(f"No 0D vessels found for centerline branch {branch_id}")
        path_value = float(path_array.GetTuple1(point_id))
        segment = find_segment(branches[branch_id], path_value)
        result = results[segment["vessel_name"]]
        local = (path_value - segment["path_start"]) / segment["length"]
        local = min(1.0, max(0.0, local))
        pressure_dyn[point_id] = (
            result["pressure_in"]
            + local * (result["pressure_out"] - result["pressure_in"])
        )
        flow[point_id] = result["flow_in"] + local * (
            result["flow_out"] - result["flow_in"]
        )
        vessel_ids[point_id] = segment["vessel_id"]
        segment_ids[point_id] = segment["segment_id"]
        known_pressure[point_id] = True

    fill_junction_pressures(pressure_dyn, known_pressure, graph_adjacency(data))
    pressure_mmhg = [value * MMHG_PER_DYN_CM2 for value in pressure_dyn]
    add_array(data, "Pressure_dyn_cm2", pressure_dyn)
    add_array(data, "Pressure_mmHg", pressure_mmhg)
    add_array(data, "Flow_ml_s", flow)
    add_array(data, "ROM_VesselId", vessel_ids, integer=True)
    add_array(data, "ROM_SegmentId", segment_ids, integer=True)
    data.GetPointData().SetActiveScalars("Pressure_mmHg")
    return pressure_mmhg, flow, vessel_ids, segment_ids


def write_polydata(path: Path, data: vtk.vtkPolyData) -> None:
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(data)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write {path}")


def write_nodes(
    path: Path,
    data: vtk.vtkPolyData,
    pressure: list[float],
    flow: list[float],
    vessel_ids: list[int],
    segment_ids: list[int],
) -> None:
    branch_array = data.GetPointData().GetArray("BranchId")
    path_array = data.GetPointData().GetArray("Path")
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "point_id",
                "x",
                "y",
                "z",
                "branch_id",
                "segment_id",
                "vessel_id",
                "path_cm",
                "pressure_mmHg",
                "flow_ml_s",
            ]
        )
        for point_id in range(data.GetNumberOfPoints()):
            x, y, z = data.GetPoint(point_id)
            writer.writerow(
                [
                    point_id,
                    x,
                    y,
                    z,
                    int(branch_array.GetTuple1(point_id)),
                    segment_ids[point_id],
                    vessel_ids[point_id],
                    path_array.GetTuple1(point_id),
                    pressure[point_id],
                    flow[point_id],
                ]
            )


COLOR_STOPS = (
    (0.00, (26, 51, 140)),
    (0.38, (10, 148, 163)),
    (0.72, (242, 199, 46)),
    (1.00, (184, 15, 31)),
)


def pressure_color(value: float, low: float, high: float) -> tuple[int, int, int]:
    fraction = 0.5 if high == low else (value - low) / (high - low)
    fraction = min(1.0, max(0.0, fraction))
    for (first_x, first_color), (second_x, second_color) in zip(
        COLOR_STOPS, COLOR_STOPS[1:]
    ):
        if fraction <= second_x:
            local = (fraction - first_x) / (second_x - first_x)
            return tuple(
                round(first + local * (second - first))
                for first, second in zip(first_color, second_color)
            )
    return COLOR_STOPS[-1][1]


def principal_axes(data: vtk.vtkPolyData):
    points = np.array([data.GetPoint(i) for i in range(data.GetNumberOfPoints())])
    center = points.mean(axis=0)
    _, _, vectors = np.linalg.svd(points - center, full_matrices=False)
    diagonal = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0)))
    return center, vectors, diagonal


def render(path: Path, data: vtk.vtkPolyData, pressure: list[float]) -> None:
    finite = [value for value in pressure if math.isfinite(value)]
    low, high = min(finite), max(finite)
    center, axes, _ = principal_axes(data)
    points = np.array([data.GetPoint(i) for i in range(data.GetNumberOfPoints())])
    centered = points - center
    projections = (
        np.column_stack((centered @ axes[0], centered @ axes[1])),
        np.column_stack((centered @ axes[0], centered @ axes[2])),
    )

    width, height = 1900, 980
    image = Image.new("RGB", (width, height), (251, 251, 251))
    draw = ImageDraw.Draw(image)
    font_path = "/System/Library/Fonts/Supplemental/Arial.ttf"
    bold_path = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
    regular = ImageFont.truetype(font_path, 24)
    small = ImageFont.truetype(font_path, 19)
    bold = ImageFont.truetype(bold_path, 29)
    panel_font = ImageFont.truetype(bold_path, 23)

    draw.text(
        (48, 26),
        f"P15 steady 0D pressure projected onto centerline  "
        f"({low:.3f}-{high:.3f} mmHg)",
        fill=(22, 22, 22),
        font=bold,
    )
    panel_boxes = ((45, 98, 910, 915), (955, 98, 1728, 915))

    def screen_coordinates(projected, box):
        left, top, right, bottom = box
        padding = 46
        minimum = projected.min(axis=0)
        maximum = projected.max(axis=0)
        spans = np.maximum(maximum - minimum, 1e-12)
        scale = min(
            (right - left - 2 * padding) / spans[0],
            (bottom - top - 2 * padding) / spans[1],
        )
        used = spans * scale
        x_offset = left + (right - left - used[0]) / 2.0
        y_offset = top + (bottom - top - used[1]) / 2.0
        screen = np.empty_like(projected)
        screen[:, 0] = x_offset + (projected[:, 0] - minimum[0]) * scale
        screen[:, 1] = bottom - (
            y_offset - top + (projected[:, 1] - minimum[1]) * scale
        )
        return screen

    for panel_index, (projected, box) in enumerate(zip(projections, panel_boxes)):
        screen = screen_coordinates(projected, box)
        draw.rectangle(box, outline=(205, 205, 205), width=2)
        draw.text(
            (box[0] + 18, box[1] + 15),
            "Principal view" if panel_index == 0 else "Orthogonal view",
            fill=(35, 35, 35),
            font=panel_font,
        )

        edges = []
        for cell_id in range(data.GetNumberOfCells()):
            ids = data.GetCell(cell_id).GetPointIds()
            for index in range(ids.GetNumberOfIds() - 1):
                first, second = ids.GetId(index), ids.GetId(index + 1)
                edges.append((first, second))
        for first, second in edges:
            coordinates = [tuple(screen[first]), tuple(screen[second])]
            draw.line(coordinates, fill=(55, 55, 55), width=10)
        for first, second in edges:
            value = (pressure[first] + pressure[second]) / 2.0
            coordinates = [tuple(screen[first]), tuple(screen[second])]
            draw.line(
                coordinates,
                fill=pressure_color(value, low, high),
                width=7,
            )

    bar_left, bar_top, bar_right, bar_bottom = 1790, 185, 1828, 815
    for y in range(bar_top, bar_bottom + 1):
        fraction = 1.0 - (y - bar_top) / (bar_bottom - bar_top)
        value = low + fraction * (high - low)
        draw.line(
            ((bar_left, y), (bar_right, y)),
            fill=pressure_color(value, low, high),
            width=1,
        )
    draw.rectangle(
        (bar_left, bar_top, bar_right, bar_bottom), outline=(70, 70, 70), width=2
    )
    draw.text((1748, 130), "Pressure", fill=(30, 30, 30), font=regular)
    draw.text((1762, 156), "[mmHg]", fill=(30, 30, 30), font=small)
    for index in range(5):
        fraction = index / 4.0
        y = round(bar_bottom - fraction * (bar_bottom - bar_top))
        value = low + fraction * (high - low)
        draw.line(((bar_right, y), (bar_right + 8, y)), fill=(40, 40, 40), width=2)
        draw.text((bar_right + 12, y - 11), f"{value:.3f}", fill=(30, 30, 30), font=small)

    draw.text(
        (48, 936),
        "Constant inflow | Resistance outlets | Total physical flow: 30.699 mL/s",
        fill=(65, 65, 65),
        font=small,
    )
    image.save(path, "PNG")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--centerline", type=Path, default=DEFAULT_CENTERLINE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--nodes", type=Path, default=DEFAULT_NODES)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    args = parser.parse_args()

    model = json.loads(args.model.resolve().read_text())
    results = read_vessel_results(args.result.resolve())
    data = read_polydata(args.centerline.resolve())
    pressure, flow, vessel_ids, segment_ids = project(model, results, data)
    write_polydata(args.output.resolve(), data)
    write_nodes(args.nodes.resolve(), data, pressure, flow, vessel_ids, segment_ids)
    render(args.figure.resolve(), data, pressure)

    finite_pressure = [value for value in pressure if math.isfinite(value)]
    print(f"Wrote centerline: {args.output.resolve()}")
    print(f"Wrote node table: {args.nodes.resolve()}")
    print(f"Wrote figure: {args.figure.resolve()}")
    print(
        f"Projected {data.GetNumberOfPoints()} points; pressure range "
        f"{min(finite_pressure):.6f}-{max(finite_pressure):.6f} mmHg"
    )


if __name__ == "__main__":
    main()
