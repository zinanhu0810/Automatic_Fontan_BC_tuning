#!/usr/bin/env python3
"""Calculate final-window Fontan pressure and pulmonary flow split."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import fmean, pstdev


HERE = Path(__file__).resolve().parent
DEFAULT_RESULTS = HERE / "48-procs"
MMHG_PER_DYN_CM2 = 1.0 / 1333.22

INLET_FACES = (
    "cap_fontan_x",
    "cap_L_HV",
    "cap_M_HV",
    "cap_R_HV",
    "cap_rt_svc_x",
)


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
            raise ValueError(
                f"Expected {len(names)} fields but found {len(fields)} in {path}"
            )
        rows.append({name: float(value) for name, value in zip(names, fields)})
    if not rows:
        raise ValueError(f"No boundary-integral samples found in {path}")
    return rows


def describe(values: list[float]) -> dict[str, float]:
    return {
        "mean": fmean(values),
        "standard_deviation": pstdev(values),
        "minimum": min(values),
        "maximum": max(values),
        "range": max(values) - min(values),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--last", type=int, default=100)
    parser.add_argument(
        "--output", type=Path, default=HERE / "steady_metrics_summary.json"
    )
    parser.add_argument(
        "--window-csv", type=Path, default=HERE / "steady_metrics_last_100.csv"
    )
    args = parser.parse_args()
    if args.last <= 0:
        raise ValueError("--last must be positive")

    results = args.results.resolve()
    pressure_rows = read_boundary_integral(results / "B_NS_Pressure_average.txt")
    flow_rows = read_boundary_integral(results / "B_NS_Velocity_flux.txt")
    if len(pressure_rows) != len(flow_rows):
        raise ValueError("Pressure and flow files contain different sample counts")

    number = min(args.last, len(flow_rows))
    pressure_rows = pressure_rows[-number:]
    flow_rows = flow_rows[-number:]
    for pressure, flow in zip(pressure_rows, flow_rows):
        if pressure["step"] != flow["step"] or pressure["time"] != flow["time"]:
            raise ValueError("Pressure and flow samples are not aligned")

    names = set(flow_rows[0])
    left_faces = sorted(name for name in names if name.startswith("cap_l_pa_"))
    right_faces = sorted(name for name in names if name.startswith("cap_r_pa_"))
    left_faces.append("cap_pa_x_2")
    right_faces.append("cap_pa_x")
    expected = set(INLET_FACES) | set(left_faces) | set(right_faces)
    missing = expected - names
    if missing:
        raise ValueError(f"Missing expected boundary faces: {sorted(missing)}")

    window_rows = []
    for pressure, flow in zip(pressure_rows, flow_rows):
        left = sum(flow[name] for name in left_faces)
        right = sum(flow[name] for name in right_faces)
        pulmonary = left + right
        inlet = sum(flow[name] for name in INLET_FACES)
        window_rows.append(
            {
                "step": int(flow["step"]),
                "time_s": flow["time"],
                "fontan_pressure_mmHg": pressure["cap_fontan_x"]
                * MMHG_PER_DYN_CM2,
                "lpa_flow_ml_s": left,
                "rpa_flow_ml_s": right,
                "pulmonary_flow_ml_s": pulmonary,
                "lpa_fraction": left / pulmonary,
                "rpa_fraction": right / pulmonary,
                "open_boundary_residual_ml_s": inlet + pulmonary,
            }
        )

    args.window_csv = args.window_csv.with_name(
        f"steady_metrics_last_{number}.csv"
    )
    with args.window_csv.resolve().open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=window_rows[0].keys())
        writer.writeheader()
        writer.writerows(window_rows)

    metrics = {
        "source_directory": str(results),
        "number_of_samples": number,
        "first_step": window_rows[0]["step"],
        "last_step": window_rows[-1]["step"],
        "first_time_s": window_rows[0]["time_s"],
        "last_time_s": window_rows[-1]["time_s"],
        "pressure_units": "mmHg",
        "flow_units": "mL/s",
        "fontan_pressure": describe(
            [row["fontan_pressure_mmHg"] for row in window_rows]
        ),
        "lpa_flow": describe([row["lpa_flow_ml_s"] for row in window_rows]),
        "rpa_flow": describe([row["rpa_flow_ml_s"] for row in window_rows]),
        "total_pulmonary_flow": describe(
            [row["pulmonary_flow_ml_s"] for row in window_rows]
        ),
        "lpa_fraction": describe([row["lpa_fraction"] for row in window_rows]),
        "rpa_fraction": describe([row["rpa_fraction"] for row in window_rows]),
        "open_boundary_residual": describe(
            [row["open_boundary_residual_ml_s"] for row in window_rows]
        ),
        "left_pulmonary_faces": left_faces,
        "right_pulmonary_faces": right_faces,
        "flow_sign_convention": "positive outward from the 3D domain",
    }
    args.output.resolve().write_text(json.dumps(metrics, indent=2) + "\n")

    print(f"Window: steps {metrics['first_step']}-{metrics['last_step']} ")
    print(f"        time {metrics['first_time_s']:.3f}-{metrics['last_time_s']:.3f} s")
    print(
        f"Mean Fontan pressure: {metrics['fontan_pressure']['mean']:.6f} mmHg"
    )
    print(f"Mean LPA flow: {metrics['lpa_flow']['mean']:.6f} mL/s")
    print(f"Mean RPA flow: {metrics['rpa_flow']['mean']:.6f} mL/s")
    print(
        f"Pulmonary split: L {100.0 * metrics['lpa_fraction']['mean']:.6f}% : "
        f"R {100.0 * metrics['rpa_fraction']['mean']:.6f}%"
    )
    print(
        f"Mean open-boundary residual: "
        f"{metrics['open_boundary_residual']['mean']:.6e} mL/s"
    )
    print(f"Wrote: {args.output.resolve()}")
    print(f"Wrote: {args.window_csv.resolve()}")


if __name__ == "__main__":
    main()
