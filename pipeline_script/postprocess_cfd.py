#!/usr/bin/env python3
"""Postprocess final-cycle svMultiPhysics CFD metrics and convergence."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import fmean


MMHG_PER_DYN_CM2 = 1.0 / 1333.22
RELATIVE_TOLERANCE = 0.05


def rows(path: Path) -> list[dict[str, float]]:
    lines = path.read_text().splitlines()
    header_index = next(i for i, line in enumerate(lines) if line.startswith("step"))
    header = lines[header_index].split()
    return [dict(zip(header, map(float, line.split()))) for line in lines[header_index + 1:] if line.split()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--config", type=Path, default=Path("step0_config.json"))
    parser.add_argument(
        "--results",
        type=Path,
        default=Path("."),
        help="Directory containing B_NS_Pressure_average.txt and B_NS_Velocity_flux.txt; defaults to the case directory.",
    )
    parser.add_argument("--output", type=Path, default=Path("."))
    args = parser.parse_args()
    case = args.case_directory.resolve()
    config_path = args.config if args.config.is_absolute() else case / args.config
    results = args.results if args.results.is_absolute() else case / args.results
    output = args.output if args.output.is_absolute() else case / args.output
    output.mkdir(parents=True, exist_ok=True)
    config = json.loads(config_path.read_text())
    pressure = {row["step"]: row for row in rows(results / "B_NS_Pressure_average.txt")}
    flow = {row["step"]: row for row in rows(results / "B_NS_Velocity_flux.txt")}
    common = sorted(set(pressure) & set(flow))
    count = int(config["simulation"].get("steps_per_cycle", 1000))
    if len(common) < count:
        raise ValueError("Less than one complete cardiac cycle is available.")
    final_steps = common[-count:]
    fontan = config["metrics"]["fontan_pressure_face"]
    lpa, rpa = config["outlets"]["lpa"], config["outlets"]["rpa"]
    t0 = pressure[final_steps[0]]["time"]
    final = [{"step": int(step), "time_in_cycle_s": pressure[step]["time"] - t0, "fontan_pressure_mmHg": pressure[step][fontan] * MMHG_PER_DYN_CM2, "lpa_flow_ml_s": sum(flow[step][face] for face in lpa), "rpa_flow_ml_s": sum(flow[step][face] for face in rpa)} for step in final_steps]
    mean_pressure = fmean(item["fontan_pressure_mmHg"] for item in final)
    pressure_amplitude = max(item["fontan_pressure_mmHg"] for item in final) - min(item["fontan_pressure_mmHg"] for item in final)
    mean_lpa, mean_rpa = fmean(item["lpa_flow_ml_s"] for item in final), fmean(item["rpa_flow_ml_s"] for item in final)
    total = mean_lpa + mean_rpa
    target_pressure = float(config["targets"]["mean_fontan_pressure_mmhg"])
    target_lpa = total * float(config["targets"]["lpa_fraction_percent"]) / 100.0
    target_rpa = total - target_lpa
    tolerance = RELATIVE_TOLERANCE
    table = []
    for label, value, target, unit in (("Mean Fontan pressure", mean_pressure, target_pressure, "mmHg"), ("Mean LPA flow", mean_lpa, target_lpa, "mL/s"), ("Mean RPA flow", mean_rpa, target_rpa, "mL/s")):
        error = abs(value / target - 1.0)
        table.append({"metric": label, "current": value, "target": target, "difference": value - target, "relative_error": error, "status": "PASS" if error <= tolerance else "FAIL", "units": unit})
    if "fontan_pressure_amplitude_mmhg" in config["targets"]:
        amplitude_target = float(config["targets"]["fontan_pressure_amplitude_mmhg"])
        error = abs(pressure_amplitude / amplitude_target - 1.0)
        table.append({"metric": "Fontan pressure amplitude", "current": pressure_amplitude, "target": amplitude_target, "difference": pressure_amplitude - amplitude_target, "relative_error": error, "status": "PASS" if error <= tolerance else "FAIL", "units": "mmHg"})
    cycles = []
    for start in range(0, len(common) - count + 1, count):
        steps = common[start:start + count]
        cycles.append({"cycle": len(cycles) + 1, "mean_fontan_pressure_mmHg": fmean(pressure[step][fontan] * MMHG_PER_DYN_CM2 for step in steps)})
    for index, item in enumerate(cycles): item["change_from_previous_cycle_mmHg"] = None if index == 0 else item["mean_fontan_pressure_mmHg"] - cycles[index - 1]["mean_fontan_pressure_mmHg"]
    for path, values in ((output / "clinical_target_comparison.csv", table),):
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0])); writer.writeheader(); writer.writerows(values)
    markdown = ["# Initial CFD Clinical Comparison", "", "| Metric | Current | Target | Difference | Relative error | Status |", "| --- | ---: | ---: | ---: | ---: | --- |"]
    for item in table: markdown.append(f"| {item['metric']} | {item['current']:.3f} {item['units']} | {item['target']:.3f} {item['units']} | {item['difference']:+.3f} {item['units']} | {item['relative_error'] * 100:.2f}% | {item['status']} |")
    last_change = cycles[-1]["change_from_previous_cycle_mmHg"]
    markdown.extend(["", "## Convergence", ""])
    if last_change is None:
        markdown.append("Only one complete cardiac cycle is available; cycle-to-cycle Fontan-pressure convergence cannot be assessed.")
    elif abs(last_change) < 0.1:
        markdown.append(f"The last-cycle mean Fontan pressure change is **{last_change:+.3f} mmHg** (smaller than 0.1 mmHg); the simulation is converged.")
    else:
        markdown.append(f"The last-cycle mean Fontan pressure change is **{last_change:+.3f} mmHg** (not smaller than 0.1 mmHg); the simulation is not converged.")
    (output / "clinical_target_comparison.md").write_text("\n".join(markdown) + "\n")
    print("\n".join(markdown))


if __name__ == "__main__":
    main()
