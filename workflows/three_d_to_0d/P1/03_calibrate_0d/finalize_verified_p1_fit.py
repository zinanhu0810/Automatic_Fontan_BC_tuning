#!/usr/bin/env python3
"""Write reproducible P1 Step 3 artifacts for the verified local optimum."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import calibrate_p1_unsteady as calibration


SCALES = (9.4, 1.005, 1.221)
OUTPUT = calibration.OUTPUT


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    solver, target = calibration.find_solver(), calibration.target_metrics()
    baseline_model = calibration.prepare_model(calibration.SOURCE_MODEL, (1.0, 1.0, 1.0))
    calibrated_model = calibration.prepare_model(calibration.SOURCE_MODEL, SCALES)
    baseline_json, baseline_csv = OUTPUT / "P1_0d_rigid_baseline.json", OUTPUT / "P1_0d_rigid_baseline.csv"
    calibrated_json, calibrated_csv = OUTPUT / "P1_0d_calibrated.json", OUTPUT / "P1_0d_calibrated.csv"
    calibration.run_solver(solver, baseline_model, baseline_json, baseline_csv)
    calibration.run_solver(solver, calibrated_model, calibrated_json, calibrated_csv)
    baseline = calibration.result_metrics(baseline_model, baseline_csv)
    calibrated = calibration.result_metrics(calibrated_model, calibrated_csv)
    rows = []
    for label, key, unit, scale in (
        ("Mean Fontan pressure", "mean_fontan_pressure_mmHg", "mmHg", 1.0),
        ("LPA flow split", "lpa_split_fraction", "%", 100.0),
        ("Fontan pressure amplitude", "fontan_pressure_amplitude_mmHg", "mmHg", 1.0),
    ):
        rows.append({"metric": label, "unit": unit, "3D_target": target[key] * scale, "0D_baseline": baseline[key] * scale, "0D_calibrated": calibrated[key] * scale, "calibrated_relative_error": abs(calibrated[key] / target[key] - 1.0)})
    calibration.write_csv(OUTPUT / "calibration_metrics.csv", rows)
    calibration.write_markdown_summary(OUTPUT / "calibration_summary.md", target, baseline, calibrated, SCALES)
    calibration.write_fontan_waveform_svg(OUTPUT / "final_cycle_fontan_pressure_3d_vs_0d.svg", calibrated_model, calibrated_csv)
    (OUTPUT / "calibration_targets.json").write_text(json.dumps(target, indent=2) + "\n")
    summary = {
        "target_source": str(calibration.THREE_D_RESULTS),
        "calibration_method": "verified local response scan; vessel-only fit against final 3D cycle",
        "parameters": {"global_vessel_resistance_scale": SCALES[0], "lpa_to_rpa_terminal_resistance_scale": SCALES[1], "global_vessel_inductance_scale": SCALES[2], "geometric_vessel_capacitance": 0.0, "boundary_conditions_preserved": True},
        "baseline": baseline, "calibrated": calibrated,
        "objective": calibration.objective(calibrated, target), "evaluations": "verified local response scan",
        "acceptance_criterion": {"relative_error_tolerance": calibration.MATCH_RELATIVE_TOLERANCE},
        "acceptance": {key: abs(calibrated[key] / target[key] - 1.0) <= calibration.MATCH_RELATIVE_TOLERANCE for key in ("mean_fontan_pressure_mmHg", "lpa_split_fraction", "fontan_pressure_amplitude_mmHg")},
    }
    (OUTPUT / "calibration_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
