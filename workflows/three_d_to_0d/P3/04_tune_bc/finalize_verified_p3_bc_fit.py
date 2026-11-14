#!/usr/bin/env python3
"""Write P3 Step 4 artifacts for the verified RCR tuning point."""

from __future__ import annotations

import json
import math

import tune_p3_bc_nelder_mead as tuning


POSITION = (math.log(0.893), -0.03)


def main() -> None:
    tuning.OUTPUT.mkdir(parents=True, exist_ok=True)
    source = json.loads(tuning.SOURCE_MODEL.read_text())
    targets = json.loads(tuning.TARGET_FILE.read_text())
    solver = tuning.find_solver()
    baseline_model = tuning.apply_rcr_scales(source, (0.0, 0.0))
    tuned_model = tuning.apply_rcr_scales(source, POSITION)
    baseline_json, baseline_csv = tuning.OUTPUT / "P3_0d_rcr_baseline.json", tuning.OUTPUT / "P3_0d_rcr_baseline.csv"
    tuned_json, tuned_csv = tuning.OUTPUT / "P3_0d_tuned_bc.json", tuning.OUTPUT / "P3_0d_tuned_bc.csv"
    tuning.run_solver(solver, baseline_model, baseline_json, baseline_csv)
    tuning.run_solver(solver, tuned_model, tuned_json, tuned_csv)
    baseline, tuned = tuning.metrics_from_result(baseline_model, baseline_csv), tuning.metrics_from_result(tuned_model, tuned_csv)
    total_flow = tuned["mean_lpa_flow_ml_s"] + tuned["mean_rpa_flow_ml_s"]
    target_lpa, target_rpa = targets["lpa_flow_fraction"] * total_flow, targets["rpa_flow_fraction"] * total_flow
    errors = {"fontan_pressure": abs(tuned["mean_fontan_pressure_mmHg"] / targets["mean_fontan_pressure_mmHg"] - 1.0), "lpa_flow": abs(tuned["mean_lpa_flow_ml_s"] / target_lpa - 1.0), "rpa_flow": abs(tuned["mean_rpa_flow_ml_s"] / target_rpa - 1.0)}
    scales = {"common_rcr_resistance_scale": math.exp(POSITION[0]), "lpa_rcr_resistance_scale": math.exp(POSITION[0] + POSITION[1]), "rpa_rcr_resistance_scale": math.exp(POSITION[0] - POSITION[1])}
    summary = {"source_model": str(tuning.SOURCE_MODEL), "method": "verified local Nelder-Mead response point", "targets": targets, "baseline": baseline, "tuned": tuned, "scales": scales, "relative_errors": errors, "acceptance": {"fontan_pressure": errors["fontan_pressure"] <= targets["mean_fontan_pressure_relative_tolerance"], "lpa_flow": errors["lpa_flow"] <= targets["pulmonary_flow_relative_tolerance"], "rpa_flow": errors["rpa_flow"] <= targets["pulmonary_flow_relative_tolerance"]}}
    (tuning.OUTPUT / "tuning_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = ["# P3 Clinical RCR BC Tuning", "", "| Metric | Clinical target | Calibrated-geometry baseline | Tuned 0D | Relative error | Status |", "| --- | ---: | ---: | ---: | ---: | --- |"]
    for label, key, target, error_key in (("Mean Fontan pressure", "mean_fontan_pressure_mmHg", targets["mean_fontan_pressure_mmHg"], "fontan_pressure"), ("LPA flow", "mean_lpa_flow_ml_s", target_lpa, "lpa_flow"), ("RPA flow", "mean_rpa_flow_ml_s", target_rpa, "rpa_flow")):
        unit = "mmHg" if error_key == "fontan_pressure" else "mL/s"
        lines.append(f"| {label} | {target:.4f} {unit} | {baseline[key]:.4f} {unit} | {tuned[key]:.4f} {unit} | {errors[error_key] * 100.0:.3f}% | {'PASS' if summary['acceptance'][error_key] else 'FAIL'} |")
    lines.extend(["", "## Tuned RCR resistance scales", "", f"- Common scale: `{scales['common_rcr_resistance_scale']:.6g}`", f"- LPA RCR scale: `{scales['lpa_rcr_resistance_scale']:.6g}`", f"- RPA RCR scale: `{scales['rpa_rcr_resistance_scale']:.6g}`", "- RCR capacitance, distal pressure, prescribed inflows, and all geometric vessel parameters are unchanged."])
    (tuning.OUTPUT / "tuning_summary.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
