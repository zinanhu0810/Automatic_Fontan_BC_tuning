# P15 Step 3: Steady 3D-to-0D Calibration

This folder calibrates the P15 geometric 0D model against the converged steady 3D CFD result.

## Method

1. Average the final 100 samples of the 3D boundary-integral pressure and flow files.
2. Recover each directed 0D vessel flow from the measured boundary fluxes and junction continuity.
3. Probe the final 3D pressure field onto each 0D vessel endpoint on the ROM centerline.
4. Build the complete `y` and `dy` observation dictionaries expected by `svzerodcalibrator`.
5. Run `svzerodcalibrator`, then run `svzerodsolver` with the calibrated model.
6. Compare mean Fontan pressure and pulmonary flow split with the steady 3D targets.

The current CFD source is steady. Therefore, only `R_poiseuille` is identifiable. Vessel `C`, `L`, and `stenosis_coefficient` are preserved exactly from the baseline geometric model. An unsteady 3D run is required before calibrating branch-specific `L` or `C`.

The C++ calibrator is unconstrained. If a local centerline static-pressure recovery produces a nonpositive fitted resistance, the simulation-ready JSON restores that branch's original positive geometric resistance. The raw calibrator output and both resistance columns remain available for audit.

After that physicality repair, the driver applies one deterministic resistance correction to the prescribed-flow root Fontan branch. This aligns mean Fontan pressure without changing the pulmonary flow split; the correction is recorded in `calibration_summary.json`.

## Run

```bash
python3 calibrate_p15.py
```

## Main outputs

- `P15_steady_calibration_input.json`: full calibrator input with 3D observations
- `P15_0d_calibrated.json`: simulation-ready calibrated 0D model
- `P15_0d_calibrated.csv`: calibrated standalone 0D result
- `calibration_metrics.csv`: 3D target vs baseline and calibrated 0D metrics
- `vessel_parameter_comparison.csv`: branch-level parameter changes
- `calibration_summary.png`: compact visual validation
- `calibration_summary.json`: machine-readable run summary
- `P15_0d_calibrated_pressure_centerlines.vtp`: calibrated pressure and flow projected onto the ROM centerline
- `P15_0d_calibrated_pressure_centerline_nodes.csv`: projected values at all centerline points
- `P15_0d_calibrated_pressure_centerlines.png`: two-view pressure preview

The script is adapted to the observation naming and calibrator invocation used in Natalia Rubio's `learn_lpns/zerod_calibration` workflow.

## Current P15 result

| Metric | Steady 3D target | Baseline 0D | Calibrated 0D | Final relative error |
|---|---:|---:|---:|---:|
| Mean Fontan pressure | 15.623 mmHg | 15.029 mmHg | 15.623 mmHg | <0.001% |
| LPA flow | 12.8467 mL/s | 12.9385 mL/s | 12.8502 mL/s | 0.027% |
| RPA flow | 17.8521 mL/s | 17.7603 mL/s | 17.8486 mL/s | 0.019% |
| LPA split | 41.8477% | 42.1465% | 41.8588% | 0.027% |

The raw C++ fit reached its 2,000-iteration cap because the steady problem contains unidentifiable dynamic directions. Six raw local resistances were negative and were restored to their positive geometric values in the applied JSON. Use `P15_0d_calibrated_raw.json` only for diagnostics; use `P15_0d_calibrated.json` for simulation.
