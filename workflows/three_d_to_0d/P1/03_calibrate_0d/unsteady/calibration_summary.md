# P1 Unsteady 0D Calibration

The 0D geometric vessels were fit to the final 3D cycle. All prescribed FLOW and RCR boundary conditions are unchanged; vessel capacitances are zero for the rigid-wall comparison.

| Metric | 3D target | Rigid 0D baseline | Calibrated 0D | Relative error | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| Mean Fontan pressure | 11.3588 mmHg | 11.0377 mmHg | 11.3585 mmHg | 0.003% | PASS |
| LPA flow split | 47.7379 % | 47.9416 % | 47.7473 % | 0.020% | PASS |
| Fontan pressure amplitude | 1.5053 mmHg | 1.3266 mmHg | 1.5042 mmHg | 0.075% | PASS |

## Vessel-only scales

- Global vessel resistance: `9.4`
- LPA/RPA terminal-vessel resistance balance: `1.005`
- Global vessel inductance: `1.221`
- Geometric vessel capacitance: `0.0`
- Acceptance criterion: relative error `<= 0.1%` for each fitted metric.
