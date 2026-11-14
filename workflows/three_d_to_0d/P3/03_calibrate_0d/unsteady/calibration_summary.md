# P3 Unsteady 0D Calibration

The 0D geometric vessels were fit to the final 3D cycle. All prescribed FLOW and RCR boundary conditions are unchanged; vessel capacitances are zero for the rigid-wall comparison.

| Metric | 3D target | Rigid 0D baseline | Calibrated 0D | Relative error | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| Mean Fontan pressure | 9.9607 mmHg | 9.0588 mmHg | 9.9614 mmHg | 0.007% | PASS |
| LPA flow split | 52.3608 % | 53.4561 % | 52.3747 % | 0.027% | PASS |
| Fontan pressure amplitude | 2.3330 mmHg | 1.0672 mmHg | 2.3315 mmHg | 0.064% | PASS |

## Vessel-only scales

- Global vessel resistance: `16.37`
- LPA/RPA terminal-vessel resistance balance: `0.533`
- Global vessel inductance: `2.879`
- Geometric vessel capacitance: `0.0`
- Acceptance criterion: relative error `<= 0.1%` for each fitted metric.
