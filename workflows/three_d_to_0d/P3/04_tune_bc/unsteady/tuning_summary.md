# P3 Clinical RCR BC Tuning

| Metric | Clinical target | Calibrated-geometry baseline | Tuned 0D | Relative error | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| Mean Fontan pressure | 9.0000 mmHg | 9.9614 mmHg | 8.9840 mmHg | 0.178% | PASS |
| LPA flow | 34.8815 mL/s | 34.0842 mL/s | 35.3411 mL/s | 1.318% | PASS |
| RPA flow | 30.1960 mL/s | 30.9933 mL/s | 29.7364 mL/s | 1.522% | PASS |

## Tuned RCR resistance scales

- Common scale: `0.892649`
- LPA RCR scale: `0.852692`
- RPA RCR scale: `0.934478`
- RCR capacitance, distal pressure, prescribed inflows, and all geometric vessel parameters are unchanged.
