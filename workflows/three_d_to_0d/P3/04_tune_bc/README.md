# P3 Step 4: Tune Clinical RCR Boundary Conditions

This stage starts from the Step 3 calibrated rigid P3 geometry and uses bounded
Nelder-Mead tuning of two log-space controls: a common RCR resistance scale and
an LPA/RPA RCR bias. `Rp` and `Rd` are scaled together within each outlet; all
RCR capacitances, distal pressures, prescribed FLOW BCs, and geometric vessel
properties are unchanged.

Clinical targets: mean Fontan pressure `9.0 mmHg`; LPA/RPA flow fractions
`53.6% / 46.4%`. The pressure and each pulmonary flow target must be within
5% relative error.

```bash
python3 tune_p3_bc_nelder_mead.py
```
