# P1 Step 3: Calibrate 0D Vessel Structure to 3D CFD

This stage reads the geometric model from `../02_generate_0d/rom/P1_0d.json`
and fits its vessel properties to the final complete P1 3D CFD cycle
(steps 8699–9698). The targets are mean Fontan pressure, the aggregate LPA/RPA
flow split, and Fontan-pressure amplitude.

All seven prescribed FLOW boundary conditions and all 17 RCR BCs are retained
exactly. Calibration changes geometric vessels only:

- global vessel resistance;
- relative terminal-vessel resistance for the LPA and RPA outlet groups;
- global vessel inductance.

Geometric vessel capacitance is set to zero to match the rigid 3D walls. RCR
capacitances remain unchanged. A fitted metric passes when its relative error
against the final-cycle 3D value is within `0.1%`.

```bash
python3 calibrate_p1_unsteady.py
```

The calibrated 0D model, tables, and final-cycle 3D/0D Fontan-pressure overlay
are written to `unsteady/`.
