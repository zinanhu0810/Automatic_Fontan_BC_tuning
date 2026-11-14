# P18 Step 3: Calibrate Vessels to Unsteady 3D CFD

This stage fits the geometric 0D vessels to the final P18 3D cardiac cycle.
The targets are mean Fontan pressure, LPA/RPA flow split, and Fontan-pressure
amplitude.

The six prescribed FLOW BCs and two RCR BCs are preserved exactly from the
Step 2 P18 model. Calibration changes only geometric vessel parameters:

- a global vessel resistance scale;
- a relative LPA/RPA terminal-vessel resistance scale;
- a global vessel inductance scale.

Geometric vessel capacitance is set to zero to represent the rigid 3D walls.
The downstream RCR capacitances remain unchanged.

Run:

```bash
python3 calibrate_p18_unsteady.py
```

Outputs are written to `unsteady/`.
