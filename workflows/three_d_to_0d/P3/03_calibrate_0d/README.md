# P3 Step 3: Calibrate 0D Vessel Structure to 3D CFD

This stage fits the geometric 0D vessels to the final P3 3D cycle (steps
5001–6000): mean Fontan pressure, aggregate LPA/RPA split, and Fontan-pressure
amplitude. All five prescribed FLOW BCs and all 18 RCR BCs remain unchanged.

Only vessel resistance, relative LPA/RPA terminal-vessel resistance, and vessel
inductance are varied. Geometric vessel capacitance is zero for the rigid-wall
3D comparison, while RCR capacitances are unchanged. The acceptance threshold
is `0.1%` relative error for every fitted metric.
