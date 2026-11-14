# P18 Step 4: Tune Clinical Boundary Conditions

This stage starts with the rigid geometric 0D model calibrated in Step 3 and
tunes only the two pulmonary RCR boundary conditions to the P18 clinical
targets. The model uses the same prescribed inflows and vessel R/L values as
Step 3.

The two Nelder-Mead controls are in log space:

- common RCR resistance scale, applied to both `Rp` and `Rd` values;
- LPA/RPA resistance bias, which changes the pulmonary flow split.

RCR capacitance, distal pressure, all prescribed FLOW BCs, and all geometric
vessel parameters are unchanged. The clinical pressure criterion is 5%
relative error. Flow criteria compare LPA and RPA flow separately, each within
5% of its target flow value.

Run:

```bash
python3 tune_p18_bc_nelder_mead.py
```

Outputs are written to `unsteady/`. The model for the next coupled 3D run is
`unsteady/P18_0d_tuned_bc.json`.
