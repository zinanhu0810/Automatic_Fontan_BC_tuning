# P15 3D-to-0D Generation Report

- Status: **PASS**
- Centerline root: `cap_fontan_x`
- Geometry-derived vessels: 29
- Junctions: 6
- Prescribed-flow inlets: 5
- RCR outlets: 19
- Mean prescribed inflow: 30.6988 mL/s
- Final svZeroD file: `P15_0d.json`

## Baseline solver check

- Status: **PASS**
- Mean RCR outflow: 30.6982 mL/s
- Mean external-flow residual: 4.63748e-06 mL/s
- Relative external-flow residual: 1.13693e-07

## Multi-inlet convention

The root cap is attached to the root vessel inlet and uses positive flow into the network. Other physical inlets terminate branches of the rooted centerline tree and therefore retain the negative native CFD surface-flux sign.
