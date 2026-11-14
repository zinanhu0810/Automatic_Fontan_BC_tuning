# P18 3D-to-0D generation report

- Status: **PASS**
- Root cap: `cap_AZG`
- Geometry-derived vessels: 10
- Junctions: 3
- FLOW boundaries: 6
- RCR boundaries: 2
- Resistance boundaries: 0
- Mean prescribed inflow: 54.2403 mL/s
- Final model: `P18_0d.json`
- Solver result: `P18_0d.csv`

## Flow sign convention

The XML surface flux is treated as positive outward from the 3D domain. The selected root cap is reversed to positive into the rooted 0D network. Other prescribed-flow caps retain their XML sign because they attach at terminal vessel outlets.
