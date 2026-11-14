# P1 3D-to-0D generation report

- Status: **PASS**
- Root cap: `cap_Fontan`
- Geometry-derived vessels: 31
- Junctions: 8
- FLOW boundaries: 7
- RCR boundaries: 17
- Resistance boundaries: 0
- Mean prescribed inflow: 49.8777 mL/s
- Final model: `P1_0d.json`
- Solver result: `P1_0d.csv`

## Flow sign convention

The XML surface flux is treated as positive outward from the 3D domain. The selected root cap is reversed to positive into the rooted 0D network. Other prescribed-flow caps retain their XML sign because they attach at terminal vessel outlets.
