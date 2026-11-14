# P3 3D-to-0D generation report

- Status: **PASS**
- Root cap: `cap_Fontan`
- Geometry-derived vessels: 29
- Junctions: 7
- FLOW boundaries: 5
- RCR boundaries: 18
- Resistance boundaries: 0
- Mean prescribed inflow: 65.0786 mL/s
- Final model: `P3_0d.json`
- Solver result: `P3_0d.csv`

## Flow sign convention

The XML surface flux is treated as positive outward from the 3D domain. The selected root cap is reversed to positive into the rooted 0D network. Other prescribed-flow caps retain their XML sign because they attach at terminal vessel outlets.
