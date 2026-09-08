# Closed-loop genBC and svZeroD comparison

This workflow translates the original closed-loop Fontan LPN from genBC into a
corresponding svZeroD JSON model and compares the results without modifying the
source parameters.

## Contents

- `passive_fluid_truncate`: legacy genBC input files and Fortran source needed
  to audit the original equations, parameters, and initial conditions
- `sv0D_LPN/build_genbc_sv0d_case.py`: generates the svZeroD JSON and comparison
  notebook
- `sv0D_LPN/compare_genBC_sv0D.ipynb`: generated comparison of genBC Python,
  Python hybrid, and `svzerodsolver`
- `sv0D_LPN/*summary*`: compact tables and slide-ready SVG summaries from the
  reference run
- `sv0D_LPN/*result.csv` and `comparison_timeseries.csv`: compact reference
  time series used by the notebook

The large `Fontan_genBC_1to1_sv0D_result.csv` file is generated locally and is
not tracked.

## Generate the model and notebook

From this directory:

```bash
python3 sv0D_LPN/build_genbc_sv0d_case.py
```

This writes:

- `sv0D_LPN/Fontan_genBC_1to1_sv0D.json`
- `sv0D_LPN/compare_genBC_sv0D.ipynb`

Run the notebook from this directory so its relative paths resolve correctly:

```bash
jupyter lab sv0D_LPN/compare_genBC_sv0D.ipynb
```

The notebook expects `svzerodsolver` on `PATH`. The generated model uses
`LinearElastanceChamber` and `PiecewiseValve`, matching the chamber and valve
formulations used by the Python hybrid comparison.

