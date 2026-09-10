# Automatic Fontan Tuning

Reproducible examples for translating and calibrating Fontan boundary conditions
with SimVascular workflow.


## Planned BC Tuning Pipeline

The same workflow is intended to support both open- and closed-loop 0D
boundary conditions; the current focus is open-loop cases.

### 1. Run an initial 3D CFD simulation

Run one steady or unsteady 3D simulation using analytical boundary-condition
parameter estimates derived from clinical data. The clinical targets are Fontan
pressure and pulmonary flow split. The resistance/RCR value at each outlet is 
calculated based on cap area (Murray's law).

If initial 3D CFD satisiyf the criteria (within 5%), the boundary conditions are 
ready. Some times, due to the branching, junction resistance, several iterations
may required and used 0D method (step 2-5) will accelerate the whole tuning procedure.


### 2. Generate a patient-specific 0D model from the 3D geometry

Generate a geometric 0D model for each patient using SimVascular's ROM
simulation code. This is being adapted here for command-line use from
[`sv_rom_simulation`](https://github.com/SimVascular/SimVascular/tree/master/Python/site-packages/sv_rom_simulation).
The approach should be evaluated on increasingly complex 3D models.

### 3. Calibrate vessel RCL parameters against the initial 3D simulation

The P15 implementation is in
[`workflows/three_d_to_0d/P15/step3_calibration`](workflows/three_d_to_0d/P15/step3_calibration).
`calibrate_p15.py` averages the final 3D pressure and flow samples, reconstructs
directed vessel flows from the boundary fluxes and junction continuity, and
probes the 3D pressure field at vessel endpoints on the ROM centerline. It then
builds the `y` and `dy` observations required by `svzerodcalibrator`, calibrates
the branch parameters, runs the calibrated model with `svzerodsolver`, and
compares Fontan pressure and pulmonary flow split against the 3D reference.

For a steady 3D reference, only vessel resistance is identifiable. The script preserves
those geometric values, repairs nonpositive fitted resistances using the
baseline model, and records every applied correction in the summary outputs.
The folder includes the calibrator input, raw and simulation-ready models,
parameter comparison, error metrics, and pressure projected onto the
centerline.

### 4. Tune physiological boundary conditions in 0D


### 5. Validate with the 3D model in the loop


## Software

Install the Python packages in `requirements.txt`. The workflows additionally
require SimVascular's `sv_rom_simulation` runtime and the `svzerodsolver` and
`svzerodcalibrator` executables on `PATH` for the stages that invoke them.

Raw 3D processor output, restart files, compiled binaries, solver logs, and
large reproducible time-series files are intentionally excluded. Compact
reference metrics and visualization files are retained.

## Other scripts

### Closed-loop genBC to svZeroD

`workflows/closed_loop_0d` contains the legacy Fortran genBC source, the
one-to-one svZeroD model builder, and a notebook comparing three formulations:

- original genBC equations implemented in Python
- a Python hybrid using the svZeroD chamber and valve laws
- the full `svzerodsolver` JSON model

The accompanying `sv0D_LPN_Tuner` runs bounded automatic tuning of chamber,
valve, and optional vascular parameters against a configurable target vector.

See [`workflows/closed_loop_0d/README.md`](workflows/closed_loop_0d/README.md).

### Patient-specific 3D to 0D

`workflows/three_d_to_0d` contains a command-line wrapper around
SimVascular's `sv_rom_simulation` package. The P15 example includes the surface
model, centerlines, generated 0D network, a resistance-only steady case, steady
3D validation summaries, branch-resistance calibration, and pressure projected
onto the centerline.

See [`workflows/three_d_to_0d/README.md`](workflows/three_d_to_0d/README.md).
