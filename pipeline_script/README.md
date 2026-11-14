# Step 0: Initial CFD Setup

This folder contains portable, command-line Step 0 scripts. It replaces the
patient-specific XML creators and postprocessors used in the early pipelines.
It has no GUI dependency and uses only the Python standard library.

## Why Step 0 Makes Sense

Step 0 is an **initial clinical-target simulation**, not a final calibrated
simulation. It uses measured inflow, cap areas, a target mean Fontan pressure,
and target LPA/RPA flow split to create a defensible first RCR estimate.
RCR capacitor pressures are initialized from the periodic prescribed inflow,
which normally avoids the long zero-pressure transient of an all-zero start.

It does not guarantee that the first CFD run meets every target: 3D geometric
losses and the area-only distribution are deliberately left to the later 0D
structure calibration and BC-tuning stages. A supplied
`fontan_pressure_amplitude_mmhg` is reported by the postprocessor, but it is
not used to analytically size the initial RCR capacitance. Capacitance is a
dynamic parameter and should be calibrated after the initial run.

## Required Case Inputs

```
case/
  step0_config.json
  cap_info.txt                         # cap name, area [cm^2]
  mesh/
    mesh-complete.mesh.vtu
    mesh-surfaces/cap_*.vtp
    walls*.vtp
  *.flow                               # for each unsteady inlet
```

Every `cap*.vtp` must be declared exactly once as an inlet, LPA outlet, or RPA
outlet. The explicit role lists are intentional: cap names are not reliable
enough to infer topology across patients. `metrics.fontan_pressure_face` is a
representative Fontan/inlet face whose boundary-average pressure is reported.

Copy [`step0_config.example.json`](step0_config.example.json) to the case as
`step0_config.json`, then replace all example cap names and targets. Waveforms
may have any original duration; the generator rescales their time coordinate to
`simulation.cardiac_period_s` while preserving flow values and sign. The solver
expects the established convention: measured inflow stored as negative flux.

## Generate solver.xml

```bash
python3 pipeline_script/generate_xml.py /path/to/case
```

This creates `solver.xml` and one `*_physical.flow` file per inlet. RCR values
use the clinical mean pressure drop and area-weighted outlet flow allocation.
`Rp = R / 11` and `Rd = 10 R / 11` by default, but
`rcr.proximal_resistance_fraction` is configurable.

The script writes phase-aware periodic RCR `Initial_pressure` values directly
into `solver.xml`; a separate initialization script is therefore unnecessary.
For a steady case, replace each inlet `flow_file` with `steady_flow_ml_s`.

## Postprocess CFD Results

After copying `B_NS_Pressure_average.txt` and `B_NS_Velocity_flux.txt` into the
case directory:

```bash
python3 pipeline_script/postprocess_cfd.py /path/to/case
```

For a solver output kept in a separate directory, pass it explicitly:

```bash
python3 pipeline_script/postprocess_cfd.py /path/to/case --results /path/to/48-procs
```

The script uses the final full cardiac cycle and produces clinical comparison
files. The Markdown table ends with a convergence statement based on the
absolute change in mean Fontan pressure between the last two cycles: a change
smaller than `0.1 mmHg` is converged. The default clinical pass criterion is 5%
relative error for mean Fontan pressure and **each pulmonary branch flow**,
rather than five percentage points of flow split.

# Step 1: 3D-to-0D Conversion

[`generate_0D.py`](generate_0D.py) converts the Step 0 surface mesh, cap
surfaces, and `solver.xml` into a geometric svZeroD model. It obtains the
network root from `metrics.fontan_pressure_face` in `step0_config.json` and
restores the prescribed-flow and RCR boundary conditions from `solver.xml`.

The script delegates geometric extraction to the shared
[`sv_rom_0d_toolkit`](../../sv_rom_0d_toolkit), which in turn requires a local
SimVascular installation containing `sv_rom_simulation`. The model must be a
rooted centerline tree; closed loops and nonstandard BC types require an
extension.

```bash
# Check the Step 1 inputs and prepare ROM files without launching SimVascular.
python3 pipeline_script/generate_0D.py /path/to/case --prepare-only

# Generate /path/to/case/0D/CASE_0d.json.
python3 pipeline_script/generate_0D.py /path/to/case \
  --simvascular /Applications/SimVascular.app/Contents/Resources/simvascular
```

Use `--root-surface CAP_NAME` to select a different prescribed-flow root. Step
1 only creates the 0D model; run `svzerodsolver` later when results are needed.

## Step 1 Alternative: VMTK Without SimVascular

[`generate_0d_withoutsim.py`](generate_0d_withoutsim.py) is a standalone
alternative that uses VMTK to extract centerlines and branch geometry directly
from the surface mesh. It writes `0d/centerlines.vtp` and the final
boundary-faithful `0d/CASE_0d.json`, without SimVascular, ROM input files,
raw JSON, or a generation manifest. Geometric vessel capacitance defaults to
zero, consistent with the rigid-wall CFD model; RCR capacitance in `solver.xml`
is retained unchanged.

### Input and Output Layout

The positional argument is the **case directory**, not an individual mesh or
JSON file:

```bash
python3 pipeline_script/generate_0d_withoutsim.py CASE_DIR
```

`CASE_DIR` must contain the Step 0 configuration and solver input, plus the
surface model and cap surfaces named by `step0_config.json`:

```text
CASE_DIR/
  step0_config.json
  solver.xml
  mesh/
    mesh-complete.exterior.vtp
    mesh-surfaces/
      cap_*.vtp
      wall*.vtp
  *.flow                         # Files referenced by solver.xml for FLOW BCs
```

The script creates the output directory automatically:

```text
CASE_DIR/
  0d/
    centerlines.vtp
    CASE_0d.json
```

`CASE` is taken from `case_name` in `step0_config.json`. For example, passing
`P3_BC_tuning_pipeline` creates
`P3_BC_tuning_pipeline/0d/P3_0d.json`. To use a different output directory,
pass `--output-directory`:

```bash
python3 pipeline_script/generate_0d_withoutsim.py CASE_DIR \
  --output-directory 02_generate_0d
```

This repository's `vmtk` Conda environment is x86_64, so invoke it through
Rosetta on this Mac:

```bash
arch -x86_64 /Users/zinanhu/anaconda3/envs/vmtk/bin/python \
  pipeline_script/generate_0d_withoutsim.py /path/to/case
```

Use `--overwrite-centerlines` to regenerate an existing centerline file, or
`--vessel-capacitance VALUE` to override the rigid-wall default. The converter
requires one connected rooted centerline tree and fails if a cap does not map
uniquely to a tree leaf.

# Step 2: Calibrate the 0D Structure

[`calibrate_unsteady.py`](calibrate_unsteady.py) fits the final 0D JSON to the
last complete 3D CFD cycle. Its inputs are `0d/CASE_0d.json`,
`B_NS_Pressure_average.txt`, `B_NS_Velocity_flux.txt`, and `step0_config.json`.
The script sets every geometric vessel capacitance to zero for rigid-wall CFD,
then fits global vessel resistance, LPA/RPA resistance balance, and global
vessel inductance to mean Fontan pressure, pulmonary split, and Fontan-pressure
amplitude. All FLOW/RCR/Resistance boundary conditions remain unchanged.

```bash
python3 pipeline_script/calibrate_unsteady.py /path/to/case
```

It writes `calibrated_0d/CASE_0d_calibrated.json`, a concise Markdown summary,
and a final-cycle 3D-vs-0D Fontan-pressure overlay. Candidate solver CSV files
are temporary and are not kept. The acceptance criterion is `0.1%` relative
error for each fitted metric.

# Step 3: Tune Terminal Boundary Conditions in 0D

[`tune_rcr_bc.py`](tune_rcr_bc.py) uses the calibrated 0D model to tune the
terminal pulmonary RCR or resistance boundary conditions to the clinical targets
in `step0_config.json`. It adjusts a common terminal-resistance scale and an
LPA/RPA resistance-balance scale. Prescribed inflows, distal pressures, and the
calibrated geometric vessel parameters remain unchanged.

The input is the case directory and the calibrated Step 2 model:

```bash
python3 pipeline_script/tune_rcr_bc.py CASE_DIR \
  --model calibrated_0d/CASE_0d_calibrated.json
```

The script requires `svzerodsolver` to be available on `PATH`; otherwise pass
its location explicitly with `--svzerodsolver /path/to/svzerodsolver`.

It creates:

```text
CASE_DIR/
  tuned_bc/
    CASE_0d_tuned_bc.json
    CASE_0d_tuned_bc.csv
    tuning_history.csv
    tuning_summary.md
```

The tuned JSON is the source of truth for the updated terminal BC values.

# Step 4: Create and Run the Tuned 3D CFD Case

[`apply_tuned_bc_to_xml.py`](apply_tuned_bc_to_xml.py) transfers terminal RCR
or resistance values from the Step 3 JSON into a new 3D `solver.xml`. It only
changes terminal BC parameters; mesh paths, prescribed inflows, solver settings,
and all other XML content are retained from the initial CFD input.

```bash
python3 pipeline_script/apply_tuned_bc_to_xml.py \
  CASE_DIR/01-first-pass/solver.xml \
  CASE_DIR/tuned_bc/CASE_0d_tuned_bc.json \
  CASE_DIR/05-rerun/solver.xml
```

Create `CASE_DIR/05-rerun/`, copy or link the mesh and prescribed-flow files
needed by the XML, submit the new `solver.xml` to the 3D solver, then
postprocess the rerun:

```bash
python3 pipeline_script/postprocess_cfd.py CASE_DIR \
  --results CASE_DIR/05-rerun/48-procs \
  --output CASE_DIR/05-rerun
```

Compare the final-cycle mean Fontan pressure and the individual LPA/RPA flows
with the clinical targets. The default acceptance criterion is 5% relative error
for mean Fontan pressure and each pulmonary branch flow.
