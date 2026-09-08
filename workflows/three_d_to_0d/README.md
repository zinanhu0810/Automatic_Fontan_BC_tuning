# Patient-specific 3D-to-0D ROM pipeline

This directory contains a command-line wrapper around SimVascular's
`sv_rom_simulation` package. It extracts centerlines and generates a geometric
0D `BloodVessel` network without using the SimVascular GUI.

The upstream implementation is maintained at:

<https://github.com/SimVascular/SimVascular/tree/master/Python/site-packages/sv_rom_simulation>

## Run P15

From this directory:

```bash
python3 generate_patient_0d.py P15
```

Add `--run-baseline` to run the generated model with the locally installed
`svzerodsolver` and verify finite output plus mean external-flow conservation.

The command reads `P15/solver.xml`, identifies prescribed-flow and RCR faces,
runs SimVascular's centerline/ROM code, and writes:

- `P15/rom/centerlines.vtp`: extracted and grouped centerline tree
- `P15/rom/P15_0d_raw.json`: direct single-inlet SimVascular ROM output
- `P15/rom/P15_0d.json`: multi-inlet Fontan svZeroD model
- `P15/rom/P15_0d_baseline.csv`: optional baseline solver output
- `P15/rom/generation_manifest.json`: source and boundary-condition provenance
- `P15/rom/generation_report.md`: compact validation summary

The root cap only defines the orientation of the centerline tree. P15 uses
`cap_fontan_x`. Other physical inflows are converted from placeholder terminal
RCRs to signed `FLOW` boundary conditions after ROM generation. All pulmonary
RCR values are copied directly from the fluid equation in `solver.xml`.

## Add another patient

Place the patient's closed exterior surface, individual cap/wall surfaces,
`solver.xml`, and waveform files in a patient directory. Add a
`patient_0d_config.json` using P15 as the template, then run the same command.

Use `--prepare-only` to inspect generated ROM inputs before centerline
extraction. Use `--reuse-centerlines` after a successful first run to avoid
recomputing centerlines while changing ROM parameters.

## Included validation and calibration

- `P15/steady` contains the steady 3D setup, postprocessor, and compact metrics
  from the final 100 samples. Raw processor output is not tracked.
- `P15/steady_0d` converts the geometric ROM to constant inflow plus resistance
  outlets and projects the resulting pressure field onto the centerline.
- `P15/step3_calibration` contains the `svZeroDCalibrator` input, calibration
  driver, calibrated model, compact metrics, and centerline visualization.

Re-running the 3D postprocessing or calibration requires raw svMultiPhysics
output under `P15/steady/48-procs`. Those large files are intentionally outside
the public repository; the derived summaries and calibrator observations are
included for audit.

## Current model boundary

This step creates the untuned geometry-derived 0D network. It does not yet fit
the vessel coefficients to the 3D pressure/flow waveforms or tune downstream
boundary conditions. Those are separate calibration stages so geometry,
3D-to-0D fitting, and physiological BC tuning remain independently auditable.
