# Fontan Closed-Loop 0D/3D Coupled Example

This is a compact, working example of a closed-loop Fontan circulation coupled to a rigid-wall 3D CFD domain through svZeroD and svMultiPhysics.

## Contents

- `0d/Fontan_full_0d.json`: standalone full closed-loop Fontan LPN.
- `0d/Fontan_3d_coupled.json`: coupled 0D model with the Fontan/central PA region represented by the 3D domain.
- `mesh/`: volume mesh and coupled face surfaces.
- `figures/fontan_full_0d_lpn_schematic.png`: full closed-loop LPN schematic.
- `solver.template.xml`: coupled svMultiPhysics input template.
- `run_sherlock.sh`: Slurm launcher.
- `scripts/postprocess_coupled.py`: final-cycle pressure and flow summary.
- `reference_results/`: validated summary metrics; transient VTK output is excluded.

## Coupling Faces

The JSON `external_solver_coupling_blocks` and the XML face names must match exactly: `lpa`, `rpa`, `rpa2`, `rpa3`, `rpa4`, `rpa5`, `svc`, and `ivc`.

## Run on Sherlock

```bash
export SV_MULTIPHYSICS_BIN=/path/to/svmultiphysics
export SVZERO_INTERFACE_LIB=/path/to/libsvzero_interface.so
sbatch run_sherlock.sh
```

The launcher writes `solver.xml` from the template using the local svZeroD interface-library path. The solver itself and the shared coupling library are not distributed in this example.

## Postprocess

```bash
python3 scripts/postprocess_coupled.py 24-procs --period 0.666667
```

The postprocessor requires NumPy and converts pressure from dyn/cm2 to mmHg.
