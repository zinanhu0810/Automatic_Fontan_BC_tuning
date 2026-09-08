# Automatic Fontan Tuning

Reproducible examples for translating and calibrating Fontan lumped-parameter
networks with SimVascular tools.

## Workflows

### Closed-loop genBC to svZeroD

`workflows/closed_loop_0d` contains the legacy Fortran genBC source, the
one-to-one svZeroD model builder, and a notebook comparing three formulations:

- original genBC equations implemented in Python
- a Python hybrid using the svZeroD chamber and valve laws
- the full `svzerodsolver` JSON model

See [`workflows/closed_loop_0d/README.md`](workflows/closed_loop_0d/README.md).

### Patient-specific 3D to 0D

`workflows/three_d_to_0d` contains a command-line wrapper around
SimVascular's `sv_rom_simulation` package. The P15 example includes the surface
model, centerlines, generated 0D network, a resistance-only steady case, steady
3D validation summaries, branch-resistance calibration, and pressure projected
onto the centerline.

See [`workflows/three_d_to_0d/README.md`](workflows/three_d_to_0d/README.md).

## Software

Install the Python packages in `requirements.txt`. The workflows additionally
require SimVascular's `sv_rom_simulation` runtime and the `svzerodsolver` and
`svzerodcalibrator` executables on `PATH` for the stages that invoke them.

Raw 3D processor output, restart files, compiled binaries, solver logs, and
large reproducible time-series files are intentionally excluded. Compact
reference metrics and visualization files are retained.

## Public-release check

The repository is organized for publication, but the owner must confirm that
the coded P15 geometry and associated boundary/summary data are approved for
public release. See [`docs/PUBLICATION_CHECKLIST.md`](docs/PUBLICATION_CHECKLIST.md)
before pushing to a public GitHub repository.

