# P18 Step 2: Generate Geometric 0D Model

This stage generates a vascular svZeroD model from the P18 Step 1 3D mesh with
the shared `sv_rom_0d_toolkit`. It references the source files in
`../01_initial_3d/unsteady` and does not duplicate the mesh or waveform files.

`cap_AZG` is used only as the rooted centerline orientation. All six
prescribed physical inflows and the two RCR outlet boundary conditions are
restored from the Step 1 `solver.xml` in the final P18 JSON.

## Reproduce

```bash
python3 ../../sv_rom_0d_toolkit/generate_0d.py . \
  --config rom_0d_config.json \
  --simvascular /Applications/SimVascular.app/Contents/Resources/simvascular \
  --run-solver \
  --svzerodsolver /Users/zinanhu/anaconda3/bin/svzerodsolver
```

The generated model and solver result are written under `rom/`.
