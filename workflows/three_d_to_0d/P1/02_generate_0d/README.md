# P1 Step 2: Generate Geometric 0D Model

This stage constructs an svZeroD model from the P1 3D geometry using the shared
`sv_rom_0d_toolkit`. The rooted network uses `cap_Fontan`; this selects only
the centerline orientation. All seven prescribed flow waveforms and all 17 RCR
boundary-condition values are restored exactly from the Step 1 `solver.xml`.

The original P1 source waveforms span one second, while the CFD cycle is
`0.64516 s`. First generate local, rescaled copies for the 0D model:

```bash
python3 prepare_physical_inflows.py
python3 ../../sv_rom_0d_toolkit/generate_0d.py . --config rom_0d_config.json \
  --simvascular /Applications/SimVascular.app/Contents/Resources/simvascular \
  --run-solver --svzerodsolver /Users/zinanhu/anaconda3/bin/svzerodsolver
```

Expected outputs are under `rom/`, including `P1_0d.json`, `P1_0d.csv`,
`centerlines.vtp`, and `generation_report.md`.
