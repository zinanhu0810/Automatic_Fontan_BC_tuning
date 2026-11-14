# P3 Step 2: Generate Geometric 0D Model

This stage uses the shared `sv_rom_0d_toolkit` and the P3 Step 1 `solver.xml`.
The selected root (`cap_Fontan`) sets only the centerline orientation. Every
prescribed FLOW waveform and RCR boundary condition is restored exactly in the
resulting svZeroD JSON.

```bash
python3 ../../sv_rom_0d_toolkit/generate_0d.py . --config rom_0d_config.json \
  --simvascular /Applications/SimVascular.app/Contents/Resources/simvascular \
  --run-solver --svzerodsolver /Users/zinanhu/anaconda3/bin/svzerodsolver
```
