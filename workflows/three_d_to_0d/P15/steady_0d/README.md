# P15 steady 0D case

This case is the standalone 0D analogue of the P15 steady 3D simulation.

- The patient-specific geometry-derived `BloodVessel` network is copied without
  changing its vessel RCL parameters.
- Each prescribed FLOW waveform is replaced by its trapezoidal cycle mean.
- Each outlet RCR boundary is replaced by a RESISTANCE boundary with
  `R = Rp + Rd` and the same distal pressure.
- The flow signs follow the centerline orientation used by `P15_0d.json`.

Build and run from `workflows/three_d_to_0d`:

```bash
python3 P15/steady_0d/build_steady_0d.py
```

The command writes the converted JSON, the svZeroD result CSV, a boundary-level
summary CSV, and mass-conservation metrics in `simulation_summary.json`.

Project the solved pressure and flow fields onto the ROM centerline:

```bash
python3 P15/steady_0d/project_0d_to_centerline.py
```

This writes a ParaView-ready VTP with point-data arrays, a node-level CSV, and
a two-view PNG rendering. Pressure is linearly interpolated through each 0D
vessel and harmonically interpolated across centerline junction regions.
