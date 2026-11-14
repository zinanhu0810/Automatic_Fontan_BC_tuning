# P18 Step 5: Rerun 3D with Tuned RCR BCs

This self-contained directory is the P18 unsteady 3D rerun case. It uses the
same mesh, six prescribed physical inflows, time step, number of cycles, and
solver controls as `01_initial_3d/unsteady`.

`build_tuned_unsteady_case.py` reads the Step 4 tuned 0D model and writes
`solver.xml`. Only four XML values are changed from the initial 3D case:

| Face | Proximal resistance (`Rp`) | Distal resistance (`Rd`) |
| --- | ---: | ---: |
| `cap_pa_x` (LPA) | 37.8034915080021 | 378.034915080021 |
| `cap_pa_x_2` (RPA) | 42.8903904362189 | 428.903904362189 |

RCR capacitance remains `1.5e-4 s`; distal pressure and initialized capacitor
pressure are unchanged. The tuned values are expressed in the same cgs units
as the original `solver.xml`.

Regenerate the XML after rerunning Step 4:

```bash
python3 build_tuned_unsteady_case.py
sbatch run_sherlock.sh
```

After the solver finishes, place the solver output in `results/48-procs/` for
the final clinical-target postprocessing step.
