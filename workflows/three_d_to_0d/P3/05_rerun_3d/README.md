# P3 Step 5: Rerun 3D With Tuned RCR BCs

This case reruns the P3 unsteady 3D CFD model using the RCR resistances tuned
in Step 4. It retains the initial-case mesh, prescribed inflows, wall boundary
condition, time step (`0.0009091 s`), six-cycle duration, RCR capacitance,
distal pressure, and initialized capacitor pressure.

`build_tuned_unsteady_case.py` is the only source of the Step-5 XML changes.
It reads `04_tune_bc/unsteady/P3_0d_tuned_bc.json`, checks that all 18
pulmonary RCR outlets correspond one-to-one with the initial 3D XML, and
updates only `Proximal_resistance` and `Distal_resistance`.

Files to run:

```bash
python3 build_tuned_unsteady_case.py
sbatch run_sherlock.sh
```

The values written to `solver.xml` are listed in `tuned_rcr_summary.csv` in
the cgs units expected by svMultiPhysics.
