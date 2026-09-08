# sv0D Full 0D Fontan Tuner

This folder sets up automatic tuning for the full 0D sv0D Fontan LPN.

The baseline is copied from:

`../sv0D_LPN/Fontan_genBC_1to1_sv0D.json`

The first tuning stage focuses on the heart layer:

- `left_atrium` `LinearElastanceChamber`: `Emax`, `Epass`
- `left_ventricle` `LinearElastanceChamber`: `Emax`, `Epass`
- `mitral` `PiecewiseValve`: `Rmin`, `Rmax`
- `aortic` `PiecewiseValve`: `Rmin`, `Rmax`

Optional vascular parameters are also listed in `tune_config.json`, but are inactive by default.

## Targets

Targets are stored in `clinical_targets.json`.

The current target vector is derived from the final cardiac cycle of the full 0D genBC baseline in `../sv0D_LPN/comparison_timeseries.csv`:

| Metric | Target |
| --- | ---: |
| SV stroke volume | 37.842 mL |
| Aortic output | 3.395 L/min |
| Collateral flow | 0.716 L/min |
| Systemic flow | 2.679 L/min |
| IVC flow | 1.691 L/min |
| SVC flow | 0.995 L/min |
| Pulmonary split | L 27.654%, R 72.346% |
| Fontan pressure | 10.278 mmHg |
| LA mean pressure | 7.261 mmHg |
| LV peak pressure | 112.545 mmHg |
| Aortic pressure | 94.440/62.411 mmHg |

## Quick Commands

Evaluate the existing full 0D result without rerunning svZeroD:

```bash
python3 tune_full0d.py --evaluate-existing ../sv0D_LPN/Fontan_genBC_1to1_sv0D_result.csv --output baseline_existing_metrics.csv
```

Run one solver evaluation using the baseline JSON and the tuning simulation overrides:

```bash
python3 tune_full0d.py --evaluate-only --groups chamber,valve
```

Run a short chamber/valve tuning test:

```bash
python3 tune_full0d.py --groups chamber,valve --maxiter 5
```

Run a longer chamber/valve tuning pass:

```bash
python3 tune_full0d.py --groups chamber,valve --maxiter 80
```

If chamber/valve tuning cannot satisfy flow split or venous-flow targets, include optional vascular knobs:

```bash
python3 tune_full0d.py --groups chamber,valve,vascular_optional --maxiter 100
```

## Outputs

Tuning writes into `runs_genbc_targets/`:

- `tuning_history.csv`: every solver candidate and metric value
- `best_tuned.json`: best JSON seen so far
- `best_metrics.csv`: target comparison for `best_tuned.json`
- `optimizer_result.json`: final optimizer summary
- `candidate_####.json` and `candidate_####.csv`: individual runs

The compact `reference_results/` directory contains the best model, metric
comparison, optimizer summary, and tuning history from the documented run.
Individual candidate files are intentionally not included.

## Notes

The local `pysvzerod` module exposes a `calibrate` binding, but there is no
separate `svzerod_tuner` command. This setup therefore uses a transparent
bounded pattern-search optimizer around `svzerodsolver`. The solver executable
must be available on `PATH`.
