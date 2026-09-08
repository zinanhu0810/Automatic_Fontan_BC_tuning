#!/usr/bin/env python3
"""Tune the full 0D Fontan sv0D model against clinical targets.

The tuner uses positive log-multipliers around a baseline JSON. This keeps
chamber elastances and valve/vessel resistances positive while allowing
bounded optimization with scipy.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any


P_CONV = 1333.22
ML_PER_S_TO_L_PER_MIN = 0.06

REQUIRED_SERIES = {
    "Vc:left_ventricle",
    "flow:aorta_arch:J3",
    "flow:collateral_left:J21",
    "flow:collateral_right:J23",
    "flow:ivc:J6",
    "flow:svc:J6",
    "flow:J7:left_pul_art",
    "flow:J7:right_pul_art",
    "pressure:Fontan:J7",
    "pressure:J2:left_atrium",
    "pressure:left_ventricle:aortic",
    "pressure:aorta_arch:J3",
}


@dataclass
class ActiveParameter:
    config: dict[str, Any]
    baseline_value: float
    log_bounds: tuple[float, float]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def resolve_path(base_dir: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base_dir / path


def find_named(items: list[dict[str, Any]], key_options: tuple[str, ...], name: str) -> dict[str, Any]:
    for item in items:
        for key in key_options:
            if item.get(key) == name:
                return item
    raise KeyError(f"Could not find {name!r} with keys {key_options}")


def get_parameter_value(model: dict[str, Any], param: dict[str, Any]) -> float:
    kind = param["kind"]
    if kind == "chamber":
        chamber = find_named(model["chambers"], ("name",), param["target"])
        return float(chamber["values"][param["field"]])
    if kind == "valve":
        valve = find_named(model["valves"], ("name",), param["target"])
        return float(valve["params"][param["field"]])
    if kind == "vessel":
        vessel = find_named(model["vessels"], ("vessel_name", "name"), param["target"])
        return float(vessel["zero_d_element_values"][param["field"]])
    raise ValueError(f"Unsupported parameter kind: {kind}")


def set_parameter_value(model: dict[str, Any], param: dict[str, Any], value: float) -> None:
    kind = param["kind"]
    if kind == "chamber":
        chamber = find_named(model["chambers"], ("name",), param["target"])
        chamber["values"][param["field"]] = float(value)
        return
    if kind == "valve":
        valve = find_named(model["valves"], ("name",), param["target"])
        valve["params"][param["field"]] = float(value)
        return
    if kind == "vessel":
        vessel = find_named(model["vessels"], ("vessel_name", "name"), param["target"])
        vessel["zero_d_element_values"][param["field"]] = float(value)
        return
    raise ValueError(f"Unsupported parameter kind: {kind}")


def active_parameters(model: dict[str, Any], config: dict[str, Any], groups: set[str] | None) -> list[ActiveParameter]:
    active: list[ActiveParameter] = []
    for param in config["parameters"]:
        if groups and param.get("group") not in groups:
            continue
        if not groups and not param.get("active", False):
            continue
        baseline_value = get_parameter_value(model, param)
        lo, hi = param["bounds_multiplier"]
        if baseline_value <= 0.0:
            raise ValueError(f"{param['name']} has non-positive baseline value {baseline_value}")
        active.append(
            ActiveParameter(
                config=param,
                baseline_value=baseline_value,
                log_bounds=(math.log(float(lo)), math.log(float(hi))),
            )
        )
    return active


def apply_log_multipliers(
    baseline_model: dict[str, Any],
    active: list[ActiveParameter],
    x_log: list[float],
    sim_overrides: dict[str, Any],
) -> dict[str, Any]:
    model = json.loads(json.dumps(baseline_model))
    model.setdefault("simulation_parameters", {}).update(sim_overrides)
    for spec, log_multiplier in zip(active, x_log):
        value = spec.baseline_value * math.exp(float(log_multiplier))
        set_parameter_value(model, spec.config, value)
    return model


def run_solver(solver: Path, json_path: Path, result_path: Path, timeout_s: int) -> None:
    completed = subprocess.run(
        [str(solver), str(json_path), str(result_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout_s,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"svzerodsolver failed with code {completed.returncode}\n"
            f"Command: {solver} {json_path} {result_path}\n"
            f"Output:\n{completed.stdout[-4000:]}"
        )


def csv_max_time(result_csv: Path) -> float:
    max_time = -math.inf
    with result_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            max_time = max(max_time, float(row["time"]))
    if not math.isfinite(max_time):
        raise ValueError(f"No time values found in {result_csv}")
    return max_time


def collect_last_cycle_series(result_csv: Path, cardiac_period: float) -> dict[str, list[float]]:
    max_time = csv_max_time(result_csv)
    t_start = max_time - cardiac_period
    series: dict[str, list[float]] = {name: [] for name in REQUIRED_SERIES}
    with result_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            if name not in series:
                continue
            if float(row["time"]) >= t_start:
                series[name].append(float(row["y"]))
    missing = [name for name, values in series.items() if not values]
    if missing:
        raise ValueError(f"Missing last-cycle series in {result_csv}: {missing}")
    return series


def compute_metrics_from_csv(result_csv: Path, cardiac_period: float) -> dict[str, float]:
    s = collect_last_cycle_series(result_csv, cardiac_period)

    def mean(name: str) -> float:
        return float(sum(s[name]) / len(s[name]))

    def max_value(name: str) -> float:
        return float(max(s[name]))

    def min_value(name: str) -> float:
        return float(min(s[name]))

    lpa = mean("flow:J7:left_pul_art")
    rpa = mean("flow:J7:right_pul_art")
    ao = mean("flow:aorta_arch:J3") * ML_PER_S_TO_L_PER_MIN
    collateral = (
        mean("flow:collateral_left:J21") + mean("flow:collateral_right:J23")
    ) * ML_PER_S_TO_L_PER_MIN

    return {
        "SVSV_mL": max_value("Vc:left_ventricle") - min_value("Vc:left_ventricle"),
        "AO_L_min": ao,
        "collateral_L_min": collateral,
        "system_L_min": ao - collateral,
        "IVC_L_min": mean("flow:ivc:J6") * ML_PER_S_TO_L_PER_MIN,
        "SVC_L_min": mean("flow:svc:J6") * ML_PER_S_TO_L_PER_MIN,
        "pulmonary_L_fraction": lpa / (lpa + rpa),
        "Fontan_mean_mmHg": mean("pressure:Fontan:J7") / P_CONV,
        "LA_mean_mmHg": mean("pressure:J2:left_atrium") / P_CONV,
        "LV_peak_mmHg": max_value("pressure:left_ventricle:aortic") / P_CONV,
        "Aorta_peak_mmHg": max_value("pressure:aorta_arch:J3") / P_CONV,
        "Aorta_min_mmHg": min_value("pressure:aorta_arch:J3") / P_CONV,
    }


def objective_value(metrics: dict[str, float], targets: dict[str, Any]) -> float:
    total = 0.0
    for name, target_spec in targets.items():
        target = float(target_spec["target"])
        weight = float(target_spec.get("weight", 1.0))
        value = float(metrics[name])
        total += weight * ((value - target) / (abs(target) + 1.0e-12)) ** 2
    return float(total)


def write_metrics_csv(path: Path, metrics: dict[str, float], targets: dict[str, Any]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "simulated", "target", "relative_error", "weight"])
        for name, target_spec in targets.items():
            target = float(target_spec["target"])
            value = float(metrics[name])
            rel_error = abs(value - target) / (abs(target) + 1.0e-12)
            writer.writerow([name, value, target, rel_error, target_spec.get("weight", 1.0)])


def write_history_header(path: Path, active: list[ActiveParameter], targets: dict[str, Any]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["iteration", "objective"]
            + [f"multiplier:{p.config['name']}" for p in active]
            + [f"metric:{name}" for name in targets]
        )


def append_history(
    path: Path,
    iteration: int,
    objective: float,
    x_log: list[float],
    metrics: dict[str, float],
    targets: dict[str, Any],
) -> None:
    with path.open("a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [iteration, objective]
            + [math.exp(float(v)) for v in x_log]
            + [metrics[name] for name in targets]
        )


def evaluate_candidate(
    x_log: list[float],
    *,
    baseline_model: dict[str, Any],
    active: list[ActiveParameter],
    config: dict[str, Any],
    targets: dict[str, Any],
    run_dir: Path,
    iteration: int,
    timeout_s: int,
) -> tuple[float, dict[str, float], Path]:
    model = apply_log_multipliers(
        baseline_model,
        active,
        x_log,
        config.get("simulation_overrides", {}),
    )
    json_path = run_dir / f"candidate_{iteration:04d}.json"
    result_path = run_dir / f"candidate_{iteration:04d}.csv"
    write_json(json_path, model)
    run_solver(Path(config["solver"]), json_path, result_path, timeout_s)
    cardiac_period = float(model["simulation_parameters"]["cardiac_period"])
    metrics = compute_metrics_from_csv(result_path, cardiac_period)
    return objective_value(metrics, targets), metrics, json_path


def clip(value: float, lo: float, hi: float) -> float:
    return min(max(value, lo), hi)


def pattern_search(
    objective,
    x0: list[float],
    bounds: list[tuple[float, float]],
    maxiter: int,
    initial_step: float,
    min_step: float,
    shrink: float,
) -> dict[str, Any]:
    x = [clip(v, lo, hi) for v, (lo, hi) in zip(x0, bounds)]
    best_value = objective(x)
    step = [float(initial_step)] * len(x)
    nfev = 1
    success = False
    message = "maximum iterations reached"

    for sweep in range(1, maxiter + 1):
        improved = False
        for i in range(len(x)):
            for direction in (-1.0, 1.0):
                trial = list(x)
                lo, hi = bounds[i]
                trial[i] = clip(trial[i] + direction * step[i], lo, hi)
                if trial[i] == x[i]:
                    continue
                value = objective(trial)
                nfev += 1
                if value < best_value:
                    x = trial
                    best_value = value
                    improved = True
                    break
        if not improved:
            step = [v * shrink for v in step]
            if max(step) < min_step:
                success = True
                message = "step size below min_step"
                break
        print(
            f"sweep={sweep} best={best_value:.6g} "
            f"max_step={max(step):.4g}"
        )

    return {
        "x": x,
        "fun": best_value,
        "nfev": nfev,
        "success": success,
        "message": message,
    }


def run_tuning(args: argparse.Namespace) -> None:
    base_dir = Path(__file__).resolve().parent
    config = load_json(resolve_path(base_dir, args.config))
    targets_config = load_json(resolve_path(base_dir, config["targets_json"]))
    targets = targets_config["targets"]
    baseline_path = resolve_path(base_dir, config["baseline_json"])
    baseline_model = load_json(baseline_path)
    groups = set(args.groups.split(",")) if args.groups else None
    active = active_parameters(baseline_model, config, groups)
    if not active:
        raise ValueError("No active parameters selected.")

    run_dir = resolve_path(base_dir, config.get("work_dir", "runs"))
    run_dir.mkdir(parents=True, exist_ok=True)
    history_path = run_dir / "tuning_history.csv"
    write_history_header(history_path, active, targets)

    best = {"objective": math.inf, "metrics": None, "json_path": None, "x_log": None}
    counter = {"iteration": 0}

    def objective(x_log: list[float]) -> float:
        counter["iteration"] += 1
        try:
            value, metrics, json_path = evaluate_candidate(
                x_log,
                baseline_model=baseline_model,
                active=active,
                config=config,
                targets=targets,
                run_dir=run_dir,
                iteration=counter["iteration"],
                timeout_s=args.timeout,
            )
        except Exception as exc:
            print(f"[{counter['iteration']:04d}] failed: {exc}")
            return args.failure_penalty

        append_history(history_path, counter["iteration"], value, x_log, metrics, targets)
        if value < best["objective"]:
            best.update(objective=value, metrics=metrics, json_path=json_path, x_log=list(x_log))
            shutil.copyfile(json_path, run_dir / "best_tuned.json")
            write_metrics_csv(run_dir / "best_metrics.csv", metrics, targets)
        print(f"[{counter['iteration']:04d}] objective={value:.6g}")
        return value

    x0 = [0.0] * len(active)
    bounds = [spec.log_bounds for spec in active]

    if args.evaluate_only:
        value = objective(x0)
        print(f"Baseline objective: {value:.6g}")
        return

    optimizer = config.get("optimizer", {})
    maxiter = args.maxiter if args.maxiter is not None else int(optimizer.get("maxiter", 80))
    result = pattern_search(
        objective,
        x0,
        bounds=bounds,
        maxiter=maxiter,
        initial_step=float(optimizer.get("initial_step", 0.35)),
        min_step=float(optimizer.get("min_step", 0.025)),
        shrink=float(optimizer.get("shrink", 0.5)),
    )

    summary_path = run_dir / "optimizer_result.json"
    summary = {
        "success": bool(result["success"]),
        "message": str(result["message"]),
        "objective": float(result["fun"]),
        "nfev": int(result["nfev"]),
        "active_parameters": [p.config["name"] for p in active],
        "multipliers": {p.config["name"]: math.exp(float(v)) for p, v in zip(active, result["x"])},
        "best_objective_seen": float(best["objective"]),
    }
    write_json(summary_path, summary)
    print(f"Wrote {summary_path}")
    print(f"Wrote {run_dir / 'best_tuned.json'}")
    print(f"Wrote {run_dir / 'best_metrics.csv'}")


def evaluate_existing(args: argparse.Namespace) -> None:
    base_dir = Path(__file__).resolve().parent
    config = load_json(resolve_path(base_dir, args.config))
    targets_config = load_json(resolve_path(base_dir, config["targets_json"]))
    targets = targets_config["targets"]
    baseline_model = load_json(resolve_path(base_dir, config["baseline_json"]))
    cardiac_period = float(baseline_model["simulation_parameters"].get("cardiac_period", 0.666667))
    metrics = compute_metrics_from_csv(resolve_path(base_dir, args.evaluate_existing), cardiac_period)
    out_path = resolve_path(base_dir, args.output)
    write_metrics_csv(out_path, metrics, targets)
    print(f"Objective: {objective_value(metrics, targets):.6g}")
    print(f"Wrote {out_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="tune_config.json", help="Tuner config JSON.")
    parser.add_argument("--groups", default=None, help="Comma-separated active parameter groups to use.")
    parser.add_argument("--maxiter", type=int, default=None, help="Override optimizer maxiter.")
    parser.add_argument("--timeout", type=int, default=600, help="Solver timeout per candidate in seconds.")
    parser.add_argument("--failure-penalty", type=float, default=1.0e6)
    parser.add_argument("--evaluate-only", action="store_true", help="Run only the baseline candidate.")
    parser.add_argument("--evaluate-existing", default=None, help="Evaluate an existing sv0D result CSV.")
    parser.add_argument("--output", default="baseline_metrics.csv", help="Output CSV for --evaluate-existing.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.evaluate_existing:
        evaluate_existing(args)
    else:
        run_tuning(args)


if __name__ == "__main__":
    main()
