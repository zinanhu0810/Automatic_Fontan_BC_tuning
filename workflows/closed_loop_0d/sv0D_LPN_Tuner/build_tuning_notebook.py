#!/usr/bin/env python3
"""Build the sv0D full 0D tuning notebook."""

import json
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK_PATH = ROOT / "run_full0d_tuning.ipynb"


def cell_source(source: str) -> list[str]:
    lines = source.splitlines(keepends=True)
    return lines if lines else [""]


def markdown_cell(source: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": uuid.uuid4().hex[:8],
        "metadata": {},
        "source": cell_source(source),
    }


def code_cell(source: str) -> dict:
    return {
        "cell_type": "code",
        "id": uuid.uuid4().hex[:8],
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": cell_source(source),
    }


SETUP_CODE = r'''
from pathlib import Path
import csv
import json
import subprocess
import sys

from IPython.display import HTML, display

CASE_DIR = Path.cwd()
if CASE_DIR.name != "sv0D_LPN_Tuner":
    CASE_DIR = Path("sv0D_LPN_Tuner")

TARGETS = json.loads((CASE_DIR / "clinical_targets.json").read_text())["targets"]
CONFIG = json.loads((CASE_DIR / "tune_config.json").read_text())

print(f"Tuner folder: {CASE_DIR.resolve()}")
print(f"Baseline JSON: {(CASE_DIR / CONFIG['baseline_json']).resolve()}")
print(f"Solver: {CONFIG['solver']}")
'''.strip()


EVALUATE_EXISTING_CODE = r'''
existing_result = CASE_DIR.parent / "sv0D_LPN" / "Fontan_genBC_1to1_sv0D_result.csv"
baseline_metrics_csv = CASE_DIR / "baseline_existing_metrics.csv"
cmd = [
    sys.executable,
    str(CASE_DIR / "tune_full0d.py"),
    "--evaluate-existing",
    str(existing_result),
    "--output",
    str(baseline_metrics_csv),
]
completed = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
print(completed.stdout)
if completed.returncode != 0:
    raise RuntimeError("baseline metric extraction failed")

def read_csv_rows(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))

def fmt_value(value):
    value = float(value)
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.2f}"
    return f"{value:.3f}"

def fmt_pct(value):
    return f"{100.0 * float(value):.1f}%"

rows = read_csv_rows(baseline_metrics_csv)
html = """
<h3>Baseline sv0D vs clinical targets</h3>
<table>
<thead><tr><th>metric</th><th>simulated</th><th>target</th><th>relative error</th><th>weight</th></tr></thead>
<tbody>
""" + "\n".join(
    f"<tr><td>{r['metric']}</td><td>{fmt_value(r['simulated'])}</td><td>{fmt_value(r['target'])}</td><td>{fmt_pct(r['relative_error'])}</td><td>{r['weight']}</td></tr>"
    for r in rows
) + "</tbody></table>"
display(HTML(html))
'''.strip()


PARAMETER_SPACE_CODE = r'''
active = [p for p in CONFIG["parameters"] if p.get("active")]
optional = [p for p in CONFIG["parameters"] if not p.get("active")]

def parameter_table(params, title):
    rows = "\n".join(
        f"<tr><td>{p['name']}</td><td>{p['group']}</td><td>{p['kind']}</td><td>{p['target']}</td><td>{p['field']}</td><td>{p['bounds_multiplier']}</td></tr>"
        for p in params
    )
    display(HTML(
        f"<h3>{title}</h3><table><thead><tr><th>name</th><th>group</th><th>kind</th><th>target</th><th>field</th><th>bounds multiplier</th></tr></thead><tbody>{rows}</tbody></table>"
    ))

parameter_table(active, "Default active chamber/valve parameters")
parameter_table(optional, "Optional vascular parameters")
'''.strip()


SMOKE_TEST_CODE = r'''
RUN_SOLVER_SMOKE_TEST = False

if RUN_SOLVER_SMOKE_TEST:
    cmd = [sys.executable, str(CASE_DIR / "tune_full0d.py"), "--evaluate-only", "--groups", "chamber,valve"]
    completed = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    print(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError("solver smoke test failed")
else:
    print("Set RUN_SOLVER_SMOKE_TEST = True to run one baseline svZeroD solver evaluation.")
'''.strip()


TUNING_CODE = r'''
RUN_SHORT_TUNING = False

if RUN_SHORT_TUNING:
    cmd = [
        sys.executable,
        str(CASE_DIR / "tune_full0d.py"),
        "--groups",
        "chamber,valve",
        "--maxiter",
        "5",
    ]
    completed = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    print(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError("short tuning run failed")
else:
    print("Set RUN_SHORT_TUNING = True to run a short chamber/valve tuning pass.")
'''.strip()


HISTORY_CODE = r'''
run_dir = CASE_DIR / CONFIG.get("work_dir", "runs")
history_csv = run_dir / "tuning_history.csv"
best_csv = run_dir / "best_metrics.csv"

if history_csv.exists():
    history_rows = read_csv_rows(history_csv)
    display(HTML(f"<p>History rows: {len(history_rows)}</p>"))
    display(HTML("<pre>" + "\n".join(",".join(str(value) for value in row.values()) for row in history_rows[-5:]) + "</pre>"))
else:
    print("No tuning history yet.")

if best_csv.exists():
    rows = read_csv_rows(best_csv)
    html = """
    <h3>Best tuned candidate vs clinical targets</h3>
    <table>
    <thead><tr><th>metric</th><th>simulated</th><th>target</th><th>relative error</th><th>weight</th></tr></thead>
    <tbody>
    """ + "\n".join(
        f"<tr><td>{r['metric']}</td><td>{fmt_value(r['simulated'])}</td><td>{fmt_value(r['target'])}</td><td>{fmt_pct(r['relative_error'])}</td><td>{r['weight']}</td></tr>"
        for r in rows
    ) + "</tbody></table>"
    display(HTML(html))
else:
    print("No best_metrics.csv yet.")
'''.strip()


def build_notebook() -> dict:
    return {
        "cells": [
            markdown_cell(
                "# Full 0D sv0D Tuning\n\n"
                "This notebook wraps `tune_full0d.py` to evaluate and tune the full 0D Fontan sv0D model against clinical targets."
            ),
            code_cell(SETUP_CODE),
            code_cell(EVALUATE_EXISTING_CODE),
            code_cell(PARAMETER_SPACE_CODE),
            code_cell(SMOKE_TEST_CODE),
            code_cell(TUNING_CODE),
            code_cell(HISTORY_CODE),
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main() -> None:
    NOTEBOOK_PATH.write_text(json.dumps(build_notebook(), indent=2) + "\n")
    print(f"Wrote {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()
