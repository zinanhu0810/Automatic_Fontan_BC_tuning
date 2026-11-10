#!/usr/bin/env python3
"""Summarize the final cardiac cycle of a coupled svMultiPhysics run."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

P_CONV = 1333.22
RPA_FACES = ["rpa", "rpa2", "rpa3", "rpa4", "rpa5"]


def read_boundary_integral(path: Path) -> dict[str, np.ndarray]:
    lines = path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.strip().startswith("step"))
    names = lines[header].split()
    data = np.loadtxt(lines[header + 1 :])
    return {name: data[:, i] for i, name in enumerate(names)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--period", type=float, default=0.666667)
    args = parser.parse_args()
    pressure = read_boundary_integral(args.results_dir / "B_NS_Pressure_average.txt")
    flux = read_boundary_integral(args.results_dir / "B_NS_Velocity_flux.txt")
    start = min(pressure["time"].max(), flux["time"].max()) - args.period
    p_mask, q_mask = pressure["time"] >= start, flux["time"] >= start
    lpa_q = flux["lpa"][q_mask]
    rpa_q = np.sum(np.vstack([flux[name][q_mask] for name in RPA_FACES]), axis=0)
    print("Final-cycle coupled results")
    print(f"Mean Fontan pressure: {(np.mean(pressure['svc'][p_mask]) + np.mean(pressure['ivc'][p_mask])) / (2 * P_CONV):.3f} mmHg")
    print(f"Mean LPA flow: {np.mean(lpa_q):.3f} mL/s")
    print(f"Mean RPA flow: {np.mean(rpa_q):.3f} mL/s")
    print(f"LPA flow split: {100 * np.mean(lpa_q / (lpa_q + rpa_q)):.2f}%")


if __name__ == "__main__":
    main()
