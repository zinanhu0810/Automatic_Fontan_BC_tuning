#!/usr/bin/env python3
"""Rescale P1 source flow waveforms to the 0.64516 s cardiac period."""

from __future__ import annotations

from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "01_initial_3d" / "unsteady"
OUTPUT = HERE / "inflows"
CARDIAC_PERIOD_S = 0.64516
INFLOWS = (
    "cap_Fontan",
    "cap_HV_L",
    "cap_HV_ML",
    "cap_HV_MR",
    "cap_HV_R",
    "cap_SVC_L",
    "cap_SVC_R",
)


def samples(path: Path) -> list[tuple[float, float]]:
    values = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        try:
            time, flow = float(fields[0]), float(fields[1])
        except ValueError:
            continue
        if not values and time > 2.0 and flow > 2.0:
            continue
        values.append((time, flow))
    if len(values) < 2 or values[0][0] != 0.0 or values[-1][0] <= 0.0:
        raise ValueError(f"Unexpected waveform format: {path}")
    return values


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name in INFLOWS:
        waveform = samples(SOURCE / f"{name}.flow")
        source_period = waveform[-1][0]
        path = OUTPUT / f"{name}_physical.flow"
        with path.open("w") as stream:
            stream.write(f"{len(waveform)} 16\n")
            for time, flow in waveform:
                stream.write(f"{time * CARDIAC_PERIOD_S / source_period:.10f}\t{flow:.10f}\n")
        print(path)


if __name__ == "__main__":
    main()
