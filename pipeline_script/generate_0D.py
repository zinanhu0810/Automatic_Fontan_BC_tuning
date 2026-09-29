#!/usr/bin/env python3
"""Generate a boundary-faithful svZeroD model from a Step 0 CFD case.

This driver derives ROM input from step0_config.json and solver.xml, then calls
the shared SimVascular-ROM implementation in sv_rom_0d_toolkit.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
TOOLKIT = HERE.parents[1] / "sv_rom_0d_toolkit" / "generate_0d.py"
ROM_DEFAULTS = {
    "material_model": "OLUFSEN", "element_size": 0.01,
    "min_num_elements": 5, "seg_min_num": 1, "seg_size": 9999,
    "olufsen_material_k1": 0.0, "olufsen_material_k2": -22.5267,
    "olufsen_material_k3": 10000000.0, "olufsen_material_exp": 1.0,
    "olufsen_material_pressure": 0.0,
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--step0-config", type=Path, default=Path("step0_config.json"))
    parser.add_argument("--root-surface", help="Override metrics.fontan_pressure_face as the root")
    parser.add_argument("--output-directory", default="0D")
    parser.add_argument("--units", choices=("cm", "mm"), default="cm")
    parser.add_argument("--element-size", type=float, default=0.01)
    parser.add_argument("--simvascular", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--reuse-centerlines", action="store_true")
    return parser.parse_args()


def build_rom_config(case: Path, step0: dict, args: argparse.Namespace) -> dict:
    for key in ("case_name", "mesh", "inlets", "metrics"):
        if key not in step0:
            raise ValueError(f"step0_config.json is missing '{key}'.")
    mesh = step0["mesh"]
    mesh_dir = mesh.get("directory", "mesh")
    surface_dir = mesh.get("surface_directory", "mesh-surfaces")
    root = args.root_surface or step0["metrics"].get("fontan_pressure_face")
    inlets = {item["name"] for item in step0["inlets"]}
    if root not in inlets:
        raise ValueError(f"Root '{root}' must be a prescribed-flow inlet: {sorted(inlets)}")
    required = {
        "solver.xml": case / "solver.xml",
        "mesh-complete.exterior.vtp": case / mesh_dir / "mesh-complete.exterior.vtp",
        "mesh-surfaces": case / mesh_dir / surface_dir,
    }
    for label, path in required.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing required {label}: {path}")
    options = dict(ROM_DEFAULTS)
    options["element_size"] = args.element_size
    options.update(step0.get("rom_options", {}))
    return {
        "case_name": step0["case_name"], "solver_xml": "solver.xml",
        "surface_model": str(required["mesh-complete.exterior.vtp"].relative_to(case)),
        "boundary_surfaces_directory": str(required["mesh-surfaces"].relative_to(case)),
        "root_surface": root, "output_directory": args.output_directory,
        "units": args.units, "wall_surface_prefixes": ["wall"],
        "rom_options": options,
    }


def main() -> None:
    args = arguments()
    case = args.case_directory.expanduser().resolve()
    config_path = args.step0_config if args.step0_config.is_absolute() else case / args.step0_config
    if not TOOLKIT.is_file():
        raise FileNotFoundError(f"Shared ROM generator not found: {TOOLKIT}")
    with config_path.open() as stream:
        step0 = json.load(stream)
    rom_config = build_rom_config(case, step0, args)
    with tempfile.NamedTemporaryFile("w", suffix=".json", prefix="rom_0d_", delete=False) as stream:
        json.dump(rom_config, stream, indent=2); stream.write("\n")
        temporary_config = Path(stream.name)
    command = [sys.executable, str(TOOLKIT), str(case), "--config", str(temporary_config)]
    if args.simvascular: command.extend(("--simvascular", str(args.simvascular.expanduser().resolve())))
    if args.prepare_only: command.append("--prepare-only")
    if args.reuse_centerlines: command.append("--reuse-centerlines")
    try:
        subprocess.run(command, check=True)
        output = Path(args.output_directory)
        output = output if output.is_absolute() else case / output
        (output / "generation_manifest.json").unlink(missing_ok=True)
    finally:
        temporary_config.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
