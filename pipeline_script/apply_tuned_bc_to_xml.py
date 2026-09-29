#!/usr/bin/env python3
"""Apply terminal RESISTANCE or RCR values from an svZeroD JSON model to solver.xml."""

from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver_xml", type=Path, help="Original svMultiPhysics solver.xml")
    parser.add_argument("model_json", type=Path, help="Tuned svZeroD JSON model")
    parser.add_argument("output_xml", type=Path, help="Updated solver.xml")
    return parser.parse_args()


def main() -> None:
    args = arguments()
    model = json.loads(args.model_json.read_text())
    tuned = {
        boundary["bc_name"]: boundary
        for boundary in model["boundary_conditions"]
        if boundary["bc_type"] in {"RESISTANCE", "RCR"}
    }
    tree = ET.parse(args.solver_xml)
    fluid = next(
        equation for equation in tree.getroot().findall("Add_equation")
        if equation.attrib.get("type", "").lower() == "fluid"
    )
    applied: set[str] = set()
    for bc in fluid.findall("Add_BC"):
        name = bc.attrib.get("name", "")
        dependence = (bc.findtext("Time_dependence") or "").lower()
        key = f"RCR_{name}" if dependence == "rcr" else f"RESISTANCE_{name}" if dependence == "resistance" else None
        if key not in tuned:
            continue
        source = tuned[key]
        if source["bc_type"] == "RESISTANCE":
            if dependence != "resistance":
                raise ValueError(f"{name} is not a Resistance BC in solver.xml.")
            bc.find("Value").text = f"{source['bc_values']['R']:.12g}"
            distal = bc.find("Distal_pressure")
            if distal is not None:
                distal.text = f"{source['bc_values'].get('Pd', 0.0):.12g}"
        else:
            if dependence != "rcr":
                raise ValueError(f"{name} is not an RCR BC in solver.xml.")
            values = bc.find("RCR_values")
            for field, key_name in (("Proximal_resistance", "Rp"), ("Capacitance", "C"), ("Distal_resistance", "Rd"), ("Distal_pressure", "Pd")):
                values.find(field).text = f"{source['bc_values'][key_name]:.12g}"
        applied.add(key)
    missing = set(tuned) - applied
    if missing:
        raise ValueError(f"Tuned boundary conditions missing from solver.xml: {sorted(missing)}")
    ET.indent(tree, space="  ")
    args.output_xml.parent.mkdir(parents=True, exist_ok=True)
    tree.write(args.output_xml, encoding="utf-8", xml_declaration=True)
    print(f"Wrote {args.output_xml} with {len(applied)} tuned terminal BCs.")


if __name__ == "__main__":
    main()
