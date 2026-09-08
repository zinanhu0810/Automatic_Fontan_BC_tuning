#!/usr/bin/env python3
"""Run SimVascular's ROM generator from its embedded Python interpreter.

This file intentionally stays Python 3.5 compatible because older SimVascular
application bundles include Python 3.5.
"""

from __future__ import print_function

import json
import os
import sys


def main():
    if len(sys.argv) != 2:
        raise RuntimeError("Usage: sv_rom_generate.py ROM_RUNTIME_CONFIG.json")

    config_path = os.path.abspath(sys.argv[1])
    with open(config_path, "r") as stream:
        config = json.load(stream)

    from sv_rom_simulation import generate_1d_mesh

    result = generate_1d_mesh.run(**config["rom_arguments"])
    output_file = config["expected_output_file"]
    if not os.path.isfile(output_file):
        raise RuntimeError(
            "SimVascular ROM generation did not create: {0}. Result: {1}".format(
                output_file, result
            )
        )

    print(result)
    print("Generated raw ROM: {0}".format(output_file))
    return 0


if __name__ == "__main__":
    sys.exit(main())
