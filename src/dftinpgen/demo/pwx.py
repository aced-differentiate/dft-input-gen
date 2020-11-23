"""Demo generating input files for doing a calculation with pw.x."""

import sys
import json
import pkg_resources
import argparse

from ase import io as ase_io

from dftinpgen.utils import read_crystal_structure
from dftinpgen.qe.pwx import PwxInputGenerator


def get_parser():
    description = """Input file generation for pw.x."""
    parser = argparse.ArgumentParser(description=description)

    # Required:
    crystal_structure = "(REQUIRED) File with the input crystal structure."
    parser.add_argument(
        "-i",
        "--crystal-structure",
        type=read_crystal_structure,
        help=crystal_structure,
        required=True,
    )

    # Optional:
    calculation_presets = "Preset group of tags and default values to use"
    parser.add_argument(
        "-pre",
        "--calculation-presets",
        choices=["scf", "relax", "vc-relax"],
        default=None,
        help=calculation_presets,
    )

    custom_settings_file = "JSON file with custom DFT settings to use"
    parser.add_argument(
        "-file",
        "--custom-settings-file",
        default=None,
        help=custom_settings_file,
    )

    custom_settings_dict = """JSON string with a dictionary of custom DFT
    settings to use. Example: '{"pseudo_dir": "/path/to/pseudo_dir/"}'"""
    parser.add_argument(
        "-dict",
        "--custom-settings-dict",
        default="{}",
        type=json.loads,
        help=custom_settings_dict,
    )

    specify_potentials = "Specify a potential for every chemical species?"
    parser.add_argument(
        "-pot",
        "--specify-potentials",
        type=bool,
        default=False,
        help=specify_potentials,
    )

    write_location = "Directory to write the input file(s) in"
    parser.add_argument("-loc", "--write-location", help=write_location)

    pwx_input_file = "Name of the pw.x input file"
    parser.add_argument("-o", "--pwx-input-file", help=pwx_input_file)

    return parser


def generate_pwx_input_files(sys_args):
    parser = get_parser()
    args = parser.parse_args(sys_args)

    pwig = PwxInputGenerator(
        crystal_structure=args.crystal_structure,
        calculation_presets=args.calculation_presets,
        custom_sett_file=args.custom_settings_file,
        custom_sett_dict=args.custom_settings_dict,
        specify_potentials=args.specify_potentials,
        write_location=args.write_location,
        pwx_input_file=args.pwx_input_file,
    )
    pwig.write_input_files()


if __name__ == "__main__":
    """
    When run as a script, this module will generate input files to use with
    pw.x, for a specified crystal structure, calculation presets, and any
    custom DFT settings on top of preset defaults.

    For a list of optional arguments, see `get_parser()` or run this script
    with "-h" as an argument.
    """
    generate_pwx_input_files(sys.argv[1:])  # pragma: no cover