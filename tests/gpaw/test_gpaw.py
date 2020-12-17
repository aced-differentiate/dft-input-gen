"""Unit tests for the `GPAWInputGenerator` class."""

import os
import pytest

from ase import io as ase_io

from dftinputgen.gpaw.gpaw import GPAWInputGenerator

test_data_dir = os.path.join(os.path.dirname(__file__), "files")
cu_bulk_struct = ase_io.read(os.path.join(test_data_dir, "cu_bulk.traj"))
with open(os.path.join(test_data_dir, "TEST_cu_bulk_opt.py"), "r") as fr:
    cu_bulk_opt_in = fr.read()


def test_dft_package_name():
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct)
    assert gig.dft_package == "GPAW"


def test_gpaw_input_file():
    # default "[calculation_presets]_in.py"
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct, calculation_presets="bulk_opt")
    assert gig.gpaw_input_file == "bulk_opt_in.py"
    # otherwise use any user-input file name
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct, calculation_presets="bulk_opt")
    gig.gpaw_input_file = "test.py"
    assert gig.gpaw_input_file == "test.py"


def test_bulk_opt_calculation_presets_settings():
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct, calculation_presets="bulk_opt")
    cs = gig.calculation_settings
    assert cs["calculation"] == "bulk_opt"
    assert cs["kpts"]["size"] == [12,12,12]
    assert cs["xc"] == "BEEF-vdW"
    assert cs["h"] == 0.16
    assert cs["occupations"]["name"] == "fermi-dirac"
    assert cs["occupations"]["width"] == 0.05

