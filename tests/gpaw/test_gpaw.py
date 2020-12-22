"""Unit tests for the `GPAWInputGenerator` class."""

import os
import pytest

from ase import io as ase_io

from dftinputgen.gpaw.gpaw import GPAWInputGenerator, GPAWInputGeneratorError

test_data_dir = os.path.join(os.path.dirname(__file__), "files")
cu_bulk_struct = ase_io.read(os.path.join(test_data_dir, "cu_bulk.traj"))
with open(os.path.join(test_data_dir, "TEST_bulk_opt.py"), "r") as fr:
    bulk_opt_in = fr.read()
with open(os.path.join(test_data_dir, "TEST_relax.py"), "r") as fr:
    relax_in = fr.read()
with open(os.path.join(test_data_dir, "TEST_relax_custom.py"), "r") as fr:
    relax_custom_in = fr.read()
with open(os.path.join(test_data_dir, "TEST_defaults.py"), "r") as fr:
    defaults_in = fr.read()


def test_dft_package_name():
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct)
    assert gig.dft_package == "GPAW"


def test_gpaw_input_file():
    # default "[calculation_presets]_in.py"
    gig = GPAWInputGenerator(
        crystal_structure=cu_bulk_struct, calculation_presets="bulk_opt"
    )
    assert gig.gpaw_input_file == "bulk_opt_in.py"
    # otherwise use any user-input file name
    gig = GPAWInputGenerator(
        crystal_structure=cu_bulk_struct, calculation_presets="bulk_opt"
    )
    gig.gpaw_input_file = "test.py"
    assert gig.gpaw_input_file == "test.py"


def test_bulk_opt_calculation_presets_settings():
    # Tests bulk_opt settings (fcc/bcc)
    gig = GPAWInputGenerator(
        crystal_structure=cu_bulk_struct, calculation_presets="bulk_opt"
    )
    cs = gig.calculation_settings
    assert cs["calculation"] == "bulk_opt"
    assert cs["kpts"]["size"] == [12, 12, 12]
    assert cs["xc"] == "BEEF-vdW"
    assert cs["h"] == 0.16
    assert cs["occupations"]["name"] == "fermi-dirac"
    assert cs["occupations"]["width"] == 0.05


def test_bulk_opt_hcp_calculation_presets_settings():
    # Tests bulk_opt_hcp settings (hcp)
    gig = GPAWInputGenerator(
        crystal_structure=cu_bulk_struct, calculation_presets="bulk_opt_hcp"
    )
    cs = gig.calculation_settings
    assert cs["calculation"] == "bulk_opt_hcp"
    assert cs["kpts"]["size"] == [12, 12, 6]
    assert cs["xc"] == "BEEF-vdW"
    assert cs["h"] == 0.16
    assert cs["occupations"]["name"] == "fermi-dirac"
    assert cs["occupations"]["width"] == 0.05


def test_molecule_calculation_presets_settings():
    # Tests molecule settings
    gig = GPAWInputGenerator(
        crystal_structure=cu_bulk_struct, calculation_presets="molecule"
    )
    cs = gig.calculation_settings
    assert cs["xc"] == "BEEF-vdW"
    assert cs["h"] == 0.16
    assert cs["occupations"]["name"] == "fermi-dirac"
    assert cs["occupations"]["width"] == 0.05
    assert cs["calculation"] == "relax"


def test_relax_calculation_presets_settings():
    # Tests relax presets
    gig = GPAWInputGenerator(
        crystal_structure=cu_bulk_struct, calculation_presets="relax"
    )
    cs = gig.calculation_settings
    assert cs["xc"] == "BEEF-vdW"
    assert cs["h"] == 0.16
    assert cs["occupations"]["name"] == "fermi-dirac"
    assert cs["occupations"]["width"] == 0.05
    assert cs["calculation"] == "relax"
    assert cs["poissonsolver"]["dipolelayer"] == "xy"


def test_calc_obj_as_str():
    # Tests generation of calculator object as string
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct)
    assert gig.calc_obj_as_str == "slab.calc = GPAW(\n\n)"
    gig.calculation_presets = "bulk_opt"
    calc_obj = "\n".join(bulk_opt_in.splitlines()[35:-1])
    assert gig.calc_obj_as_str == calc_obj


def test_gpaw_input_as_str():
    # Test defaults
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct)
    assert gig.gpaw_input_as_str == "\n".join(defaults_in.splitlines()[1:])
    # Test generation of full bulk opt script
    gig = GPAWInputGenerator(
        crystal_structure=cu_bulk_struct, calculation_presets="bulk_opt"
    )
    assert gig.gpaw_input_as_str == "\n".join(bulk_opt_in.splitlines()[1:])
    gig = GPAWInputGenerator(
        crystal_structure=cu_bulk_struct, calculation_presets="relax"
    )
    assert gig.gpaw_input_as_str == "\n".join(relax_in.splitlines()[1:])


def test_write_gpaw_input():
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct)
    gig.calculation_presets = "relax"
    # no `write_location` input: error
    with pytest.raises(GPAWInputGeneratorError, match="Location to write"):
        gig.write_gpaw_input()
    # no input filename: error
    with pytest.raises(GPAWInputGeneratorError, match="file to write"):
        gig.write_gpaw_input(write_location="/path/to/write_location")

    import tempfile

    _tmp_file = tempfile.NamedTemporaryFile(mode="w", delete=True)
    filename = _tmp_file.name
    write_location = os.path.dirname(filename)
    gig.custom_sett_dict = {
        "kpts": {"size": [6, 6, 1]},
        "xc": "PBE",
        "h": 0.18,
    }
    gig.write_gpaw_input(write_location=write_location, filename=filename)
    with open(filename, "r") as fr:
        assert fr.read() == "\n".join(relax_custom_in.splitlines()[1:])


def test_write_input_files():
    import tempfile

    _tmp_file = tempfile.NamedTemporaryFile(mode="w", delete=True)
    filename = _tmp_file.name
    write_location = os.path.dirname(filename)
    gig = GPAWInputGenerator(crystal_structure=cu_bulk_struct)
    gig.calculation_presets = "relax"
    gig.custom_sett_dict = {
        "kpts": {"size": [6, 6, 1]},
        "xc": "PBE",
        "h": 0.18,
    }
    gig.write_location = write_location
    gig.gpaw_input_file = filename
    gig.write_input_files()
    with open(filename, "r") as fr:
        assert fr.read() == "\n".join(relax_custom_in.splitlines()[1:])
