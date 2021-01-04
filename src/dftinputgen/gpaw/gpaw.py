import os
import json
import six
import itertools

from dftinputgen.data import STANDARD_ATOMIC_WEIGHTS
from dftinputgen.utils import get_elem_symbol
from dftinputgen.gpaw.settings import GPAW_TAGS
from dftinputgen.gpaw.settings.calculation_presets import GPAW_PRESETS

from dftinputgen.base import DftInputGenerator
from dftinputgen.base import DftInputGeneratorError


class GPAWInputGeneratorError(DftInputGeneratorError):
    pass


class GPAWInputGenerator(DftInputGenerator):
    """
    Constructor for the GPAW class.

    Parameters
    ----------

    crystal_structure: :class:`ase.Atoms` object
        :class:`ase.Atoms` object from `ase.io.read([crystal structure
        file])`.

    calculation_presets: str, optional
        The "base" calculation settings to use--must be one of the
        pre-defined groups of tags and values provided for pw.x.

        Pre-defined settings for some common calculation types are in
        INSTALL_PATH/qe/settings/calculation_presets/

    custom_sett_file: str, optional
        Location of a JSON file with custom calculation settings as a
        dictionary of tags and values.

        NB: Custom settings specified here always OVERRIDE those in
        `calculation_presets` in case of overlap.

    custom_sett_dict: dict, optional
        Dictionary with custom calculation settings as tags and values/

        NB: Custom settings specified here always OVERRIDE those in
        `calculation_presets` and `custom_sett_file`.

    write_location: str, optional
        Path to the directory in which to write the input file(s).

        Default: Current working directory.

    gpaw_input_file: str, optional
        Name of the file in which to write the GPAW python script

        Default: "[`calculation_presets`]_in.py" if `calculation_presets` is
        specified by the user, else "gpaw_in.py".

    struct_filename: str, optional
        Name of the structure file readable by `ase.io.read` that the written
        gpaw script will call from.

        Default: "input.traj"

    overwrite_files: bool, optional
        To overwrite files or not, that is the question.

        Default: True

    **kwargs:
        Arbitrary keyword arguments.

    """

    def __init__(
        self,
        crystal_structure=None,
        calculation_presets=None,
        custom_sett_file=None,
        custom_sett_dict=None,
        write_location=None,
        gpaw_input_file=None,
        struct_filename=None,
        overwrite_files=None,
        **kwargs,
    ):
        """"""
        super(GPAWInputGenerator, self).__init__(
            crystal_structure=crystal_structure,
            calculation_presets=calculation_presets,
            custom_sett_file=custom_sett_file,
            custom_sett_dict=custom_sett_dict,
            write_location=write_location,
            overwrite_files=overwrite_files,
            **kwargs,
        )

        self._calculation_settings = self._get_calculation_settings()

        self._gpaw_input_file = self._get_default_input_filename()
        self.gpaw_input_file = gpaw_input_file

        self._struct_filename = "input.traj"
        self.struct_filename = struct_filename

    @property
    def dft_package(self):
        return "GPAW"

    @property
    def gpaw_input_file(self):
        """Name of the gpaw input file to write to."""
        return self._gpaw_input_file

    @gpaw_input_file.setter
    def gpaw_input_file(self, gpaw_input_file):
        if gpaw_input_file is not None:
            self._gpaw_input_file = gpaw_input_file

    @property
    def struct_filename(self):
        """Name of the structure file that the gpaw script will read from."""
        return self._struct_filename

    @struct_filename.setter
    def struct_filename(self, struct_filename):
        if struct_filename is not None:
            self._struct_filename = struct_filename

    @property
    def calculation_settings(self):
        """Dictionary of all calculation settings to use as input for gpaw."""
        return self._get_calculation_settings()

    def _get_calculation_settings(self):
        """Load all calculation settings: user-input and auto-determined."""
        calc_sett = {}
        if self.calculation_presets is not None:
            calc_sett.update(GPAW_PRESETS[self.calculation_presets])
        if self.custom_sett_from_file is not None:
            calc_sett.update(self.custom_sett_from_file)
        if self.custom_sett_dict is not None:
            calc_sett.update(self.custom_sett_dict)
        return calc_sett

    def _get_default_input_filename(self):
        if self.calculation_presets is None:
            return "gpaw_in.py"
        return "{}_in.py".format(self.calculation_presets)

    @property
    def calculator_object_as_str(self):
        top = "slab.calc = GPAW("

        calc_sett = self.calculation_settings

        params = []
        for p in GPAW_TAGS["parameters"]:
            if p in calc_sett:
                if type(calc_sett[p]) is str:
                    params.append(f"{p}='{str(calc_sett[p])}'")
                else:
                    params.append(f"{p}={str(calc_sett[p])}")

        return "\n".join([top, ",\n".join(params), ")"])

    @property
    def gpaw_input_as_str(self):
        header = """from gpaw import GPAW
from ase.io import read
from ase.io import write
from ase.optimize import BFGS
from ase.eos import EquationOfState
from fractions import Fraction
import numpy as np
import glob
"""

        read_init_traj = """
a = glob.glob('{}')
slab = read(a[-1])
""".format(
            self.struct_filename
        )

        calc_sett = self.calculation_settings
        try:
            calc_type = calc_sett["calculation"]
        except KeyError:
            # if no input settings found will return calc object without any settings defined
            calc_type = None
        if calc_type == "relax":
            define_relax_fn = """
def relax(atoms, fmax=0.05, step=0.04):
    name = atoms.get_chemical_formula(mode='hill')
    atoms.calc.set(txt='output.txt')
    atoms.calc.attach(atoms.calc.write, 5, 'output.gpw')
    dyn = BFGS(atoms=atoms, trajectory='output.traj', logfile='qn.log', maxstep=step)
    dyn.run(fmax=fmax)
"""
            return "\n".join(
                [
                    header,
                    read_init_traj,
                    define_relax_fn,
                    self.calculator_object_as_str,
                    "relax(slab)",
                ]
            )

        elif calc_type == "bulk_opt":
            define_bulk_opt_fn = """
def bulk_opt(atoms, step=0.05):
    cell = atoms.get_cell()
    name = atoms.get_chemical_formula(mode='hill')
    vol=atoms.get_volume()
    volumes =[]
    energies=[]
    for x in np.linspace(1-2*step,1+2*step,5):
        atoms.set_cell(cell*x, scale_atoms=True)
        atoms.calc.set(txt=name+'_'+str(x)+'.txt')
        energies.append(atoms.get_potential_energy())
        volumes.append(atoms.get_volume())
    eos = EquationOfState(volumes, energies)
    v0,e0,B= eos.fit()
    atoms.set_cell((v0/vol)**Fraction('1/3')*cell,scale_atoms=True)
    x0=(v0/vol)**Fraction('1/3')
    atoms.calc.set(txt='output.txt')
    dyn=BFGS(atoms=atoms,trajectory='output.traj',logfile = 'qn.log')
    dyn.run(fmax=0.05)
    atoms.calc.write('output.gpw')
"""
            return "\n".join(
                [
                    header,
                    read_init_traj,
                    define_bulk_opt_fn,
                    self.calculator_object_as_str,
                    "bulk_opt(slab)",
                ]
            )
        elif calc_type == "bulk_opt_hcp":
            define_bulk_opt_hcp_fn = """
def optimize_bulk_hcp(atoms, step=0.05, nstep=5):
    # Get initial cell
    cell = atoms.get_cell()
    a_in = cell[0]
    b_in = cell[1]
    c_in = cell[2]
    name = atoms.get_chemical_formula(mode='hill')
    vol=atoms.get_volume()

    # Optimize in ab first
    volumes_ab =[]
    energies_ab=[]
    scale_factors = np.linspace(1-2*step,1+2*step,nstep)
    for x in scale_factors:
        atoms.set_cell([a_in*x,b_in*x,c_in], scale_atoms=True)
        atoms.calc.set(txt=name+'_'+'ab'+'_'+str(x)+'.txt')
        energies_ab.append(atoms.get_potential_energy())
        volumes_ab.append(atoms.get_volume())

    # Fit in ab
    eos = EquationOfState(volumes_ab, energies_ab)
    v0,e0,B= eos.fit()
    a_fit = a_in*np.sqrt(v0/vol)
    b_fit = b_in*np.sqrt(v0/vol)

    # Generate new cell with fit ab parameters
    atoms.set_cell([a_fit,b_fit,c_in], scale_atoms=True)
    vol_ab = atoms.get_volume()

    # Optimize c
    energies_c = []
    volumes_c = []
    for x in scale_factors:
        atoms.set_cell([a_fit,b_fit,c_in*x], scale_atoms=True)
        atoms.calc.set(txt=name+'_'+'c'+'_'+str(x)+'.txt')
        energies_c.append(atoms.get_potential_energy())
        volumes_c.append(atoms.get_volume())

    # Fit in c
    eos_c = EquationOfState(volumes_c,energies_c,eos='birchmurnaghan')
    v0, e0, B = eos_c.fit()
    c_scaling = v0/vol_ab
    c_fit = c_in*c_scaling

    # Update cell with optimized parameters and relax atoms
    atoms.set_cell([a_fit,b_fit,c_fit], scale_atoms=True)
    atoms.calc.set(txt='output.txt')
    dyn=BFGS(atoms=atoms,trajectory='output.traj',logfile = 'qn.log')
    dyn.run(fmax=0.05)
    atoms.calc.write("output.gpw")
"""
            return "\n".join(
                [
                    header,
                    read_init_traj,
                    define_bulk_opt_hcp_fn,
                    self.calculator_object_as_str,
                    "bulk_opt_hcp(slab)",
                ]
            )

        # if not a relax or bulk_opt calculation, defaults to getting total energy of static structure
        return "\n".join(
            [
                header,
                read_init_traj,
                self.calculator_object_as_str,
                "slab.get_total_energy()",
            ]
        )

    def write_gpaw_input(self, write_location=None, filename=None):
        if write_location is None:
            msg = "Location to write files not specified"
            raise GPAWInputGeneratorError(msg)
        if filename is None:
            msg = "Name of the input file to write into not specified"
            raise GPAWInputGeneratorError(msg)
        with open(os.path.join(write_location, filename), "w") as fw:
            fw.write(self.gpaw_input_as_str)

    def write_input_files(self):
        self.write_gpaw_input(
            write_location=self.write_location, filename=self.gpaw_input_file,
        )
