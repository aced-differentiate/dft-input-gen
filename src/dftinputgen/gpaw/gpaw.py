import os
import json
import six
import itertools

from dftinpgen.data import STANDARD_ATOMIC_WEIGHTS
from dftinpgen.utils import get_elem_symbol
from dftinpgen.gpaw.settings import GPAW_TAGS
from dftinpgen.gpaw.settings.base_recipes import GPAW_BASE_RECIPES

from dftinpgen.base import DftInputGenerator
from dftinpgen.base import DftInputGeneratorError

class GPAWInputGeneratorError(DftInputGeneratorError):
    pass

class GPAWInputGenerator(DftInputGenerator):
    """Base class to generate input python scripts for GPAW """

    def __init__(self, crystal_structure=None, base_recipe=None,
                 custom_sett_file=None, custom_sett_dict=None,
                 write_location=None, overwrite_files=None, **kwargs):
        """
        """
        super(GPAWInputGenerator, self).__init__(
            crystal_structure=crystal_structure,
            base_recipe=base_recipe,
            custom_sett_file=custom_sett_file,
            custom_sett_dict=custom_sett_dict,
            write_location=write_location,
            overwrite_files=overwrite_files,
            **kwargs)

    @property
    def dft_package(self):
        return 'GPAW'


    @property
    def calculation_settings(self):
        calc_sett = {}
        if self.base_recipe is not None:
            calc_sett.update(GPAW_BASE_RECIPES[self.base_recipe])
        if self.custom_sett_file is not None:
            with open(self.custom_sett_file, 'r') as fr:
                calc_sett_update(json.load(fr))
        if self.custom_sett_dict is not None:
            calc_sett.update(self.custom_sett_dict)
        return calc_sett


    def _get_default_input_filename(self):
        return '{}_in.py'.format(self.base_recipe) \
            if self.base_recipe is not None else 'gpaw_in.py'

    @property
    def calc_obj_as_str(self):
        top = "slab.calc = GPAW("
        
        calc_sett = self.calculation_settings

        params = []
        for p in GPAW_TAGS['parameters']:
            if p in calc_sett:
                if type(calc_sett[p]) is str:
                    params.append(f"{p}='{str(calc_sett[p])}'")
                else:
                    params.append(f"{p}={str(calc_sett[p])}")

        return '\n'.join([top,',\n'.join(params),')'])


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
a = glob.glob('input.traj')
slab = read(a[-1])
"""

        calc_sett = self.calculation_settings
#
#        gpaw_skel=f"""slab.calc=GPAW(
#            xc='{calc_sett['xc']}',
#            h={calc_sett['h']},
#            occupations={str(calc_sett['occupations'])},
#            poissonsolver={str(calc_sett['poissonsolver'])}
#        )
#        """

        calc_type = calc_sett['calculation']
        if calc_type == 'relax':
            # This definition will become unnecessary when defined in DFT FLOW
            define_relax_fn="""
def relax(atoms, fmax=0.05, step=0.04):
    name = atoms.get_chemical_formula(mode='hill')
    atoms.calc.set(txt='output.txt')
    atoms.calc.attach(atoms.calc.write, 5, 'output.gpw')
    dyn = BFGS(atoms=atoms, trajectory='output.traj', logfile='qn.log', maxstep=step)
    dyn.run(fmax=fmax)
"""
            return '\n'.join([header,read_init_traj,define_relax_fn,self.calc_obj_as_str,'relax(slab)'])

        elif calc_type == 'bulk_opt' or calc_type == 'bulk_opt_hcp':
            define_bulk_opt_fn="""
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
            return '\n'.join([header,read_init_traj,define_bulk_opt_fn,self.calc_obj_as_str,'bulk_opt(slab)'])
        # if not a relax or bulk_opt calculation, defaults to getting total energy of static structure
        return '\n'.join([header,read_init_traj,define_relax_fn,calc_obj_as_str,'slab.get_total_energy()'])

    def write_gpaw_input(self, write_location=None, filename=None):
        if not self.gpaw_input_as_str.strip():
            msg = 'Nothing to write (probably no input settings found.)'
            raise GPAWInputGeneratorError(msg)
        if write_location is None:
            msg = 'Location to write files not specified'
            raise GPAWInputGeneratorError(msg)
        if filename is None:
            msg = 'Name of the input file to write into not specified'
            raise GPAWInputGeneratorError(msg)
        with open(os.path.join(write_location, filename), 'w') as fw:
            fw.write(self.gpaw_input_as_str)

    def write_input_files(self):
        self.write_gpaw_input(
            write_location = self.write_location,
            filename = self._get_default_input_filename(),
        )
