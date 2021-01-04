# fmt: off
from gpaw import GPAW
from ase.io import read
from ase.io import write
from ase.optimize import BFGS
from ase.eos import EquationOfState
from fractions import Fraction
import numpy as np
import glob

from gpaw import restart

a = glob.glob('C3N4.gpw')
slab,calc = restart(a[-1])


def relax(atoms, fmax=0.05, step=0.04):
    name = atoms.get_chemical_formula(mode='hill')
    atoms.calc.set(txt='output.txt')
    atoms.calc.attach(atoms.calc.write, 5, 'output.gpw')
    dyn = BFGS(atoms=atoms, trajectory='output.traj', logfile='qn.log', maxstep=step)
    dyn.run(fmax=fmax)

slab.calc = GPAW(
h=0.16,
kpts={'size': [4, 4, 1]},
occupations={'name': 'fermi-dirac', 'width': 0.05},
poissonsolver={'dipolelayer': 'xy'},
xc='BEEF-vdW'
)
relax(slab)
