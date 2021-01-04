# fmt: off
from gpaw import GPAW
from ase.io import read
from ase.io import write
from ase.optimize import BFGS
from ase.eos import EquationOfState
from fractions import Fraction
import numpy as np
import glob


a = glob.glob('input.traj')
slab = read(a[-1])


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

slab.calc = GPAW(
h=0.16,
kpts={'size': [12, 12, 6]},
occupations={'name': 'fermi-dirac', 'width': 0.05},
xc='BEEF-vdW'
)
bulk_opt_hcp(slab)
