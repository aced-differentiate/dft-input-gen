# fmt: off
from gpaw import GPAW
from ase.io import read
from ase.io import write
from ase.optimize import BFGS
from ase.eos import EquationOfState
import numpy as np
import glob


a = glob.glob('input.traj')
slab = read(a[-1])

slab.calc = GPAW(

)
slab.get_total_energy()
