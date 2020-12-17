from gpaw import GPAW
from ase.io import read
from ase.io import write
from ase.optimize import BFGS
import numpy as np
import glob


a = glob.glob("input.traj")
slab = read(a[-1])


def bulk_opt(atoms, step=0.05):
    cell = atoms.get_cell()
    name = atoms.get_chemical_formula(mode="hill")
    vol = atoms.get_volume()
    volumes = []
    energies = []
    for x in np.linspace(1 - 2 * step, 1 + 2 * step, 5):
        atoms.set_cell(cell * x, scale_atoms=True)
        atoms.calc.set(txt=name + "_" + str(x) + ".txt")
        energies.append(atoms.get_potential_energy())
        volumes.append(atoms.get_volume())
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    atoms.set_cell((v0 / vol) ** Fraction("1/3") * cell, scale_atoms=True)
    x0 = (v0 / vol) ** Fraction("1/3")
    atoms.calc.set(txt=name + "_" + str(x0) + ".txt")
    dyn = BFGS(
        atoms=atoms,
        trajectory=name + "_" + str(x0) + ".traj",
        logfile="qn.log",
    )
    dyn.run(fmax=0.05)
    atoms.calc.write(name + str(x0) + ".gpw")


slab.calc = GPAW(
    h=0.16,
    kpts={"size": [12, 12, 12]},
    occupations={"name": "fermi-dirac", "width": 0.05},
    xc="BEEF-vdW",
)
bulk_opt(slab)
