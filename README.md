# dftinpgen

[![Build Status](https://travis-ci.com/CitrineInformatics/dft-input-gen.svg?token=qbMA4N9P9kHgFLrLQ51g&branch=master)](https://travis-ci.com/CitrineInformatics/dft-input-gen)

Unopinionated input file generator for DFT codes.


## Requirements

Python >=2.7 or >=3.8, with dependencies listed in
[requirements.txt](https://github.com/CitrineInformatics/dft-input-gen/blob/master/requirements.txt).


## Installation

1. Clone from Github:

```
git clone git@github.com:CitrineInformatics/dft-input-gen.git
```

2. Install package requirements:

```
cd dft-input-gen
pip -r requirements.txt
```

3. Install the package:

```
pip install -e .
```

4. (optional) Unit tests can be run using `pytest`:

```
pip -r test_requirements.txt
pytest -sv
```


## Usage

To generate input files to run an `scf` calculation using `pw.x` for a input
crystal structure in `my_crystal_structure.cif`, do:

```
from dftinpgen.utils import read_crystal_structure
from dftinpgen.qe.pwx import PwxInputGenerator

# read the input crystal into an `ase.Atoms` object
crystal_structure = read_crystal_structure("my_crystal_structure.cif")

# print formatted pw.x input to standard output
pwig = PwxInputGenerator(
   crystal_structure=crystal_structure,
   calculation_presets="scf",
)
print(pwig.pwx_input_as_str)
```

Further details of the API and examples can be found in the package
documentation.


## DFT codes supported

1. [pw.x](https://www.quantum-espresso.org/Doc/INPUT_PW.html) from the
   [Quantum Espresso package](https://www.quantum-espresso.org/)
2. (under development) Post-processing utilities for pw.x:
   [dos.x](https://www.quantum-espresso.org/Doc/INPUT_DOS.html),
   [bands.x](https://www.quantum-espresso.org/Doc/INPUT_BANDS.html),
   [projwfc.x](https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html)


## Contributing

Contributions are welcome, both issues and pull requests.
This project follows the [gitflow
workflow](https://www.atlassian.com/git/tutorials/comparing-workflows#gitflow-workflow); 
please submit all PRs to the `develop` branch.


## Documentation

Further documentation for this package, including a guide for developers, is
available in `docs/src`.
