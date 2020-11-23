import os
from setuptools import setup, find_packages


with open(os.path.join(os.path.dirname(__file__),
                       'src', 'dftinpgen', 'VERSION.txt')) as fr:
    version = fr.read().strip()


setup(
    name='dftinpgen',
    version=version,
    description='Unopinionated library to generate input files for DFT codes',
    url='https://github.com/CitrineInformatics/dft-input-gen',
    author='Vinay Hegde',
    author_email='vhegde@citrine.io',
    packages=find_packages(where='src', exclude=['docs']),
    package_dir={'': 'src'},
    install_requires=[
        'six',
        'numpy',
        'ase <= 3.17'
    ],
    include_package_data=True,
)