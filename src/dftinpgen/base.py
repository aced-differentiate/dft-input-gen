import os
import six
import abc
from abc import abstractproperty
from abc import abstractmethod

import numpy as np
import ase
from ase import io


class DftInputGeneratorError(Exception):
    """Base class for errors associated with DFT input files generation."""
    pass


@six.add_metaclass(abc.ABCMeta)
class DftInputGenerator(object):
    """
    Base class (interface) to model input generators for specific DFT
    codes after.
    """

    def __init__(self, crystal_structure=None, base_recipe=None,
                 custom_sett_file=None, custom_sett_dict=None,
                 write_location=None, overwrite_files=None, **kwargs):
        """
        Constructor.

        Parameters
        ----------

        crystal_structure: :class:`ase.Atoms` object
            :class:`ase.Atoms` object resulting from `ase.io.read([crystal
            structure file])`.

        base_recipe: str, optional
            The "base" calculation settings to use--must be one of the
            pre-defined recipes provided for `self.dft_package`.

            Pre-defined recipes are in
            [INSTALL_PATH]/[dft_package]/settings/base_recipes/[recipe].json

            For example, if `dft_package` = "vasp", `base_recipe` = "scf", the
            settings in "dftinpgen/vasp/settings/base_recipes/scf.json" are
            used.

        custom_sett_file: str, optional
            Location of a JSON file with custom calculation settings as a
            dictionary of tags and values.

            NB: Custom settings specified here always OVERRIDE those in
            `base_recipe` in case of overlap.

        custom_sett_dict: dict, optional
            Dictionary with custom calculation settings as tags and values.

            NB: Custom settings specified here always OVERRIDE those in
            `base_recipe` and `custom_sett_file`.

            Default: {}

        write_location: str, optional
            Path to the directory in which to write the input files.

            Default: current working directory.

        overwrite_files: bool, optional
            To overwrite files or not, that is the question.

            Default: True

        **kwargs:
            Arbitrary keyword arguments.

        """

        self._crystal_structure = None
        self.crystal_structure = crystal_structure

        self._base_recipe = None
        self.base_recipe = base_recipe

        self._custom_sett_file = None
        self.custom_sett_file = custom_sett_file

        self._custom_sett_dict = {}
        if custom_sett_dict is not None:
            self.custom_sett_dict = custom_sett_dict

        self._write_location = None
        self.write_location = write_location

        self._overwrite_files = True
        if overwrite_files is not None:
            self.overwrite_files = overwrite_files

    @property
    def crystal_structure(self):
        return self._crystal_structure

    @crystal_structure.setter
    def crystal_structure(self, crystal_structure):
        if not isinstance(crystal_structure, ase.Atoms):
            msg = 'Expected type `ase.Atoms`; found {}'.format(type(crystal_structure))
            raise TypeError(msg)
        self._crystal_structure = crystal_structure

    @staticmethod
    def read_crystal_structure(crystal_structure, **kwargs):
        if isinstance(crystal_structure, six.string_types):
            return io.read(crystal_structure, **kwargs)
        else:
            msg = 'Expected type str; found {}'.format(type(crystal_structure))
            raise TypeError(msg)

    @property
    def base_recipe(self):
        return self._base_recipe

    @base_recipe.setter
    def base_recipe(self, base_recipe):
        self._base_recipe = base_recipe

    @property
    def custom_sett_file(self):
        return self._custom_sett_file

    @custom_sett_file.setter
    def custom_sett_file(self, custom_sett_file):
        self._custom_sett_file = custom_sett_file

    @property
    def custom_sett_dict(self):
        return self._custom_sett_dict

    @custom_sett_dict.setter
    def custom_sett_dict(self, custom_sett_dict):
        self._custom_sett_dict = custom_sett_dict

    @property
    def write_location(self):
        return self._write_location

    @write_location.setter
    def write_location(self, write_location):
        if write_location is None:
            self._write_location = os.getcwd()
        else:
            self._write_location = write_location

    @property
    def overwrite_files(self):
        return self._overwrite_files

    @overwrite_files.setter
    def overwrite_files(self, overwrite_files):
        self._overwrite_files = overwrite_files

    def get_kpoint_grid_from_spacing(self, spacing):
        """Returns a list [k1, k2, k3] with the dimensions of a uniform
        k-point grid corresponding to the input `spacing`."""
        if not self.crystal_structure:
            msg = 'Crystal structure not specified'
            raise DftInputGeneratorError(msg)
        rcell = 2*np.pi*(np.linalg.inv(self.crystal_structure.cell).T)
        return list(map(int, np.ceil(np.linalg.norm(rcell, axis=1)/spacing)))

    @abstractproperty
    def dft_package(self):
        """Name of the DFT package input files are generated for."""
        return

    @abstractproperty
    def calculation_settings(self):
        """Dictionary of keywords and values used to generate input files."""
        return

    @abstractmethod
    def write_input_files(self):
        return
