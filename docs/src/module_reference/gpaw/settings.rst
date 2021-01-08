.. _sssec-gpaw-input-settings:

Input settings
++++++++++++++

All input settings have been parsed from the `GPAW manual`_
(last updated January 2021).
The calculator parameter names are stored in `tags_and_groups.json`_.
These parameter names are then made available to the user via
a module level variable ``GPAW_TAGS``.

The settings module also makes available a few sets of default settings 
to be used for common DFT calculations such as ``surface_relax``,
and ``bulk_opt``.
The defaults are stored in JSON files in the `calculation_presets`_ module.
These can be accessed by the user via a module level variable ``GPAW_PRESETS``.
Note that these presets are only reasonable defaults and are not meant to be
prescriptive.

.. _`GPAW manual`: https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html#parameters
.. _`tags_and_groups.json`: https://github.com/CitrineInformatics/dft-input-gen/blob/master/src/dftinputgen/qe/settings/tags_and_groups.json 
.. _`calculation_presets`: https://github.com/CitrineInformatics/dft-input-gen/tree/master/src/dftinputgen/qe/settings/calculation_presets

.. automodule:: dftinputgen.gpaw.settings
    :members:
    :undoc-members:
    :special-members:

.. automodule:: dftinputgen.gpaw.settings.calculation_presets
    :members:
    :undoc-members:
    :special-members:

