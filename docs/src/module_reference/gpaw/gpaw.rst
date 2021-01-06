.. _sssec-gpaw:

Input GPAW
++++++++++

The :class:`GPAWInputGenerator <dftinputgen.gpaw.GPAWInputGenerator>` class
(derived from :class:`DftInputGenerator <dftinputgen.base.DftInputGenerator>`)
implements functionality to generate python scripts for the `GPAW`_
package. 

Python scripts are written that use `ASE`_ with GPAW specified via a
calculator object. This calculator object is where all settings for
the DFT calculation is given. An example of the calculator object
that is automatically written is given below.

.. code-block:: python

   slab.calc = GPAW(
   h=0.16,
   kpts={'size': [4, 4, 1]},
   occupations={'name': 'fermi-dirac', 'width': 0.05},
   poissonsolver={'dipolelayer': 'xy'},
   xc='BEEF-vdW'
   )

For each calculator setting, a list of valid parameters is looked up from
a ``GPAW_TAGS`` dictionary which is then defined within the written
calculator object (more information about the valid calculator parameters
and defaults provided as ``calculation_presets`` is :ref:`here
<sssec-gpaw-input-settings>`).

Currently supported calculations include relaxations, bulk optimizations,
and total energy.

**Note:** The scripts are written assuming that the crystal structure
is stored in an `ASE readable format`_ (e.g. traj, cif, etc..) with 
the desired initial magnetic moments defined.

.. _`GPAW`: https://wiki.fysik.dtu.dk/gpaw/index.html
.. _`ASE`: https://wiki.fysik.dtu.dk/ase/
.. _`ASE readable format`: https://wiki.fysik.dtu.dk/ase/ase/io/io.html

Interfaces
==========

.. automodule:: dftinputgen.gpaw.gpaw
    :members:
    :inherited-members:
    :undoc-members:
