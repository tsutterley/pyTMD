=========
arguments
=========

- Calculates the nodal corrections for tidal constituents
- Originally based on Richard Ray's ``ARGUMENTS`` fortran subroutine

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.arguments
    pu,pf,G = pyTMD.arguments(MJD, constituents,
        deltat=DELTAT, corrections=CORRECTIONS)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/arguments.py

.. autofunction:: pyTMD.arguments

.. autofunction:: pyTMD.doodson_number

.. autofunction:: pyTMD.arguments._arguments_table

.. autofunction:: pyTMD.arguments._minor_table
