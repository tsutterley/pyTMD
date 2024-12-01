=========
arguments
=========

- Calculates the nodal corrections for tidal constituents
- Originally based on Richard Ray's ``ARGUMENTS`` fortran subroutine

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.arguments
    pu,pf,G = pyTMD.arguments.arguments(MJD, constituents,
        deltat=DELTAT, corrections=CORRECTIONS)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/arguments.py

.. autofunction:: pyTMD.arguments.arguments

.. autofunction:: pyTMD.arguments.minor_arguments

.. autofunction:: pyTMD.arguments.coefficients_table

.. autofunction:: pyTMD.arguments.doodson_number

.. autofunction:: pyTMD.arguments.nodal

.. autofunction:: pyTMD.arguments.frequency

.. autofunction:: pyTMD.arguments._arguments_table

.. autofunction:: pyTMD.arguments._minor_table

.. autofunction:: pyTMD.arguments._constituent_parameters

.. autofunction:: pyTMD.arguments._love_numbers

.. autofunction:: pyTMD.arguments._parse_tide_potential_table

.. autofunction:: pyTMD.arguments._to_doodson_number

.. autofunction:: pyTMD.arguments._to_extended_doodson

.. autofunction:: pyTMD.arguments._from_doodson_number

.. autofunction:: pyTMD.arguments._from_extended_doodson
