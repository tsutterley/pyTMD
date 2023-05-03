=====
astro
=====
- Astronomical and nutation routines
- Computes the basic astronomical mean longitudes: `S`, `H`, `P`, `N` and `PP`
- Computes astronomical phase angles for the sun and moon:  `S`, `H`, `P`, `TAU`, `ZNS` and `PS`
- Computes the solar and lunar positions in Earth-Centered Earth-Fixed (ECEF) coordinates

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.astro
    S,H,P,N,PP = pyTMD.astro.mean_longitudes(MJD, ASTRO5=True)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/astro.py

.. autofunction:: pyTMD.astro.polynomial_sum

.. autofunction:: pyTMD.astro.rotate

.. autofunction:: pyTMD.astro.mxm

.. autofunction:: pyTMD.astro.mxmxm

.. autofunction:: pyTMD.astro.mean_longitudes

.. autofunction:: pyTMD.astro.phase_angles

.. autofunction:: pyTMD.astro.delaunay_arguments

.. autofunction:: pyTMD.astro.mean_obliquity

.. autofunction:: pyTMD.astro.solar_ecef

.. autofunction:: pyTMD.astro.solar_ephemerides

.. autofunction:: pyTMD.astro.lunar_ecef

.. autofunction:: pyTMD.astro.lunar_ephemerides

.. autofunction:: pyTMD.astro._parse_table_5_2e

.. autofunction:: pyTMD.astro._parse_table_5_3a

.. autofunction:: pyTMD.astro._parse_table_5_3b
