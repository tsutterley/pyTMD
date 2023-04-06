=====
astro
=====

- Computes the basic astronomical mean longitudes: `s`, `h`, `p`, `N` and `PP`
- Note `N` is not `N'`, i.e. `N` is decreasing with time.
- Formulae for the period 1990--2010 were derived by David Cartwright

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.astro
    s,h,p,N,PP = pyTMD.astro.mean_longitudes(MJD, ASTRO5=True)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/astro.py

.. autofunction:: pyTMD.mean_longitudes

.. autofunction:: pyTMD.astro.polynomial_sum

.. autofunction:: pyTMD.solar_ecef

.. autofunction:: pyTMD.lunar_ecef
