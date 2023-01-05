======================
calc_astrol_longitudes
======================

- Computes the basic astronomical mean longitudes: `s`, `h`, `p`, `N` and `PP`
- Note `N` is not `N'`, i.e. `N` is decreasing with time.
- Formulae for the period 1990--2010 were derived by David Cartwright

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.calc_astrol_longitudes import calc_astrol_longitudes
    s,h,p,N,PP = calc_astrol_longitudes(MJD, ASTRO5=True)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/calc_astrol_longitudes.py

.. autofunction:: pyTMD.calc_astrol_longitudes

.. autofunction:: pyTMD.calc_astrol_longitudes.polynomial_sum
