==================
calc_delta_time.py
==================

- Calculates the difference between dynamic time and universal time (TT - UT1) following Richard Ray's PERTH3 algorithms

   * `http://maia.usno.navy.mil/ser7/deltat.data <http://maia.usno.navy.mil/ser7/deltat.data>`_
   * `ftp://cddis.nasa.gov/products/iers/deltat.data <ftp://cddis.nasa.gov/products/iers/deltat.data>`_

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.calc_delta_time import calc_delta_time
    deltat = calc_delta_time(delta_file, idays)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/calc_delta_time.py

.. autofunction:: pyTMD.calc_delta_time
