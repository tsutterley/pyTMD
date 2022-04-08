=================
iers_mean_pole.py
=================

- Provides the angular coordinates of the IERS Conventional Mean Pole (CMP)
- Based on the coordinates from `IERS mean pole tabulations <ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab>`_
- See `materials from Chapter 7 of the IERS Conventions <https://webtai.bipm.org/iers/convupdt/convupdt_c7.html>`_

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.iers_mean_pole import iers_mean_pole
    x,y,flag = iers_mean_pole(input_file,input_epoch,version,FILL_VALUE=FILL_VALUE)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/iers_mean_pole.py

.. autofunction:: pyTMD.iers_mean_pole
