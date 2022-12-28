===========
interpolate
===========

- Interpolators for spatial data

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.interpolate
    data = pyTMD.interpolate.bilinear(ilon, ilat, idata, lon, lat)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/interpolate.py

.. autofunction:: pyTMD.interpolate.bilinear

.. autofunction:: pyTMD.interpolate.spline

.. autofunction:: pyTMD.interpolate.regulargrid

.. autofunction:: pyTMD.interpolate.extrapolate
