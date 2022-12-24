======
io.FES
======

- Reads files for Finite Element Solution (FES) models

   * ascii
   * netCDF4
- Spatially interpolates tidal constituents to input coordinates

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.io
    amp,ph = pyTMD.io.FES.extract_constants(ilon, ilat, model_files,
       type='z',
       version=version,
       method='spline',
       compressed=True,
       scale=1.0/100.0)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/io/FES.py

.. autofunction:: pyTMD.io.FES.extract_constants

.. autofunction:: pyTMD.io.FES.read_constants

.. autofunction:: pyTMD.io.FES.interpolate_constants

.. autofunction:: pyTMD.io.FES.read_ascii_file

.. autofunction:: pyTMD.io.FES.read_netcdf_file
