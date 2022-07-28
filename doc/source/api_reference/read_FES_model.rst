read_FES_model.py
=================

- Reads files for Finite Element Solution (FES) models
   * ascii
   * netCDF4
- Spatially interpolates tidal constituents to input coordinates

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.read_FES_model import extract_FES_constants
    amp,ph = extract_FES_constants(ilon, ilat, model_files,
       type='z',
       version=version,
       method='spline',
       compressed=True,
       scale=1.0/100.0)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/read_FES_model.py

.. autofunction:: pyTMD.read_FES_model.extract_FES_constants

.. autofunction:: pyTMD.read_FES_model.read_ascii_file

.. autofunction:: pyTMD.read_FES_model.read_netcdf_file
