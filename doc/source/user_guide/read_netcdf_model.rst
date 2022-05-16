read_netcdf_model.py
====================

- Reads netCDF format tidal solutions provided by Ohio State University and ESR
- Spatially interpolates tidal constituents to input coordinates

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.read_netcdf_model import read_netcdf_model
    amp,ph,D,c = read_netcdf_model(ilon, ilat, grid_file, model_files,
       type='z', method='spline')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/read_netcdf_model.py

.. autofunction:: pyTMD.read_netcdf_model.extract_netcdf_constants

.. autofunction:: pyTMD.read_netcdf_model.read_netcdf_grid

.. autofunction:: pyTMD.read_netcdf_model.read_elevation_file

.. autofunction:: pyTMD.read_netcdf_model.read_transport_file

.. autofunction:: pyTMD.read_netcdf_model.extend_array

.. autofunction:: pyTMD.read_netcdf_model.extend_matrix
