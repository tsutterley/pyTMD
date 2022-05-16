read_tide_model.py
==================

- Reads OTIS format tidal solutions provided by Ohio State University and ESR
  * multi-constituent binary
  * ATLAS-compact binary
  * single-constituent binary
- Spatially interpolates tidal constituents to input coordinates

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.read_tide_model import extract_tidal_constants
    amp,ph,D,c = extract_tidal_constants(ilon, ilat, grid_file, model_file, EPSG,
        TYPE='z', METHOD='spline', GRID='OTIS')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/read_tide_model.py

.. autofunction:: pyTMD.read_tide_model.extract_tidal_constants

.. autofunction:: pyTMD.read_tide_model.read_tide_grid

.. autofunction:: pyTMD.read_tide_model.read_atlas_grid

.. autofunction:: pyTMD.read_tide_model.read_netcdf_grid

.. autofunction:: pyTMD.read_tide_model.read_constituents

.. autofunction:: pyTMD.read_tide_model.read_elevation_file

.. autofunction:: pyTMD.read_tide_model.read_atlas_elevation

.. autofunction:: pyTMD.read_tide_model.read_transport_file

.. autofunction:: pyTMD.read_tide_model.read_atlas_transport

.. autofunction:: pyTMD.read_tide_model.create_atlas_mask

.. autofunction:: pyTMD.read_tide_model.interpolate_atlas_model

.. autofunction:: pyTMD.read_tide_model.combine_atlas_model

.. autofunction:: pyTMD.read_tide_model.read_netcdf_file

.. autofunction:: pyTMD.read_tide_model.extend_array

.. autofunction:: pyTMD.read_tide_model.extend_matrix
