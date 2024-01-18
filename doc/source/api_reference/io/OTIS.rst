====
OTIS
====

- Reads OTIS format tidal solutions provided by Ohio State University and ESR

  * multi-constituent binary
  * ATLAS-compact binary
  * single-constituent binary
- Spatially interpolates tidal constituents to input coordinates
- Writes OTIS-format tide files

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.io
    amp,ph,D,c = pyTMD.io.OTIS.extract_constants(ilon, ilat, grid_file, model_file, EPSG,
        type='z',
        method='spline',
        grid='OTIS')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/io/OTIS.py

.. autofunction:: pyTMD.io.OTIS.extract_constants

.. autofunction:: pyTMD.io.OTIS.read_constants

.. autofunction:: pyTMD.io.OTIS.interpolate_constants

.. autofunction:: pyTMD.io.OTIS.read_otis_grid

.. autofunction:: pyTMD.io.OTIS.read_atlas_grid

.. autofunction:: pyTMD.io.OTIS.read_netcdf_grid

.. autofunction:: pyTMD.io.OTIS.read_constituents

.. autofunction:: pyTMD.io.OTIS.read_otis_elevation

.. autofunction:: pyTMD.io.OTIS.read_atlas_elevation

.. autofunction:: pyTMD.io.OTIS.read_otis_transport

.. autofunction:: pyTMD.io.OTIS.read_atlas_transport

.. autofunction:: pyTMD.io.OTIS.create_atlas_mask

.. autofunction:: pyTMD.io.OTIS.interpolate_atlas_model

.. autofunction:: pyTMD.io.OTIS.combine_atlas_model

.. autofunction:: pyTMD.io.OTIS.read_netcdf_file

.. autofunction:: pyTMD.io.OTIS.output_otis_grid

.. autofunction:: pyTMD.io.OTIS.output_otis_elevation

.. autofunction:: pyTMD.io.OTIS.output_otis_transport

.. autofunction:: pyTMD.io.OTIS.extend_array

.. autofunction:: pyTMD.io.OTIS.extend_matrix

.. autofunction:: pyTMD.io.OTIS._mask_nodes

.. autofunction:: pyTMD.io.OTIS._interpolate_to_nodes
