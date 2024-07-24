===
GOT
===

- Reads files for Richard Ray's Global Ocean Tide (GOT) models

  * ascii format
  * netcdf format
- Spatially interpolates tidal constituents to input coordinates

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.io
    amp,ph,c = pyTMD.io.GOT.extract_constants(ilon, ilat, model_files,
       method='spline')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/io/GOT.py

.. autofunction:: pyTMD.io.GOT.extract_constants

.. autofunction:: pyTMD.io.GOT.read_constants

.. autofunction:: pyTMD.io.GOT.interpolate_constants

.. autofunction:: pyTMD.io.GOT.read_ascii_file

.. autofunction:: pyTMD.io.GOT.read_netcdf_file

.. autofunction:: pyTMD.io.GOT.output_netcdf_file

.. autofunction:: pyTMD.io.GOT._extend_array

.. autofunction:: pyTMD.io.GOT._extend_matrix

.. autofunction:: pyTMD.io.GOT._crop

.. autofunction:: pyTMD.io.GOT._shift
