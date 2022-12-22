======
io.GOT
======

- Reads files for Richard Ray's Global Ocean Tide (GOT) models
- Spatially interpolates tidal constituents to input coordinates

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.io
    amp,ph,c = pyTMD.io.GOT.extract_constants(ilon,i lat, model_files,
       method='spline')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/io/GOT.py

.. autofunction:: pyTMD.io.GOT.extract_constants

.. autofunction:: pyTMD.io.GOT.read_constants

.. autofunction:: pyTMD.io.GOT.interpolate_constants

.. autofunction:: pyTMD.io.GOT.read_ascii_file

.. autofunction:: pyTMD.io.GOT.extend_array

.. autofunction:: pyTMD.io.GOT.extend_matrix
