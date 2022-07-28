read_GOT_model.py
=================

- Reads files for Richard Ray's Global Ocean Tide (GOT) models
- Spatially interpolates tidal constituents to input coordinates

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.read_GOT_model import extract_GOT_constants
    amp,ph,c = extract_GOT_constants(ilon,i lat, model_files,
       method='spline')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/read_GOT_model.pyu

.. autofunction:: pyTMD.read_GOT_model.extract_GOT_constants

.. autofunction:: pyTMD.read_GOT_model.read_GOT_grid

.. autofunction:: pyTMD.read_GOT_model.extend_array

.. autofunction:: pyTMD.read_GOT_model.extend_matrix
