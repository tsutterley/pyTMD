===========
convert_crs
===========

- Uses `pyproj <https://pyproj4.github.io/pyproj>`_ to convert points to and from the Coordinates Reference Systems (CRS) of tide models

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.convert_crs
    x, y = pyTMD.convert_crs(lon, lat, PROJ, 'F')
    lon, lat = pyTMD.convert_crs(x, y, PROJ, 'B')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/convert_crs.py

.. autofunction:: pyTMD.convert_crs

.. autofunction:: pyTMD.convert_crs.crs_from_input

.. autofunction:: pyTMD.convert_crs._EPSG3031

.. autofunction:: pyTMD.convert_crs._EPSG3413

.. autofunction:: pyTMD.convert_crs._CATS2008

.. autofunction:: pyTMD.convert_crs._EPSG3976

.. autofunction:: pyTMD.convert_crs._PSNorth

.. autofunction:: pyTMD.convert_crs._EPSG4326

.. autofunction:: pyTMD.convert_crs._custom
