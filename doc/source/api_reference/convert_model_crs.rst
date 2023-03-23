=================
convert_model_crs
=================

- Uses `pyproj <https://pyproj4.github.io/pyproj>`_ to convert points to and from tide model Coordinates Reference Systems (CRS)

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.convert_model_crs
    x, y = pyTMD.convert_model_crs(lon,lat,PROJ,'F')
    lon, lat = pyTMD.convert_model_crs(x,y,PROJ,'B')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/convert_model_crs.py

.. autofunction:: pyTMD.convert_model_crs

.. autofunction:: pyTMD.convert_model_crs.crs_from_input

.. autofunction:: pyTMD.convert_model_crs.convert_EPSG3031

.. autofunction:: pyTMD.convert_model_crs.convert_EPSG3413

.. autofunction:: pyTMD.convert_model_crs.convert_CATS2008

.. autofunction:: pyTMD.convert_model_crs.convert_EPSG3976

.. autofunction:: pyTMD.convert_model_crs.convert_PSNorth

.. autofunction:: pyTMD.convert_model_crs.convert_EPSG4326

.. autofunction:: pyTMD.convert_model_crs.convert_crs
