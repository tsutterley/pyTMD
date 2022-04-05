================
convert_ll_xy.py
================

- Uses `pyproj <https://pyproj4.github.io/pyproj>`_ to convert lat/lon points to and from projected coordinates

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.convert_ll_xy import convert_ll_xy
    x,y = convert_ll_xy(lon,lat,PROJ,'F')
    lon,lat = convert_ll_xy(x,y,PROJ,'B')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/convert_ll_xy.py

.. autofunction:: pyTMD.convert_ll_xy
