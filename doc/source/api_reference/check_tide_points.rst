=================
check_tide_points
=================

- Check if points are within a tide model domain
- Can check OTIS format tidal solutions provided by Ohio State University and ESR
- Can check Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
- Can check Finite Element Solution (FES) models provided by AVISO

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.check_tide_points import check_tide_points
    valid = check_tide_points(x, y, DIRECTORY=DIRECTORY,
        MODEL=MODEL, EPSG=3031)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/check_tide_points.py

.. autofunction:: pyTMD.check_tide_points
