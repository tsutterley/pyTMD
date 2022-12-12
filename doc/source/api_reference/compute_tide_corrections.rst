========================
compute_tide_corrections
========================

- Calculates tidal elevations for correcting elevation or imagery data
- Can use OTIS format tidal solutions provided by Ohio State University and ESR
- Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
- Can use Finite Element Solution (FES) models provided by AVISO

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.compute_tide_corrections import compute_tide_corrections
    tide = compute_tide_corrections(x, y, delta_time, DIRECTORY=DIRECTORY,
        MODEL=MODEL, EPOCH=(2000,1,1,0,0,0), EPSG=3031, TYPE='drift')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute_tide_corrections.py

.. autofunction:: pyTMD.compute_tide_corrections
