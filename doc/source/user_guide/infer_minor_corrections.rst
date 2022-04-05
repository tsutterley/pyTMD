==========================
infer_minor_corrections.py
==========================

- Calculates the tidal corrections for minor constituents
- Based on Richard Ray's PERTH3 (PREdict Tidal Heights) algorithms

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.infer_minor_corrections import infer_minor_corrections
    dh = infer_minor_corrections(t, zmajor, constituents,
        DELTAT=DELTAT, CORRECTIONS=CORRECTIONS)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/infer_minor_corrections.py

.. autofunction:: pyTMD.infer_minor_corrections
