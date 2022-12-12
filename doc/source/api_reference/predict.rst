=======
predict
=======

- Predict tidal values using harmonic constants

  * At a single time (``map``) such as for imagery
  * Multiple times and locations (``drift``) such as for airborne and satellite altimetry
  * Time series at a location (``time_series``) such as to compare with tide gauges

- Predicts tidal values from minor constituents inferred using major constituents
- Predicts long-period equilibrium ocean tides


Calling Sequence
----------------

.. code-block:: python

    import pyTMD.predict
    ht = pyTMD.predict.time_series(time, hc, con)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/predict.py

.. autofunction:: pyTMD.predict.map

.. autofunction:: pyTMD.predict.drift

.. autofunction:: pyTMD.predict.time_series

.. autofunction:: pyTMD.predict.infer_minor

.. autofunction:: pyTMD.predict.equilibrium_tide
