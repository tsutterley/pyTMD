predict_tidal_ts.py
===================

 - Predict tidal time series at a location using harmonic constants
 - Can be used to calculate errors as compared to tide gauges or to predict tides at a point  

#### Calling Sequence
```python
from pyTMD.predict_tidal_ts import predict_tidal_ts
ht = predict_tidal_ts(time,hc,con)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/pyTMD/predict_tidal_ts.py)

#### Inputs
 1. `time`: days relative to Jan 1, 1992 (48622mjd)
 2. `hc`: harmonic constant vector (complex)
 3. `constituents`: tidal constituent IDs

#### Options
 - `DELTAT`: time correction for converting to Ephemeris Time (days)
 - `CORRECTIONS`: use nodal corrections from OTIS/ATLAS or GOT models

#### Outputs
 - `ht`: time series reconstructed using the nodal corrections
