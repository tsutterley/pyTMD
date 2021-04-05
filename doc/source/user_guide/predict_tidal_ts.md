predict_tidal_ts.py
===================

 - Predict tidal time series at a location using harmonic constants
 - Can be used to calculate errors as compared to tide gauges or to predict tides at a point

#### Calling Sequence
```python
from pyTMD.predict_tidal_ts import predict_tidal_ts
ht = predict_tidal_ts(time,hc,con)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/predict_tidal_ts.py)

#### Arguments
 1. `time`: days relative to January 1, 1992 (MJD: 48622)
 2. `hc`: harmonic constant vector (complex)
 3. `constituents`: tidal constituent IDs

#### Keyword arguments
 - `DELTAT`: time correction for converting to Ephemeris Time (days)
 - `CORRECTIONS`: use nodal corrections from OTIS/ATLAS or GOT models

#### Returns
 - `ht`: time series reconstructed using the nodal corrections
