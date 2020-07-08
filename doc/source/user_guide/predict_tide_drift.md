predict_tide_drift.py
===============

 - Predict tides at multiple times and locations using harmonic constants  
 - Can be used to calculate tidal corrections for airborne and satellite altimetry  

#### Calling Sequence
```python
from pyTMD.predict_tide_drift import predict_tide_drift
ht = predict_tide_drift(time,hc,con)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/pyTMD/predict_tide_drift.py)

#### Inputs
 1. `time`: days relative to Jan 1, 1992 (48622mjd)
 2. `hc`: harmonic constant vector (complex)
 3. `constituents`: tidal constituent IDs

#### Options
 - `DELTAT`: time correction for converting to Ephemeris Time (days)
 - `CORRECTIONS`: use nodal corrections from OTIS/ATLAS or GOT models

#### Outputs
 - `ht`: time series reconstructed using the nodal corrections
