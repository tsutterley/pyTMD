predict_tide.py
===============

 - Predict tides at a single time using harmonic constants
 - Can be used to calculate tidal corrections for imagery  

#### Calling Sequence
```python
from pyTMD.predict_tide import predict_tide
ht = predict_tide(time,hc,con)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/pyTMD/predict_tide.py)

#### Inputs
 1. `time`: days relative to Jan 1, 1992 (48622mjd)
 2. `hc`: harmonic constant vector (complex)
 3. `constituents`: tidal constituent IDs

#### Options
 - `DELTAT`: time correction for converting to Ephemeris Time (days)
 - `CORRECTIONS`: use nodal corrections from OTIS/ATLAS or GOT models

#### Outputs
 - `ht`: tide values reconstructed using the nodal corrections
