infer_minor_corrections.py
==========================

 - Calculates the tidal corrections for minor constituents
 - Based on Richard Ray's PERTH3 (PREdict Tidal Heights) algorithms

#### Calling Sequence
```python
from pyTMD.infer_minor_corrections import infer_minor_corrections
dh = infer_minor_corrections(t, zmajor, constituents,
    DELTAT=DELTAT, CORRECTIONS=CORRECTIONS)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/infer_minor_corrections.py)

#### Inputs
 1. `t`: days relative to Jan 1, 1992 (48622mjd)
 2. `zmajor`: Complex oscillations for given constituents/points
 3. `constituents`: tidal constituent IDs

#### Options
 - `DELTAT`: time correction for converting to Ephemeris Time (days)
 - `CORRECTIONS`: use nodal corrections from OTIS/ATLAS or GOT models

#### Outputs
 - `dh`: height from minor constituents
