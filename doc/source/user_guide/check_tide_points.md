check_tide_points.py
====================

 - Check if points are within a tide model domain
 - Can check OTIS format tidal solutions provided by Ohio State University and ESR
 - Can check Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
 - Can check Finite Element Solution (FES) models provided by AVISO

#### Calling Sequence
```python
from pyTMD.check_tide_points import check_tide_points
valid = check_tide_points(x, y, DIRECTORY=DIRECTORY,
    MODEL=MODEL, EPSG=3031)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/check_tide_points.py)

#### Arguments
 1. `x`: x-coordinates in projection EPSG
 2. `y`: y-coordinates in projection EPSG

#### Keyword arguments
 - `DIRECTORY`: working data directory for tide models
 - `MODEL`: Tide model to use in correction
 - `EPSG`: input coordinate system
     * default: `3031` Polar Stereographic South, WGS84
 - `METHOD`: interpolation method
     * `'bilinear'`: quick bilinear interpolation
     * `'spline'`: scipy bivariate spline interpolation (default)
     * `'linear'`, `'nearest'`: scipy regular grid interpolations

#### Returns
 - `valid`: array describing if input coordinate is within model domain
