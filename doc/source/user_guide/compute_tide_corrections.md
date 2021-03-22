compute_tide_corrections.py
===========================

 - Calculates tidal elevations for correcting elevation or imagery data
 - Can use OTIS format tidal solutions provided by Ohio State University and ESR
 - Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
 - Can use Finite Element Solution (FES) models provided by AVISO

#### Calling Sequence
```python
from pyTMD.compute_tide_corrections import compute_tide_corrections
tide = compute_tide_corrections(x, y, delta_time, DIRECTORY=DIRECTORY,
    MODEL=MODEL, EPOCH=(2000,1,1,0,0,0), EPSG=3031, TYPE='drift')
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute_tide_corrections.py)

#### Arguments
 1. `x`: x-coordinates in projection EPSG
 2. `y`: y-coordinates in projection EPSG
 3. `delta_time`: seconds since EPOCH

#### Keyword arguments
 - `DIRECTORY`: working data directory for tide models
 - `MODEL`: Tide model to use in correction
 - `EPOCH`: time period for calculating delta times
     * default: J2000 (seconds since 2000-01-01T00:00:00)
 - `TYPE`: input data type
     * 'drift': drift buoys or satellite/airborne altimetry (time per data point)
     * 'grid': spatial grids or images (single time per image)
 - `TIME`: time type if need to compute leap seconds to convert to UTC
     * 'GPS': leap seconds needed
     * 'TAI': leap seconds needed (TAI = GPS + 19 seconds)
     * 'UTC': no leap seconds needed
 - `EPSG`: input coordinate system
     * default: 3031 Polar Stereographic South, WGS84
 - `METHOD`: interpolation method
     * `bilinear`: quick bilinear interpolation
     * `spline`: scipy bivariate spline interpolation (default)
     * `linear`, `nearest`: scipy regular grid interpolations
 - `EXTRAPOLATE`: Extrapolate with nearest-neighbors
 - `FILL_VALUE`: output invalid value

#### Returns
 - `tide`: tide height correction reconstructed using the nodal corrections
