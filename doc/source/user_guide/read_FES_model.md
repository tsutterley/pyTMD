read_FES_model.py
=================

 - Reads files for Finite Element Solution (FES) models
      * ascii
      * netCDF4
 - Spatially interpolates tidal constituents to input coordinates

#### Calling Sequence
```python
from pyTMD.read_FES_model import read_FES_model
amp,ph = read_FES_model(ilon, ilat, model_files, TYPE='z',
    VERSION=version,METHOD='spline',GZIP=True,SCALE=1.0/100.0)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/read_FES_model.py)

#### Arguments
  1. `ilon`: longitude to interpolate
  2. `ilat`: latitude to interpolate
  3. `model_files`: list of model files for each constituent

#### Keyword arguments
- `TYPE`: tidal variable to read
   * `'z'`: heights
   * `'u'`: horizontal transport velocities
   * `'v'`: vertical transport velocities
 - `VERSION`: tide model version to read
    * `'FES1999'`
    * `'FES2004'`
    * `'FES2012'`
    * `'FES2014'`
 - `METHOD`: interpolation method
    * `'bilinear'`: quick bilinear interpolation
    * `'spline'`: scipy bivariate spline interpolation
    * `'linear'`, `'nearest'`: scipy regular grid interpolations
 - `GZIP`: input files are compressed
 - `SCALE`: scaling factor for converting to output units

#### Returns
- `amplitude`: amplitudes of tidal constituents
- `phase`: phases of tidal constituents
