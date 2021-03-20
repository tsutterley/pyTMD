read_tide_model.py
==================

 - Reads OTIS format tidal solutions provided by Ohio State University and ESR
      * multi-constituent binary
      * ATLAS-compact binary
      * single-constituent binary
 - Spatially interpolates tidal constituents to input coordinates

#### Calling Sequence
```python
from pyTMD.read_tide_model import read_tide_model
amp,ph,D,c = read_tide_model(ilon, ilat, grid_file, model_file, EPSG,
    TYPE='z', METHOD='spline', GRID='OTIS')
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/read_tide_model.py)

#### Arguments
 1. `ilon`: longitude to interpolate
 2. `ilat`: latitude to interpolate
 3. `grid_file`: grid file for model
 4. `model_file`: model file(s) containing constituent data
 5. `EPSG`: projection of tide model data

#### Keyword arguments
 - `TYPE`: tidal variable to read
    * `'z'`: heights
    * `'u'`: horizontal transport velocities
    * `'U'`: horizontal depth-averaged transport
    * `'v'`: vertical transport velocities
    * `'v'`: vertical depth-averaged transport
 - `METHOD`: interpolation method
    * `'bilinear'`: quick bilinear interpolation
    * `'spline'`: scipy bivariate spline interpolation
    * `'linear'`, `'nearest'`: scipy regular grid interpolations
 - `GRID`: binary file type to read
    * `'ATLAS'`: reading a global solution with high-resolution local solutions
    * `'OTIS'`: combined global solution

#### Returns
 - `amplitude`: amplitudes of tidal constituents
 - `phase`: phases of tidal constituents
 - `D`: bathymetry of tide model
 - `constituents`: list of model constituents
