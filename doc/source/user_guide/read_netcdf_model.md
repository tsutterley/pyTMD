read_netcdf_model.py
====================

 - Reads netCDF format tidal solutions provided by Ohio State University and ESR
 - Spatially interpolates tidal constituents to input coordinates

#### Calling Sequence
```python
from pyTMD.read_netcdf_model import read_netcdf_model
amp,ph,D,c = read_netcdf_model(ilon,ilat,grid_file,model_files,TYPE='z',METHOD='spline')
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/read_netcdf_model.py)

#### Arguments
  1. `ilon`: longitude to interpolate
  2. `ilat`: latitude to interpolate
  3. `grid_file`: grid file for model
  4. `model_files`: list of model files for each constituent

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
 - `GZIP`: input netCDF4 files are compressed
 - `SCALE`: scaling factor for converting to output units

#### Returns
 - `amplitude`: amplitudes of tidal constituents
 - `phase`: phases of tidal constituents
 - `D`: bathymetry of tide model
 - `constituents`: list of model constituents
