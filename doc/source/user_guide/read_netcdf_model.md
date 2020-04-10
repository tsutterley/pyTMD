read_netcdf_model.py
====================

 - Reads netCDF format tidal solutions provided by Ohio State University and ESR
 - Spatially interpolates tidal constituents to input coordinates  

#### Calling Sequence
```python
from pyTMD.read_netcdf_model import read_netcdf_model
amp,ph,D,c = read_netcdf_model(ilon,ilat,directory,model_files,type,METHOD='spline')
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/pyTMD/read_netcdf_model.py)

#### Inputs
  1. `ilon`: longitude to interpolate
  2. `ilat`: latitude to interpolate
  3. `directory`: data directory for tide data files
  4. `grid_file`: grid file for model
  5. `model_files`: list of model files for each constituent
  6. `type`: tidal variable to run
     - 'z': heights
     - 'u': horizontal transport velocities
     - 'v': vertical transport velocities

#### Options
 - `METHOD`: interpolation method
    * `bilinear`: quick bilinear interpolation
    * `spline`: scipy bivariate spline interpolation
    * `linear`, `cubic`, `nearest`: scipy griddata interpolations
 - `GZIP`: input netCDF4 files are compressed
 - `SCALE`: scaling factor for converting to output units

#### Outputs
 - `amplitude`: amplitudes of tidal constituents
 - `phase`: phases of tidal constituents
 - `D`: bathymetry of tide model
 - `constituents`: list of model constituents
