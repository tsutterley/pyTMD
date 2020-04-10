read_tide_model.py
==================

 - Reads OTIS format tidal solutions provided by Ohio State University and ESR
 - Spatially interpolates tidal constituents to input coordinates  

#### Calling Sequence
```python
from pyTMD.read_tide_model import read_tide_model
amp,ph,D,c = read_tide_model(ilon, ilat, grid_file, model_file, EPSG, type,
    METHOD='spline', GRID='OTIS')
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/pyTMD/read_tide_model.py)

#### Inputs
 1. `ilon`: longitude to interpolate
 2. `ilat`: latitude to interpolate
 3. `grid_file`: grid file for model
 4. `model_file`: model file containing each constituent
 5. `EPSG`: projection of tide model data
 6. `type`: tidal variable to run
    - 'z': heights
    - 'u': horizontal transport velocities
    - 'v': vertical transport velocities

#### Options
 - `METHOD`: interpolation method
    * `bilinear`: quick bilinear interpolation
    * `spline`: scipy bivariate spline interpolation
    * `linear`, `cubic`, `nearest`: scipy griddata interpolations
 - `GRID`: binary file type to read
    - 'ATLAS': reading a global solution with high-resolution local solutions
    - 'OTIS': combined global solution

#### Outputs
 - `amplitude`: amplitudes of tidal constituents
 - `phase`: phases of tidal constituents
 - `D`: bathymetry of tide model
 - `constituents`: list of model constituents
