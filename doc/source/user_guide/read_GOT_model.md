read_GOT_model.py
=================

 - Reads files for Richard Ray's Global Ocean Tide (GOT) models  
 - Spatially interpolates tidal constituents to input coordinates  

#### Calling Sequence
```python
from pyTMD.read_GOT_model import read_GOT_model
amp,ph = read_GOT_model(ilon,ilat,directory,model_files,METHOD='spline')
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/pyTMD/read_GOT_model.py)

#### Inputs
  1. `ilon`: longitude to interpolate
  2. `ilat`: latitude to interpolate
  3. `directory`: data directory for tide data files
  4. `model_files`: list of model files for each constituent

#### Options
 - `METHOD`: interpolation method
    * `'bilinear'`: quick bilinear interpolation
    * `'spline'`: scipy bivariate spline interpolation
    * `'linear'`, `'nearest'`: scipy regular grid interpolations
 - `GZIP`: input files are compressed
 - `SCALE`: scaling factor for converting to output units

#### Outputs
- `amplitude`: amplitudes of tidal constituents
- `phase`: phases of tidal constituents
