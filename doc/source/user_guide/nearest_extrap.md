nearest_extrap.py
=================

 - Uses kd-trees for nearest-neighbor extrapolation of valid model data

#### Calling Sequence
```python
from pyTMD.nearest_extrap import nearest_extrap
data = nearest_extrap(ilon,ilat,idata,lon,lat)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/nearest_extrap.py)

#### Arguments
 1. `ilon`: longitude of tidal model
 2. `ilat`: latitude of tidal model
 3. `idata`: tide model data
 4. `lon`: output longitude
 5. `lat`: output latitude

#### Keyword arguments
 - `fill_value`: invalid value
 - `dtype`: output data type
 - `cutoff`: return only neighbors within distance in kilometers
 - `EPSG`: projection of tide model data

#### Returns
 - `data`: interpolated data
