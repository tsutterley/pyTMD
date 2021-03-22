bilinear_interp.py
==================

 - Bilinear interpolation of input data to output coordinates

#### Calling Sequence
```python
from pyTMD.bilinear_interp import bilinear_interp
data = bilinear_interp(ilon,ilat,idata,lon,lat)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/bilinear_interp.py)

#### Arguments
 1. `ilon`: longitude of tidal model
 2. `ilat`: latitude of tidal model
 3. `idata`: tide model data
 4. `lon`: output longitude
 5. `lat`: output latitude

#### Keyword arguments
 - `fill_value`: invalid value
 - `dtype`: output data type

#### Returns
 - `data`: interpolated data
