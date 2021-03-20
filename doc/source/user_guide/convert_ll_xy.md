convert_ll_xy.py
================

 - Uses pyproj to convert lat/lon points to and from projected coordinates

#### Calling Sequence
```python
from pyTMD.convert_ll_xy import convert_ll_xy
x,y = convert_ll_xy(lon,lat,PROJ,'F')
lon,lat = convert_ll_xy(x,y,PROJ,'B')
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/convert_ll_xy.py)

#### Arguments
 1. `i1`: longitude ('F') or projection easting x ('B')
 2. `i2`: latitude ('F') or projection northing y ('B')
 3. `PROJ`: spatial reference system code for coordinate transformations
 4. `BF`: backwards ('B') or forward ('F') translations

#### Keyword arguments
 - `EPSG`: spatial reference system code for input (F) and output (B) coordinates

#### Returns
 - `o1`: projection easting x ('F') or longitude ('B')
 - `o2`: projection northing y ('F') or latitude ('B')
