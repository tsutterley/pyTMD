convert_xy_ll.py
================

 - Converts lat/lon points to and from projected coordinates

#### Calling Sequence
```python
from pyTMD.convert_xy_ll import convert_xy_ll
x,y = convert_xy_ll(lon,lat,PROJ,'F')
lon,lat = convert_xy_ll(x,y,PROJ,'B')
```

#### Inputs
 1. `i1`: longitude ('F') or projection easting x ('B')
 2. `i2`: latitude ('F') or projection northing y ('B')
 3. `PROJ`: spatial reference system code for coordinate transformations
 4. `BF`: backwards ('B') or forward ('F') translations

#### Outputs
 - `o1`: projection easting x ('F') or longitude ('B')
 - `o2`: projection northing y ('F') or latitude ('B')
