nearest_extrap.py
=================

- Uses [kd-trees](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html) for nearest-neighbor (NN) extrapolation of valid tide model data

#### Calling Sequence
```python
from pyTMD.nearest_extrap import nearest_extrap
extrap = nearest_extrap(x,y,data,XI,YI)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/nearest_extrap.py)

#### Arguments
1. `x`: x-coordinates of tidal model
2. `y`: y: y-coordinates of tidal model
3. `data`: tide model data
4. `XI`: output x-coordinates
5. `YI`: output y-coordinates

#### Keyword arguments
- `fill_value`: invalid value
- `dtype`: output data type
- `cutoff`: return only neighbors within distance in kilometers
    * set to `np.inf` to extrapolate for all points
- `EPSG`: projection of tide model data

#### Returns
- `DATA`: extrapolated data
