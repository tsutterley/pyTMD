compute_equilibrium_tide.py
===========================

 - Calculates the long-period equilibrium ocean tides
 - Can be used to calculate tidal corrections for imagery  

#### Calling Sequence
```python
from pyTMD.compute_equilibrium_tide import compute_equilibrium_tide
lpet = compute_equilibrium_tide(time,lat)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute_equilibrium_tide.py)

#### Inputs
 1. `t`: days relative to Jan 1, 1992 (48622mjd)
 2. `lat`: latitudes in degrees

#### Outputs
 - `lpet`: long-period equilibrium tide values
