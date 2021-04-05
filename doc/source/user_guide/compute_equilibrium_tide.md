compute_equilibrium_tide.py
===========================

 - Calculates the long-period equilibrium ocean tides using fifteen spectral lines from Cartwright-Tayler-Edden tables
 - Can be used to calculate tidal corrections for imagery

#### Calling Sequence
```python
from pyTMD.compute_equilibrium_tide import compute_equilibrium_tide
lpet = compute_equilibrium_tide(time,lat)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute_equilibrium_tide.py)

#### Arguments
 1. `t`: days relative to January 1, 1992 (MJD: 48622)
 2. `lat`: latitudes in degrees

#### Returns
 - `lpet`: long-period equilibrium tide values
