load_constituent.py
====================

 - Loads parameters for a given tidal constituent

#### Calling Sequence
```python
from pyTMD.load_constituent import load_constituent
amplitude,phase,omega,alpha,species = load_constituent(c)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/load_constituent.py)

#### Arguments
 1. `c`: tidal constituent IDs

#### Returns
 - `amplitude`: amplitude of equilibrium tide in m for tidal constituent
 - `phase`: phase of tidal constituent
 - `omega`: angular frequency of constituent in radians
 - `alpha`: load love number of tidal constituent
 - `species`: spherical harmonic dependence of quadrupole potential
