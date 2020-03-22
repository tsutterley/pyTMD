load_constituent.py
====================

 - Loads parameters for a given tidal constituent

#### Calling Sequence
```
from gravity_toolkit.load_constituent import load_constituent
flag,amplitude,phase,alpha,species = load_constituent(c)
```

#### Inputs
 1. `c`: tidal constituent IDs

#### Outputs
 - `flag`: test for constituent being part of tidal program
 - `amplitude`, phase, frequency of tidal constituent
 - `alpha`: load love number of tidal constituent
 - `species`: spherical harmonic dependence of quadropole potential
