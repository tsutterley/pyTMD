load_nodal_corrections.py
=========================

 - Calculates the nodal corrections for tidal constituents
 - Based on Richard Ray's ARGUMENTS fortran subroutine

#### Calling Sequence
```python
from pyTMD.load_nodal_corrections import load_nodal_corrections
pu,pf,G = load_nodal_corrections(MJD,constituents)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/load_nodal_corrections.py)

#### Arguments
 1. `MJD`: Modified Julian Day of input date
 2. `zmajor`: Complex oscillations for given constituents/points
 3. `constituents`: tidal constituent IDs

#### Returns
 -  `pu`,`pf`: nodal corrections for the constituents
 - `G`: phase correction in degrees
