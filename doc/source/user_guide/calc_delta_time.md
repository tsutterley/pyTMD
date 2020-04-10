calc_delta_time.py
==================

 - Calculates the difference between universal time and dynamical time (TT - UT1) following Richard Ray's PERTH3 algorithms

#### Calling Sequence
```python
from pyTMD.calc_delta_time import calc_delta_time
deltat = calc_delta_time(delta_file, iMJD)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/pyTMD/calc_delta_time.py)

#### Inputs
 1. `delta_file`:  
    - http://maia.usno.navy.mil/ser7/deltat.data  
    - ftp://cddis.nasa.gov/products/iers/deltat.data  
 2. `iMJD`: Modified Julian Day of times to interpolate  

#### Outputs
 - `deltat`: (TT - UT1) in days
