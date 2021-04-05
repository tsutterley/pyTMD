calc_delta_time.py
==================

 - Calculates the difference between dynamic time and universal time (TT - UT1) following Richard Ray's PERTH3 algorithms

#### Calling Sequence
```python
from pyTMD.calc_delta_time import calc_delta_time
deltat = calc_delta_time(delta_file, idays)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/calc_delta_time.py)

#### Arguments
 1. `delta_file`:
    - [http://maia.usno.navy.mil/ser7/deltat.data](http://maia.usno.navy.mil/ser7/deltat.data)
    - [ftp://cddis.nasa.gov/products/iers/deltat.data](ftp://cddis.nasa.gov/products/iers/deltat.data)
 2. `idays`: days relative to January 1, 1992 (MJD: 48622)

#### Returns
 - `deltat`: (TT - UT1) in days
