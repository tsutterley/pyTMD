count_leap_seconds.py
=====================

 - Count number of leap seconds that have passed for each GPS time
 - Can be used for converting from GPS or TAI times (TAI = GPS + 19)

#### Calling Sequence
```python
from pyTMD.count_leap_seconds import count_leap_seconds
n_leaps = count_leap_seconds(GPS_Time)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/pyTMD/count_leap_seconds.py)

#### Inputs
 1. `GPS_Time`: seconds since January 6, 1980 at 00:00:00

#### Outputs
 - `n_leaps`: number of elapsed leap seconds
