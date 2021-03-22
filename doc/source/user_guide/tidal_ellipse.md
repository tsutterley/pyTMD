tidal_ellipse.py
================

 - Expresses the amplitudes and phases for the u and v components in terms of four ellipse parameters using [Foreman's formula](https://www.sciencedirect.com/science/article/pii/0309170889900171)

#### Calling Sequence
```python
from pyTMD.tidal_ellipse import tidal_ellipse
umajor,uminor,uincl,uphase = tidal_ellipse(u,v)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/tidal_ellipse.py)

#### Arguments
 1. `u`: zonal current (EW)
 2. `v`: meridional current (NS)

#### Returns
 - `umajor`: amplitude of the semimajor semi-axis
 - `uminor`: amplitude of the semiminor semi-axis
 - `uincl`: angle of inclination of the northern semimajor semi-axis
 - `uphase`: phase lag of the maximum current behind the maximum tidal potential
     of the individual constituent
