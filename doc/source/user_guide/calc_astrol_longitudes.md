calc_astrol_longitudes.py
=========================

 - Computes the basic astronomical mean longitudes: s, h, p, N and PP
 - Note N is not N', i.e. N is decreasing with time.
 - Formulae for the period 1990-2010 were derived by David Cartwright

#### Calling Sequence
```python
from pyTMD.calc_astrol_longitudes import calc_astrol_longitudes
s,h,p,N,PP = calc_astrol_longitudes(MJD, ASTRO5=True)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/calc_astrol_longitudes.py)

#### Arguments
 1. `MJD`: Modified Julian Day of input date

#### Keyword arguments
 - `MEEUS`: use additional coefficients from Meeus Astronomical Algorithms
 - `ASTRO5`: use Meeus Astronomical coefficients as implemented in ASTRO5

#### Returns
 - `s`: mean longitude of moon (degrees)
 - `h`: mean longitude of sun (degrees)
 - `p`: mean longitude of lunar perigee (degrees)
 - `N`: mean longitude of ascending lunar node (degrees)
 - `PP`: longitude of solar perigee (degrees)
