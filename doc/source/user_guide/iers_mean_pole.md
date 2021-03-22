iers_mean_pole.py
=================

 - Provides the angular coordinates of the IERS Conventional Mean Pole (CMP)
 - Based on the coordinates from IERS mean pole tabulations ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab
 - See ftp://tai.bipm.org/iers/convupdt/chapter7/

#### Calling Sequence
```python
from pyTMD.iers_mean_pole import iers_mean_pole
x,y,flag = iers_mean_pole(input_file,input_epoch,version,FILL_VALUE=FILL_VALUE)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/iers_mean_pole.py)

#### Arguments
 1. `input_file`: full path to mean-pole.tab file provided by IERS
 2. `input_epoch`: dates for which the angular coordinates of the Conventional Mean Pole are desired in decimal years
 3. `version`: Year of the conventional model (2003, 2010, 2015)

#### Keyword arguments
 - `FILL_VALUE`: value for invalid flags

#### Returns
 - `x`: Angular coordinate x of conventional mean pole [arcsec]
 - `y`: Angular coordinate y of conventional mean pole [arcsec]
 - `flag`: epoch is valid for version and version number is valid
