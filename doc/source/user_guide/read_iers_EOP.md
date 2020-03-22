read_iers_EOP.py
================

 - Provides the [daily earth orientation parameters (EOP) from IERS](http://www.usno.navy.mil/USNO/earth-orientation/eo-products/weekly)  
 - Data format: http://maia.usno.navy.mil/ser7/readme.finals  

#### Calling Sequence
```
from gravity_toolkit.read_iers_EOP import read_iers_EOP
MJD,x,y,flag = read_iers_EOP(input_file)
```

#### Inputs
 1. `input_file`:  full path to IERS EOP "finals" file

#### Outputs
 - `MJD`: modified Julian date of EOP measurements
 - `x`: Angular coordinate [arcsec]
 - `y`: Angular coordinate [arcsec]
 - `flag`: IERS (I) or Prediction (P) flag for polar motion values
