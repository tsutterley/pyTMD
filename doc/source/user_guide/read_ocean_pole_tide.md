read_ocean_pole_tide.py
=======================

 - Reads ocean pole load tide coefficients provided by IERS as computed by
 [Desai et al. (2002)](https://doi.org/10.1029/2001JC001224) and
 [Desai et al. (2015)](https://doi.org/10.1007/s00190-015-0848-7)

#### Calling Sequence
```
from gravity_toolkit.read_ocean_pole_tide import read_ocean_pole_tide
ht = read_ocean_pole_tide(input_file)
```

#### Inputs
  1. `input_file`: [IERS 0.5x0.5 map of ocean pole tide coefficients](ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/additional_info/opoleloadcoefcmcor.txt.gz)

#### Outputs
 - `ur`: ocean pole tide coefficients
 - `glon`: grid longitude
 - `glat`: grid latitude
