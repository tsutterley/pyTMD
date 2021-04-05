read_ocean_pole_tide.py
=======================

 - Reads ocean pole load tide coefficients provided by IERS as computed by
 [Desai et al. (2002)](https://doi.org/10.1029/2001JC001224) and
 [Desai et al. (2015)](https://doi.org/10.1007/s00190-015-0848-7)

#### Calling Sequence
```python
from pyTMD.utilities import get_data_path
from pyTMD.read_ocean_pole_tide import read_ocean_pole_tide
ocean_pole_tide_file = get_data_path(['data','opoleloadcoefcmcor.txt.gz'])
ur,un,ue,glon,glat = read_ocean_pole_tide(ocean_pole_tide_file)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/read_ocean_pole_tide.py)

#### Arguments
  1. `input_file`: [IERS map of ocean pole tide coefficients](ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/additional_info/opoleloadcoefcmcor.txt.gz)

#### Returns
 - `ur`: radial ocean pole tide coefficients
 - `un`: north ocean pole tide coefficients
 - `ue`: east ocean pole tide coefficients
 - `glon`: grid longitude
 - `glat`: grid latitude
