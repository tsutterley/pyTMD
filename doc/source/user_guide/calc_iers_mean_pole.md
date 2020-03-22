calc_iers_mean_pole.py
======================

 - Calculates the mean pole coordinates x and y are obtained by Gaussian-weighted average of the IERS C01 pole coordinates  
 - Follows ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.readme  
 - Coordinates (particularly at end of the time series) will change with more dates  

#### Inputs
 1. `input_file`: full path to eopc01.1900-now.dat file provided by IERS

#### Outputs
 - `T`: date [decimal-years]
 - `xm`: mean pole coordinate x [arcsec]
 - `ym`: mean pole coordinate y [arcsec]
