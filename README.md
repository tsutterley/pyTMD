pyTMD
=====

[![Language](https://img.shields.io/badge/python-v3.7-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/tsutterley/pyTMD/blob/master/LICENSE)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tsutterley/pyTMD/master)
[![Binder](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/tsutterley/pyTMD/master)

#### Python-based tidal prediction software that reads OTIS and GOT formatted tidal solutions for calculating ocean and load tides

- [OSU Tidal Prediction Software (OTPS)](http://volkov.oce.orst.edu/tides/otps.html)  
- [ESR Tide Model Driver (TMD) Matlab Toolbox](https://www.esr.org/research/polar-tide-models/tmd-software/)  
- [OSU List of Regional Tide Models](http://volkov.oce.orst.edu/tides/region.html)  
- [ESR List of Polar Tide Models](https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/)  
- [A Global Ocean Tide Model From TOPEX/POSEIDON Altimetry: GOT99.2](https://ntrs.nasa.gov/search.jsp?R=19990089548)  
- [Delta times from US Naval Observatory (USNO) Earth Orientation Products](http://maia.usno.navy.mil/ser7/deltat.data)  
- [Delta times from NASA Crustal Dynamics Data Information System (CDDIS)](ftp://cddis.nasa.gov/products/iers/deltat.data)  

#### Pole tide prediction software for calculating radial pole tide displacements

- [IERS Conventions (2010)](http://iers-conventions.obspm.fr/)  
- [IERS Mean Pole Location](https://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab)  
- [IERS Pole Coordinates to Calculate Mean Pole](https://hpiers.obspm.fr/iers/eop/eopc01/eopc01.1900-now.dat)  
- [IERS Daily Earth Orientation Parameters (EOP) from USNO](http://www.usno.navy.mil/USNO/earth-orientation/eo-products/weekly)  
- [IERS Daily Earth Orientation Parameters (EOP) from NASA CDDIS](ftp://cddis.nasa.gov/products/iers/finals.all)  
- [IERS Ocean Pole Load Tide Coefficients Map](http://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/additional_info/opoleloadcoefcmcor.txt.gz)

#### Dependencies
[numpy: Scientific Computing Tools For Python](https://www.numpy.org)  
[scipy: Scientific Tools for Python](https://www.scipy.org/)  
[pyproj: Python interface to PROJ library](https://pypi.org/project/pyproj/)  
[h5py: Python interface for Hierarchal Data Format 5 (HDF5)](https://www.h5py.org/)  
[netCDF4: Python interface to the netCDF C library](https://unidata.github.io/netcdf4-python/)  
[matplotlib: Python 2D plotting library](http://matplotlib.org/)  

#### Reference
T. C. Sutterley, T. Markus, T. A. Neumann, M. R. van den Broeke, J. M. van Wessem, and S. R. M. Ligtenberg, S. R. M.,
"Antarctic ice shelf thickness change from multimission lidar mapping", *The Cryosphere*,
13, 1801–1817, (2019). [doi:tc-13-1801-2019](https://doi.org/10.5194/tc-13-1801-2019)  

#### Download
The program homepage is:  
https://github.com/tsutterley/pyTMD  
A zip archive of the latest version is available directly at:  
https://github.com/tsutterley/pyTMD/archive/master.zip  
Incorporated into the NASA Cryosphere Altimetry Processing Toolkit at:  
https://github.com/fspaolo/captoolkit  

#### Disclaimer  
This program is not sponsored or maintained by the Universities Space Research Association (USRA) or NASA.  It is provided here for your convenience but _with no guarantees whatsoever_.  

#### Credits
The Tidal Model Driver (TMD) Matlab Toolbox was developed by Laurie Padman and Lana Erofeeva.  The OSU Tidal Inversion Software (OTIS) and OSU Tidal Prediction Software (OTPS) were developed by Lana Erofeeva and Gary Egbert ([copyright OSU](http://volkov.oce.orst.edu/tides/COPYRIGHT.pdf), licensed for non-commercial use).  The NASA Goddard Space Flight Center (GSFC) PREdict Tidal Heights (PERTH3) software was developed by Richard Ray and Remko Scharroo.  
