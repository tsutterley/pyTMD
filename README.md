pyTMD
=====

[![Language](https://img.shields.io/badge/python-v3.7-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/tsutterley/pyTMD/blob/main/LICENSE)
[![PyPI Version](https://img.shields.io/pypi/v/pyTMD.svg)](https://pypi.python.org/pypi/pyTMD/)
[![Documentation Status](https://readthedocs.org/projects/pytmd/badge/?version=latest)](https://pytmd.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/tsutterley/pyTMD/branch/main/graph/badge.svg)](https://codecov.io/gh/tsutterley/pyTMD)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tsutterley/pyTMD/main)
[![Binder](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/tsutterley/pyTMD/main)

#### Python-based tidal prediction software that reads OTIS, GOT and FES formatted tidal solutions for calculating ocean and load tides

- [OSU Tidal Prediction Software (OTPS)](https://www.tpxo.net/otps)  
- [ESR Tide Model Driver (TMD) Matlab Toolbox](https://www.esr.org/research/polar-tide-models/tmd-software/)  
- [OSU Global and Regional Tide Models](https://www.tpxo.net)  
- [ESR Polar Tide Models](https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/)  
- [A Global Ocean Tide Model From TOPEX/POSEIDON Altimetry: GOT99.2](https://ntrs.nasa.gov/search.jsp?R=19990089548)  
- [Finite Element Solution (FES) tide models](https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html)  
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
 - [numpy: Scientific Computing Tools For Python](https://www.numpy.org)  
 - [scipy: Scientific Tools for Python](https://www.scipy.org/)  
 - [pyproj: Python interface to PROJ library](https://pypi.org/project/pyproj/)  
 - [dateutil: powerful extensions to datetime](https://dateutil.readthedocs.io/en/stable/)  
 - [lxml: processing XML and HTML in Python](https://pypi.python.org/pypi/lxml)  
 - [PyYAML: YAML parser and emitter for Python](https://github.com/yaml/pyyaml)  
 - [gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)](https://pypi.python.org/pypi/GDAL)  
 - [h5py: Python interface for Hierarchal Data Format 5 (HDF5)](https://www.h5py.org/)  
 - [netCDF4: Python interface to the netCDF C library](https://unidata.github.io/netcdf4-python/)  
 - [matplotlib: Python 2D plotting library](https://matplotlib.org/)  
 - [cartopy: Python package designed for geospatial data processing](https://scitools.org.uk/cartopy/docs/latest/)  
 - [ipyleaflet: Jupyter / Leaflet bridge enabling interactive maps](https://github.com/jupyter-widgets/ipyleaflet)  
 - [read-ICESat-2: Python tools to read data from the NASA ICESat-2 mission](https://github.com/tsutterley/read-ICESat-2/)  
 - [read-ATM1b-QFIT-binary: Python reader for Airborne Topographic Mapper (ATM) QFIT data products](https://github.com/tsutterley/read-ATM1b-QFIT-binary)  

#### Reference
T. C. Sutterley, T. Markus, T. A. Neumann, M. R. van den Broeke, J. M. van Wessem, and S. R. M. Ligtenberg,
"Antarctic ice shelf thickness change from multimission lidar mapping", *The Cryosphere*,
13, 1801-1817, (2019). [doi:tc-13-1801-2019](https://doi.org/10.5194/tc-13-1801-2019)  

L. Padman, M. R. Siegfried, H. A. Fricker,
"Ocean Tide Influences on the Antarctic and Greenland Ice Sheets", *Reviews of Geophysics*,
56, 142-184, (2018). [doi:10.1002/2016RG000546](https://doi.org/10.1002/2016RG000546)  

#### Download
The program homepage is:  
https://github.com/tsutterley/pyTMD  
A zip archive of the latest version is available directly at:  
https://github.com/tsutterley/pyTMD/archive/main.zip  

#### Software
Matlab Tide Model Driver from Earth & Space Research is available at:  
https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5  
Fortran OSU Tidal Prediction Software OTPS is available at:  
https://www.tpxo.net/otps  
pyTMD was incorporated into the NASA Cryosphere Altimetry Processing Toolkit at:  
https://github.com/fspaolo/captoolkit  

#### Disclaimer  
This package includes software developed at NASA Goddard Space Flight Center (GSFC) and the University of Washington Applied Physics Laboratory (UW-APL).
It is not sponsored or maintained by the Universities Space Research Association (USRA), AVISO or NASA.
The software is provided here for your convenience but _with no guarantees whatsoever_.
It should not be used for coastal navigation or any application that may risk life or property.  

#### Credits
The Tidal Model Driver (TMD) Matlab Toolbox was developed by Laurie Padman, Lana Erofeeva and Susan Howard.
The OSU Tidal Inversion Software (OTIS) and OSU Tidal Prediction Software (OTPS) were developed by Lana Erofeeva and Gary Egbert ([copyright OSU](http://volkov.oce.orst.edu/tides/COPYRIGHT.pdf), licensed for non-commercial use).
The NASA Goddard Space Flight Center (GSFC) PREdict Tidal Heights (PERTH3) software was developed by Richard Ray and Remko Scharroo.  

#### License
The content of this project is licensed under the [Creative Commons Attribution 4.0 Attribution license](https://creativecommons.org/licenses/by/4.0/) and the source code is licensed under the [MIT license](LICENSE).  
