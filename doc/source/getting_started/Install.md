Setup and Installation
======================

#### Dependencies
pyTMD is dependent on open source programs:
- [GDAL](https://gdal.org/index.html)
- [GEOS](https://trac.osgeo.org/geos)
- [PROJ](https://proj.org/)
- [HDF5](https://www.hdfgroup.org)
- [netCDF](https://www.unidata.ucar.edu/software/netcdf)
- [libxml2](http://xmlsoft.org/)
- [libxslt](http://xmlsoft.org/XSLT/)

#### Installation
Presently pyTMD is available for use as a [GitHub repository](https://github.com/tsutterley/pyTMD) and from the [Python Package Index (pypi)](https://pypi.org/project/pyTMD/).
The contents of the repository can be download from GitHub as a [zipped file](https://github.com/tsutterley/pyTMD/archive/main.zip) or cloned.
To use this repository from GitHub, please fork into your own account and then clone onto your system.  
```bash
git clone https://github.com/tsutterley/pyTMD.git
```
Can then install using `setuptools`
```bash
python setup.py install
```
or `pip`
```bash
python3 -m pip install --user .
```

To install from pypi using pip
```bash
python3 -m pip install --user pyTMD
```

Executable versions of this repository can also be tested using [Binder](https://mybinder.org/v2/gh/tsutterley/pyTMD/main) and [Pangeo](https://binder.pangeo.io/v2/gh/tsutterley/pyTMD/main).
