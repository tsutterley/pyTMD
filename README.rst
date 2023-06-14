=====
pyTMD
=====

|Language|
|License|
|PyPI Version|
|Anaconda-Server|
|Documentation Status|
|codecov|
|zenodo|

.. |Language| image:: https://img.shields.io/pypi/pyversions/pyTMD?color=green
   :target: https://www.python.org/

.. |License| image:: https://img.shields.io/github/license/tsutterley/pyTMD
   :target: https://github.com/tsutterley/pyTMD/blob/main/LICENSE

.. |PyPI Version| image:: https://img.shields.io/pypi/v/pyTMD.svg
   :target: https://pypi.python.org/pypi/pyTMD/

.. |Anaconda-Server| image:: https://img.shields.io/conda/vn/conda-forge/pytmd
   :target: https://anaconda.org/conda-forge/pytmd

.. |Documentation Status| image:: https://readthedocs.org/projects/pytmd/badge/?version=latest
   :target: https://pytmd.readthedocs.io/en/latest/?badge=latest

.. |codecov| image:: https://codecov.io/gh/tsutterley/pyTMD/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/tsutterley/pyTMD

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5555395.svg
   :target: https://doi.org/10.5281/zenodo.5555395

Python-based tidal prediction software for estimating ocean, load, solid Earth and pole tides

Ocean and load tidal predictions using OTIS, GOT and FES formatted tidal solutions

- `OSU Tidal Prediction Software (OTPS) <https://www.tpxo.net/otps>`_
- `ESR Tide Model Driver (TMD) Matlab Toolbox <https://www.esr.org/research/polar-tide-models/tmd-software/>`_
- `OSU Global and Regional Tide Models <https://www.tpxo.net>`_
- `ESR Polar Tide Models <https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/>`_
- `A Global Ocean Tide Model From TOPEX/POSEIDON Altimetry: GOT99.2 <https://ntrs.nasa.gov/citations/19990089548>`_
- `Finite Element Solution (FES) tide models <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html>`_
- `Delta times from US Naval Observatory (USNO) Earth Orientation Products <http://maia.usno.navy.mil/ser7/deltat.data>`_
- `Delta times from NASA Crustal Dynamics Data Information System (CDDIS) <ftp://cddis.nasa.gov/products/iers/deltat.data>`_

Radial solid Earth and pole tide displacements following IERS conventions

- `IERS Conventions (2010) <http://iers-conventions.obspm.fr/>`_
- `IERS Mean Pole Location <https://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab>`_
- `IERS Pole Coordinates to Calculate Mean Pole <https://hpiers.obspm.fr/iers/eop/eopc01/eopc01.1900-now.dat>`_
- `IERS Daily Earth Orientation Parameters (EOP) from USNO <http://www.usno.navy.mil/USNO/earth-orientation/eo-products/weekly>`_
- `IERS Daily Earth Orientation Parameters (EOP) from NASA CDDIS <ftp://cddis.nasa.gov/products/iers/finals.all>`_
- `IERS Ocean Pole Load Tide Coefficients Map <http://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/additional_info/opoleloadcoefcmcor.txt.gz>`_

Dependencies
############

- `dateutil: powerful extensions to datetime <https://dateutil.readthedocs.io/en/stable/>`_
- `lxml: processing XML and HTML in Python <https://pypi.python.org/pypi/lxml>`_
- `netCDF4: Python interface to the netCDF C library <https://unidata.github.io/netcdf4-python/>`_
- `numpy: Scientific Computing Tools For Python <https://www.numpy.org>`_
- `pyproj: Python interface to PROJ library <https://pypi.org/project/pyproj/>`_
- `scipy: Scientific Tools for Python <https://www.scipy.org/>`_
- `setuptools_scm: manager for python package versions using scm metadata <https://pypi.org/project/setuptools-scm>`_

References
##########

    T. C. Sutterley, T. Markus, T. A. Neumann, M. R. van den Broeke, J. M. van Wessem, and S. R. M. Ligtenberg,
    "Antarctic ice shelf thickness change from multimission lidar mapping", *The Cryosphere*,
    13, 1801-1817, (2019). `doi: 10.5194/tc-13-1801-2019 <https://doi.org/10.5194/tc-13-1801-2019>`_

    L. Padman, M. R. Siegfried, H. A. Fricker,
    "Ocean Tide Influences on the Antarctic and Greenland Ice Sheets", *Reviews of Geophysics*,
    56, 142-184, (2018). `doi: 10.1002/2016RG000546 <https://doi.org/10.1002/2016RG000546>`_

Download
########

| The program homepage is:
| https://github.com/tsutterley/pyTMD
| A zip archive of the latest version is available directly at:
| https://github.com/tsutterley/pyTMD/archive/main.zip

Software
########

| Matlab Tide Model Driver from Earth & Space Research is available at:
| https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5
| Fortran OSU Tidal Prediction Software OTPS is available at:
| https://www.tpxo.net/otps
| Incorporated into the NASA Cryosphere Altimetry Processing Toolkit at:
| https://github.com/fspaolo/captoolkit

Disclaimer
##########

This package includes software developed at NASA Goddard Space Flight Center (GSFC) and the University of Washington Applied Physics Laboratory (UW-APL).
It is not sponsored or maintained by the Universities Space Research Association (USRA), AVISO or NASA.
The software is provided here for your convenience but *with no guarantees whatsoever*.
It should not be used for coastal navigation or any application that may risk life or property.

Credits
#######

This project contains work and contributions from the `scientific community <./CONTRIBUTORS.rst>`_.
The Tidal Model Driver (TMD) Matlab Toolbox was developed by Laurie Padman, Lana Erofeeva and Susan Howard.
The OSU Tidal Inversion Software (OTIS) and OSU Tidal Prediction Software (OTPS) were developed by Lana Erofeeva and Gary Egbert (`copyright OSU <http://volkov.oce.orst.edu/tides/COPYRIGHT.pdf>`_, licensed for non-commercial use).
The NASA Goddard Space Flight Center (GSFC) PREdict Tidal Heights (PERTH3) software was developed by Richard Ray and Remko Scharroo.

License
#######

The content of this project is licensed under the `Creative Commons Attribution 4.0 Attribution license <https://creativecommons.org/licenses/by/4.0/>`_ and the source code is licensed under the `MIT license <LICENSE>`_.
