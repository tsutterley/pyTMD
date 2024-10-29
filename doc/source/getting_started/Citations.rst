====================
Citation Information
====================

References
##########

This work was initially supported by an appointment to the NASA Postdoctoral
Program (NPP) at NASA Goddard Space Flight Center (GSFC), administered by
Universities Space Research Association (USRA) under contract with NASA.
It is currently supported under the NASA Cryospheric Sciences Program (Grant Number 80NSSC22K0379).
The programs included in this software have contributed most recently to the
following work:

    T. C. Sutterley, T. Markus, T. A. Neumann, M. R. van den Broeke, J. M. van Wessem, and S. R. M. Ligtenberg,
    "Antarctic ice shelf thickness change from multimission lidar mapping", *The Cryosphere*,
    13, 1801--1817, (2019). `doi: 10.5194/tc-13-1801-2019 <https://doi.org/10.5194/tc-13-1801-2019>`_


If you have used ``pyTMD`` in your work, please consider citing our library:

    T. C. Sutterley, K. Alley, K. Brunt, S. Howard, L. Padman, and M. Siegfried,
    "pyTMD: Python-based tidal prediction software", (2017).
    `doi: 10.5281/zenodo.5555395 <https://doi.org/10.5281/zenodo.5555395>`_

Contributors
############

.. include:: ../../../CONTRIBUTORS.rst

Development
###########

``pyTMD`` is an open source project.
We welcome any help in maintaining and developing the software and documentation.
Anyone at any career stage and with any level of coding experience can contribute.
Please see the `Contribution Guidelines <./Contributing.html>`_ for more information.

Problem Reports
###############

If you have found a problem in ``pyTMD``, or you would like to suggest an improvement or modification,
please submit a `GitHub issue <https://github.com/tsutterley/pyTMD/issues>`_ and we will get back to you.

Dependencies
############

This software is also dependent on other commonly used Python packages:

- `dateutil: powerful extensions to datetime <https://dateutil.readthedocs.io/en/stable/>`_
- `lxml: processing XML and HTML in Python <https://pypi.python.org/pypi/lxml>`_
- `netCDF4: Python interface to the netCDF C library <https://unidata.github.io/netcdf4-python/>`_
- `numpy: Scientific Computing Tools For Python <https://www.numpy.org>`_
- `pyproj: Python interface to PROJ library <https://pypi.org/project/pyproj/>`_
- `scipy: Scientific Tools for Python <https://www.scipy.org/>`_
- `setuptools_scm: manager for python package versions using scm metadata <https://pypi.org/project/setuptools-scm>`_
- `timescale: Python tools for time and astronomical calculations <https://pypi.org/project/timescale/>`_

Optional Dependencies
---------------------

- `cartopy: Python package designed for geospatial data processing <https://scitools.org.uk/cartopy/docs/latest/>`_
- `gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL) <https://pypi.python.org/pypi/GDAL>`_
- `h5py: Python interface for Hierarchal Data Format 5 (HDF5) <https://www.h5py.org/>`_
- `ipyleaflet: Jupyter / Leaflet bridge enabling interactive maps <https://github.com/jupyter-widgets/ipyleaflet>`_
- `ipywidgets: interactive HTML widgets for Jupyter notebooks and IPython <https://ipywidgets.readthedocs.io/en/latest/>`_
- `jplephem: Python implementation of the math for predicting raw (x,y,z) planetary positions from JPL ephemerides <https://pypi.org/project/jplephem/>`_
- `matplotlib: Python 2D plotting library <https://matplotlib.org/>`_
- `pandas: Python Data Analysis Library <https://pandas.pydata.org/>`_
- `PyYAML: YAML parser and emitter for Python <https://github.com/yaml/pyyaml>`_

Credits
#######

The Tidal Model Driver (TMD) Matlab Toolbox was developed by Laurie Padman, Lana Erofeeva and Susan Howard.
An updated version of the TMD Matlab Toolbox (TMD3) was developed by Chad Greene.
The OSU Tidal Inversion Software (OTIS) and OSU Tidal Prediction Software (OTPS) were developed by
Lana Erofeeva and Gary Egbert (`copyright OSU <http://volkov.oce.orst.edu/tides/COPYRIGHT.pdf>`_,
licensed for non-commercial use).
The NASA Goddard Space Flight Center (GSFC) PREdict Tidal Heights (PERTH3) software was developed by
Richard Ray and Remko Scharroo.
An updated and more versatile version of the NASA GSFC tidal prediction software (PERTH5) was developed by Richard Ray.

Data Citations
##############

Internally, ``pyTMD`` includes datasets from the following:

    D. E. Cartwright and R. J. Tayler, "New Computations of the Tide-generating Potential,"
    *Geophysical Journal of the Royal Astronomical Society*, 23(1), 45--73. (1971).
    `doi: 10.1111/j.1365-246X.1971.tb01803.x <https://doi.org/10.1111/j.1365-246X.1971.tb01803.x>`_

    D. E. Cartwright and A. C. Edden, "Corrected Tables of Tidal Harmonics,"
    *Geophysical Journal of the Royal Astronomical Society*, 33(3), 253--264, (1973).
    `doi: 10.1111/j.1365-246X.1973.tb03420.x <https://doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_
    
    S. Desai, J. Wahr and B. Beckley "Revisiting the pole tide for and from satellite altimetry",
    *Journal of Geodesy*, 89(12), p1233-1243, (2015).
    `doi: 10.1007/s00190-015-0848-7 <https://doi.org/10.1007/s00190-015-0848-7>`_
    
    *IERS Conventions (2010)*, G. Petit and B. Luzum (eds.), (IERS Technical Note; 36),
    Frankfurt am Main: Verlag des Bundesamts f\ |uuml|\ r Kartographie und Geod\ |auml|\ sie, 179 pp. (2010).
    `ISBN 3-89888-989-6 <https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html>`_

Disclaimer
##########

This package includes software developed at NASA Goddard Space Flight Center (GSFC) and the University
of Washington Applied Physics Laboratory (UW-APL).
It is not sponsored or maintained by the Universities Space Research Association (USRA), AVISO or NASA.
Outputs from this software should be used for scientific or technical purposes only.
This software should not be used for coastal navigation or any application that may risk life or property.

.. |auml|    unicode:: U+00E4 .. LATIN SMALL LETTER A WITH DIAERESIS
.. |uuml|    unicode:: U+00FC .. LATIN SMALL LETTER U WITH DIAERESIS

