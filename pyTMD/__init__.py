"""
A tide prediction toolkit for Python
====================================

pyTMD is a Python-based tidal prediction software for estimating ocean,
load, solid Earth and pole tides.

The package works using scientific Python packages (numpy, scipy and pyproj)
combined with data storage in netCDF4 and HDF5 and mapping using
matplotlib and cartopy

Documentation is available at https://pytmd.readthedocs.io
"""
import pyTMD.arguments
import pyTMD.astro
import pyTMD.compute
import pyTMD.ellipse
import pyTMD.interpolate
import pyTMD.predict
import pyTMD.spatial
import pyTMD.tools
import pyTMD.utilities
import pyTMD.version
from pyTMD import io
from pyTMD import solve
from pyTMD.check_points import check_points
from pyTMD.crs import (
    crs,
    datum,
    _ellipsoids
)

# Deprecated functions
from pyTMD.compute_tide_corrections import (
    compute_corrections,
    compute_tide_corrections,
    compute_LPET_corrections,
    compute_LPT_corrections,
    compute_OPT_corrections,
    compute_SET_corrections,
)
import pyTMD.eop
import pyTMD.time

# get semantic version from setuptools-scm
__version__ = pyTMD.version.version
# read model database
models = io.load_database()
