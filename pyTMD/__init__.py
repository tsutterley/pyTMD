"""
A tide prediction toolkit for Python
====================================

pyTMD contains Python tools for reading OTIS, GOT and FES formatted tidal
solutions to predict ocean and load tides

The package works using scientific Python packages (numpy, scipy and pyproj)
combined with data storage in netCDF4 and HDF5 and mapping using
matplotlib and cartopy

Documentation is available at https://pytmd.readthedocs.io
"""
import pyTMD.astro
import pyTMD.eop
import pyTMD.interpolate
import pyTMD.predict
import pyTMD.spatial
import pyTMD.time
import pyTMD.tools
import pyTMD.utilities
import pyTMD.version
from pyTMD import io
from pyTMD.arguments import arguments
from pyTMD.check_tide_points import check_tide_points
from pyTMD.compute_tide_corrections import (
    compute_corrections,
    compute_tide_corrections,
    compute_LPET_corrections,
    compute_LPT_corrections,
    compute_OPT_corrections,
    compute_SET_corrections,
)
from pyTMD.constants import (
    constants,
    _ellipsoids
)
from pyTMD.convert_crs import convert_crs
from pyTMD.load_constituent import load_constituent
from pyTMD.tidal_ellipse import tidal_ellipse

# Deprecated functions
from pyTMD.calc_astrol_longitudes import calc_astrol_longitudes
from pyTMD.convert_ll_xy import convert_ll_xy
from pyTMD.load_nodal_corrections import load_nodal_corrections

# get semantic version from setuptools-scm
__version__ = pyTMD.version.version
