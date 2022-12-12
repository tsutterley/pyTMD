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
import pyTMD.eop
import pyTMD.predict
import pyTMD.time
import pyTMD.tools
import pyTMD.utilities
import pyTMD.version
from pyTMD import io
from pyTMD.bilinear_interp import bilinear_interp
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.calc_astrol_longitudes import calc_astrol_longitudes
from pyTMD.check_tide_points import check_tide_points
from pyTMD.compute_equilibrium_tide import compute_equilibrium_tide
from pyTMD.compute_tide_corrections import compute_tide_corrections
from pyTMD.convert_ll_xy import convert_ll_xy
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.load_constituent import load_constituent
from pyTMD.load_nodal_corrections import load_nodal_corrections
from pyTMD.model import model
from pyTMD.nearest_extrap import nearest_extrap
from pyTMD.tidal_ellipse import tidal_ellipse

# get semantic version from setuptools-scm
__version__ = pyTMD.version.version
