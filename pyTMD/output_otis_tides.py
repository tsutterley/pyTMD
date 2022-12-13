#!/usr/bin/env python
u"""
output_otis_tides.py
Written by Tyler Sutterley (04/2022)
Writes OTIS-format tide files for use in the tide program
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2020: python3 compatibility updates for struct and utf8 encoding
    Updated 07/2020: added function docstrings
    Written 08/2018
"""
import warnings
import pyTMD.io

# PURPOSE: output grid file in OTIS format
def output_otis_grid(*args):
    """
    Writes OTIS-format grid files

    Parameters
    ----------
    FILE: str
        output OTIS grid file name
    xlim: float
        x-coordinate grid-cell edges of output grid
    ylim: float
        y-coordinate grid-cell edges of output grid
    hz:float
        bathymetry
    mz: int
        land/water mask
    iob: int
        open boundary index
    dt: float
        time step
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    return pyTMD.io.OTIS.output_otis_grid(*args)

# PURPOSE: output elevation file in OTIS format
def output_otis_elevation(*args):
    """
    Writes OTIS-format elevation files

    Parameters
    ----------
    FILE: str
        output OTIS elevation file name
    h: complex
        Eulerian form of tidal height oscillation
    xlim: float
        x-coordinate grid-cell edges of output grid
    ylim: float
        y-coordinate grid-cell edges of output grid
    constituents: list
        tidal constituent IDs
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    return pyTMD.io.OTIS.output_otis_elevation(*args)

# PURPOSE: output transport file in OTIS format
def output_otis_transport(*args):
    """
    Writes OTIS-format transport files

    Parameters
    ----------
    FILE: str
        output OTIS transport file name
    u: complex
        Eulerian form of tidal zonal transport oscillation
    v: complex
        Eulerian form of tidal meridional transport oscillation
    xlim: float
        x-coordinate grid-cell edges of output grid
    ylim: float
        y-coordinate grid-cell edges of output grid
    constituents: list
        tidal constituent IDs
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    return pyTMD.io.OTIS.output_otis_transport(*args)
