#!/usr/bin/env python
u"""
convert_crs.py
Written by Tyler Sutterley (03/2023)
Converts points to and from Coordinates Reference Systems (CRS)

CALLING SEQUENCE:
    x, y = convert_crs(lon, lat, PROJ, 'F')
    lon, lat = convert_crs(x, y, PROJ, 'B')

INPUTS:
    i1: longitude ('F') or projection easting x ('B')
    i2: latitude ('F') or projection northing y ('B')
    PROJ: spatial reference system code for coordinate transformations
    BF: backwards ('B') or forward ('F') translations

OPTIONS:
    EPSG: spatial reference system code for input (F) and output (B) coordinates

OUTPUTS:
    o1: projection easting x ('F') or longitude ('B')
    o2: projection northing y ('F') or latitude ('B')

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
        https://pyproj4.github.io/pyproj/

UPDATE HISTORY:
    Updated 03/2023: add basic variable typing to function inputs
        renamed coordinate reference system conversion functions
    Updated 02/2023: use named exception before passing to custom
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: added function for using custom projections
    Updated 06/2021: added 3413 for new 1km Greenland model from ESR
    Updated 08/2020: using conversion protocols following pyproj-2 updates
        https://pyproj4.github.io/pyproj/stable/gotchas.html
    Updated 07/2020: added function docstrings. changed function name
    Updated 03/2020: remove commented coordinate conversion functions
    Updated 11/2019: using pyproj for coordinate conversions
    Written 09/2017
"""
from __future__ import annotations

import logging
import numpy as np

# attempt imports
try:
    import pyproj
except (ImportError, ModuleNotFoundError) as exc:
    logging.critical("pyproj not available")

def convert_crs(
        i1: np.ndarray,
        i2: np.ndarray,
        PROJ: str,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts points to and from Coordinates Reference Systems (CRS)

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    PROJ: str
        Spatial reference system code for coordinate transformations
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    # python dictionary with named conversion functions
    conversion_functions = {}
    conversion_functions['3031'] = _EPSG3031
    conversion_functions['3413'] = _EPSG3413
    conversion_functions['CATS2008'] = _CATS2008
    conversion_functions['3976'] = _EPSG3976
    conversion_functions['PSNorth'] = _PSNorth
    conversion_functions['4326'] = _EPSG4326
    # check that PROJ for conversion was entered correctly
    # run named conversion program and return values
    try:
        o1, o2 = conversion_functions[PROJ](i1, i2, BF, EPSG=EPSG)
    except KeyError as exc:
        pass
    else:
        return (o1, o2)
    # try changing the projection using a custom projection
    # run custom conversion program and return values
    try:
        o1, o2 = _custom(i1, i2, PROJ, BF, EPSG=EPSG)
    except Exception as exc:
        pass
    else:
        return (o1, o2)
    # projection not found or available
    raise Exception(f'PROJ: {PROJ} conversion function not found')

# PURPOSE: try to get the projection information
def crs_from_input(PROJECTION: int | str):
    """
    Attempt to get the Coordinate Reference System for an input code

    Parameters
    ----------
    PROJECTION: int or str
        Coordinate Reference System code
    """
    # EPSG projection code
    try:
        crs = pyproj.CRS.from_epsg(int(PROJECTION))
    except (ValueError, pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # coordinate reference system string
    try:
        crs = pyproj.CRS.from_string(PROJECTION)
    except (ValueError, pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # no projection can be made
    raise pyproj.exceptions.CRSError

# wrapper function for models in EPSG 3031 (Antarctic Polar Stereographic)
def _EPSG3031(
        i1: np.ndarray,
        i2: np.ndarray,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts models in EPSG 3031 (Antarctic Polar Stereographic)

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    # projections for converting from input EPSG (default latitude/longitude)
    crs1 = crs_from_input(EPSG)
    crs2 = pyproj.CRS.from_user_input({'proj':'stere','lat_0':-90,'lat_ts':-71,
        'lon_0':0,'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84',
        'units':'km'})
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    # convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    # return the output variables
    return transformer.transform(i1, i2, direction=direction)

# wrapper function for models in EPSG 3413 (Sea Ice Polar Stereographic North)
def _EPSG3413(
        i1: np.ndarray,
        i2: np.ndarray,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts models in EPSG 3413 (Sea Ice Polar Stereographic North)

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    # projections for converting from input EPSG (default latitude/longitude)
    crs1 = crs_from_input(EPSG)
    crs2 = pyproj.CRS.from_user_input({'proj':'stere','lat_0':90,'lat_ts':70,
        'lon_0':-45,'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84',
        'units':'km'})
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    # convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    # return the output variables
    return transformer.transform(i1, i2, direction=direction)

# wrapper function for CATS2008 tide models
def _CATS2008(
        i1: np.ndarray,
        i2: np.ndarray,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts Circum-Antarctic Tidal Simulation models

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    # projections for converting from input EPSG (default latitude/longitude)
    crs1 = crs_from_input(EPSG)
    crs2 = pyproj.CRS.from_user_input({'proj':'stere','lat_0':-90,'lat_ts':-71,
        'lon_0':-70,'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84',
        'units':'km'})
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    # convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    # return the output variables
    return transformer.transform(i1, i2, direction=direction)

# wrapper function for models in EPSG 3976 (NSIDC Sea Ice Stereographic South)
def _EPSG3976(
        i1: np.ndarray,
        i2: np.ndarray,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts models in EPSG 3976 (Sea Ice Polar Stereographic South)

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    # projections for converting from input EPSG (default latitude/longitude)
    crs1 = crs_from_input(EPSG)
    crs2 = pyproj.CRS.from_user_input({'proj':'stere','lat_0':-90,'lat_ts':-70,
        'lon_0':0,'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84',
        'units':'km'})
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    # convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    # return the output variables
    return transformer.transform(i1, i2, direction=direction)

# wrapper function for models in (idealized) PSNorth projection
def _PSNorth(
        i1: np.ndarray,
        i2: np.ndarray,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts idealized Arctic Polar Stereographic models

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    # projections for converting to and from input EPSG
    crs1 = crs_from_input(EPSG)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # convert lat/lon to (idealized) Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
        lon, lat = transformer.transform(i1, i2, direction=direction)
        o1 = (90.0-lat)*111.7*np.cos(lon/180.0*np.pi)
        o2 = (90.0-lat)*111.7*np.sin(lon/180.0*np.pi)
    # convert (idealized) Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
        lon = np.arctan2(i2, i1)*180.0/np.pi
        lat = 90.0 - np.sqrt(i1**2+i2**2)/111.7
        ii, = np.nonzero(lon < 0)
        lon[ii] += 360.0
        o1, o2 = transformer.transform(lon, lat, direction=direction)
    # return the output variables
    return (o1, o2)

# wrapper function to pass lat/lon values or convert if EPSG
def _EPSG4326(
        i1: np.ndarray,
        i2: np.ndarray,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts models in EPSG 4326 (WGS84 Latitude/Longitude)

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    crs1 = crs_from_input(EPSG)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    # return the output variables
    return transformer.transform(i1, i2, direction=direction)

# wrapper function for using custom projections
def _custom(
        i1: np.ndarray,
        i2: np.ndarray,
        PROJ: int | str,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts models in a custom projection

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    PROJ: int or str
        Spatial reference system code for coordinate transformations
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    # projections for converting from input EPSG (default latitude/longitude)
    crs1 = crs_from_input(EPSG)
    crs2 = crs_from_input(PROJ)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # convert lat/lon to custom projection
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    # convert custom projection to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    # return the output variables
    return transformer.transform(i1, i2, direction=direction)
