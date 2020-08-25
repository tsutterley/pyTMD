#!/usr/bin/env python
u"""
convert_ll_xy.py (08/2020)
Wrapper function to convert lat/lon points to and from projected coordinates

CALLING SEQUENCE:
    x,y = convert_ll_xy(lon,lat,PROJ,'F')
    lon,lat = convert_ll_xy(x,y,PROJ,'B')

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
    Updated 08/2020: using conversion protocols following pyproj-2 updates
        https://pyproj4.github.io/pyproj/stable/gotchas.html
    Updated 07/2020: added function docstrings. changed function name
    Updated 03/2020: remove commented coordinate conversion functions
    Updated 11/2019: using pyproj for coordinate conversions
    Written 09/2017
"""
import numpy as np
import pyproj

def convert_ll_xy(i1,i2,PROJ,BF,EPSG=4326):
    """
    Converts lat/lon points to and from projected coordinates

    Arguments
    ---------
    i1: longitude ('F') or projection easting x ('B')
    i2: latitude ('F') or projection northing y ('B')
    PROJ: spatial reference system code for coordinate transformations
    BF: backwards ('B') or forward ('F') translations

    Keyword Arguments
    -----------------
    EPSG: spatial reference system code for input (F) and output (B) coordinates

    Returns
    -------
    o1: projection easting x ('F') or longitude ('B')
    o2: projection northing y ('F') or latitude ('B')
    """
    #-- python dictionary with conversion functions
    conversion_functions = {}
    conversion_functions['3031'] = xy_ll_EPSG3031
    conversion_functions['CATS2008'] = xy_ll_CATS2008
    conversion_functions['3976'] = xy_ll_EPSG3976
    conversion_functions['PSNorth'] = xy_ll_PSNorth
    conversion_functions['4326'] = xy_ll_EPSG4326
    #-- check that PROJ for conversion was entered correctly
    if PROJ not in conversion_functions.keys():
        raise Exception('PROJ:{0} conversion function not found'.format(PROJ))
    #-- run conversion program and return values
    o1,o2 = conversion_functions[PROJ](i1,i2,BF,EPSG=EPSG)
    return (o1,o2)

#-- wrapper function for models in EPSG 3031 (Antarctic Polar Stereographic)
def xy_ll_EPSG3031(i1,i2,BF,EPSG=4326):
    #-- projections for converting from input EPSG (default latitude/longitude)
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG))
    crs2 = pyproj.CRS.from_user_input({'proj':'stere','lat_0':-90,'lat_ts':-71,
        'lon_0':0,'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84',
        'units':'km'})
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    #-- convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    #-- return the output variables
    return transformer.transform(i1, i2, direction=direction)

#-- wrapper function for CATS2008 tide models
def xy_ll_CATS2008(i1,i2,BF,EPSG=4326):
    #-- projections for converting from input EPSG (default latitude/longitude)
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG))
    crs2 = pyproj.CRS.from_user_input({'proj':'stere','lat_0':-90,'lat_ts':-71,
        'lon_0':-70,'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84',
        'units':'km'})
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    #-- convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    #-- return the output variables
    return transformer.transform(i1, i2, direction=direction)

#-- wrapper function for models in EPSG 3976 (NSIDC Sea Ice Stereographic South)
def xy_ll_EPSG3976(i1,i2,BF,EPSG=4326):
    #-- projections for converting from input EPSG (default latitude/longitude)
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG))
    crs2 = pyproj.CRS.from_user_input({'proj':'stere','lat_0':-90,'lat_ts':-70,
        'lon_0':0,'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84',
        'units':'km'})
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    #-- convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    #-- return the output variables
    return transformer.transform(i1, i2, direction=direction)

#-- wrapper function for models in (idealized) PSNorth projection
def xy_ll_PSNorth(i1,i2,BF,EPSG=4326):
    #-- projections for converting to and from input EPSG
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG))
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- convert lat/lon to (idealized) Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
        lon,lat = transformer.transform(i1, i2, direction=direction)
        o1 = (90.0-lat)*111.7*np.cos(lon/180.0*np.pi)
        o2 = (90.0-lat)*111.7*np.sin(lon/180.0*np.pi)
    #-- convert (idealized) Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
        lon = np.arctan2(i2,i1)*180.0/np.pi
        lat = 90.0 - np.sqrt(i1**2+i2**2)/111.7
        ii, = np.nonzero(lon < 0)
        lon[ii] += 360.0
        o1,o2 = transformer.transform(lon, lat, direction=direction)
    #-- return the output variables
    return (o1,o2)

#-- wrapper function to pass lat/lon values or convert if EPSG
def xy_ll_EPSG4326(i1,i2,BF,EPSG=4326):
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG))
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    if (BF.upper() == 'F'):
        direction = pyproj.enums.TransformDirection.FORWARD
    elif (BF.upper() == 'B'):
        direction = pyproj.enums.TransformDirection.INVERSE
    #-- return the output variables
    return transformer.transform(i1, i2, direction=direction)
