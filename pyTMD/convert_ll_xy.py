#!/usr/bin/env python
u"""
convert_ll_xy.py (07/2020)
Wrapper function to convert lat/lon points to and from projected coordinates

CALLING SEQUENCE:
    x,y = convert_ll_xy(lon,lat,PROJ,'F')
    lon,lat = convert_ll_xy(x,y,PROJ,'B')

INPUTS:
    i1: longitude ('F') or projection easting x ('B')
    i2: latitude ('F') or projection northing y ('B')
    PROJ: spatial reference system code for coordinate transformations
    BF: backwards ('B') or forward ('F') translations

OUTPUTS:
    o1: projection easting x ('F') or longitude ('B')
    o2: projection northing y ('F') or latitude ('B')

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

UPDATE HISTORY:
    Updated 07/2020: added function docstrings. changed function name
    Updated 03/2020: remove commented coordinate conversion functions
    Updated 11/2019: using pyproj for coordinate conversions
    Written 09/2017
"""
import numpy as np
import pyproj

def convert_ll_xy(i1,i2,PROJ,BF):
    """
    Converts lat/lon points to and from projected coordinates

    Arguments
    ---------
    i1: longitude ('F') or projection easting x ('B')
    i2: latitude ('F') or projection northing y ('B')
    PROJ: spatial reference system code for coordinate transformations
    BF: backwards ('B') or forward ('F') translations

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
    conversion_functions['4326'] = pass_values
    #-- check that PROJ for conversion was entered correctly
    if PROJ not in conversion_functions.keys():
        raise Exception('PROJ:{0} conversion function not found'.format(PROJ))
    #-- run conversion program and return values
    o1,o2 = conversion_functions[PROJ](i1,i2,BF)
    return (o1,o2)

#-- wrapper function for models in EPSG 3031 (Antarctic Polar Stereographic)
def xy_ll_EPSG3031(i1,i2,BF):
    #-- projections for converting from latitude/longitude
    proj1 = pyproj.Proj("+init=EPSG:{0:d}".format(4326))
    proj2 = pyproj.Proj({'proj':'stere','lat_0':-90,'lat_ts':-71,'lon_0':0,
        'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84','units':'km'})
    #-- convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        o1,o2 = pyproj.transform(proj1, proj2, i1, i2)
    #-- convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        o1,o2 = pyproj.transform(proj2, proj1, i1, i2)
    #-- return the output variables
    return (o1,o2)

#-- wrapper function for CATS2008 tide models
def xy_ll_CATS2008(i1,i2,BF):
    #-- projections for converting from latitude/longitude
    proj1 = pyproj.Proj("+init=EPSG:{0:d}".format(4326))
    proj2 = pyproj.Proj({'proj':'stere','lat_0':-90,'lat_ts':-71,'lon_0':-70,
        'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84','units':'km'})
    #-- convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        o1,o2 = pyproj.transform(proj1, proj2, i1, i2)
    #-- convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        o1,o2 = pyproj.transform(proj2, proj1, i1, i2)
    #-- return the output variables
    return (o1,o2)

#-- wrapper function for models in EPSG 3976 (NSIDC Sea Ice Stereographic South)
def xy_ll_EPSG3976(i1,i2,BF):
    #-- projections for converting from latitude/longitude
    proj1 = pyproj.Proj("+init=EPSG:{0:d}".format(4326))
    proj2 = pyproj.Proj({'proj':'stere','lat_0':-90,'lat_ts':-70,'lon_0':0,
        'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84','units':'km'})
    #-- convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        o1,o2 = pyproj.transform(proj1, proj2, i1, i2)
    #-- convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        o1,o2 = pyproj.transform(proj2, proj1, i1, i2)
    #-- return the output variables
    return (o1,o2)

#-- wrapper function for models in PSNorth projection
def xy_ll_PSNorth(i1,i2,BF):
    # #-- projections for converting from latitude/longitude
    # proj1 = pyproj.Proj("+init=EPSG:{0:d}".format(4326))
    # proj2 = pyproj.Proj({'proj':'stere','lat_0':90,'lat_ts':90,'lon_0':270,
    #     'x_0':0.,'y_0':0.,'ellps': 'WGS84','datum': 'WGS84','units':'km'})
    #-- convert lat/lon to Polar-Stereographic x/y
    if (BF.upper() == 'F'):
        # o1,o2 = pyproj.transform(proj1, proj2, i1, i2)
        o1 = (90.0-i2)*111.7*np.cos(i1/180.0*np.pi)
        o2 = (90.0-i2)*111.7*np.sin(i1/180.0*np.pi)
    #-- convert Polar-Stereographic x/y to lat/lon
    elif (BF.upper() == 'B'):
        # o1,o2 = pyproj.transform(proj2, proj1, i1, i2)
        o1 = 90.0 - np.sqrt(i1**2+i2**2)/111.7
        o2 = np.arctan2(i2,i1)*180.0/np.pi
        ii, = np.nonzero(o1 < 0)
        o1[ii] += 360.0
    #-- return the output variables
    return (o1,o2)

#-- wrapper function to pass lat/lon values
def pass_values(i1,i2,BF):
    return (i1,i2)
