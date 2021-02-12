#!/usr/bin/env python
u"""
nearest_extrap.py (12/2020)
Uses kd-trees for nearest-neighbor extrapolation of valid model data

CALLING SEQUENCE:
    data = nearest_extrap(ilon,ilat,idata,lon,lat)

INPUTS:
    ilon: longitude of tidal model
    ilat: latitude of tidal model
    idata: tide model data
    lat: output latitude
    lon: output longitude

OPTIONS:
    fill_value: invalid value
    dtype: output data type
    cutoff: return only neighbors within distance [km]
    EPSG: projection of tide model data

OUTPUT:
    data: extrapolated data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

UPDATE HISTORY:
    Written 12/2020
"""
import numpy as np
import scipy.spatial

#-- PURPOSE: Nearest-neighbor extrapolation of valid data to output data
def nearest_extrap(ilon,ilat,idata,lon,lat,fill_value=np.nan,
    dtype=np.float,cutoff=np.inf,EPSG='4326'):
    """
    Nearest-neighbor extrapolation of valid model data

    Arguments
    ---------
    ilon: longitude of tidal model
    ilat: latitude of tidal model
    idata: tide model data
    lat: output latitude
    lon: output longitude

    Keyword arguments
    -----------------
    fill_value: invalid value
    dtype: output data type
    cutoff: return only neighbors within distance [km]
    EPSG: projection of tide model data

    Returns
    -------
    data: interpolated data
    """
    #-- grid step size of tide model
    dlon = np.abs(ilon[1] - ilon[0])
    dlat = np.abs(ilat[1] - ilat[0])
    #-- extrapolate valid data values to data
    npts = len(lon)
    #-- allocate to output extrapolate data array
    data = np.ma.zeros((npts),dtype=dtype,fill_value=fill_value)
    data.mask = np.ones((npts),dtype=np.bool)
    #-- initially set all data to fill value
    data.data[:] = data.fill_value
    #-- range of output points
    xmin,xmax = (np.min(lon),np.max(lon))
    ymin,ymax = (np.min(lat),np.max(lat))

    #-- calculate meshgrid of model coordinates
    gridlon,gridlat = np.meshgrid(ilon,ilat)
    #-- find where input grid is valid and close to output points
    indy,indx = np.nonzero((~idata.mask) & np.isfinite(idata.data) &
        (gridlon >= (xmin-2.0*dlon)) & (gridlon <= (xmax+2.0*dlon)) &
        (gridlat >= (ymin-2.0*dlat)) & (gridlat <= (ymax+2.0*dlat)))
    #-- flattened valid data array
    iflat = idata.data[indy,indx]

    #-- calculate coordinates for nearest-neighbors
    if (EPSG == '4326'):
        #-- valid grid latitude and longitude in radians
        iphi = np.pi*gridlon[indy,indx]/180.0
        itheta = np.pi*gridlat[indy,indx]/180.0
        #-- WGS84 ellipsoid parameters
        a_axis = 6378.1366
        f = 1.0/298.257223563
        ecc1 = np.sqrt((2.0*f - f**2)*a_axis**2)/a_axis
        #-- calculate Cartesian coordinates of input grid
        N = a_axis/np.sqrt(1.0-ecc1**2.0*np.sin(itheta)**2.0)
        xflat = N*np.cos(itheta)*np.cos(iphi)
        yflat = N*np.cos(itheta)*np.sin(iphi)
        zflat = (N*(1.0-ecc1**2.0))*np.sin(itheta)
        tree = scipy.spatial.cKDTree(np.c_[xflat,yflat,zflat])
        #-- calculate Cartesian coordinates of output coordinates
        ns = a_axis/np.sqrt(1.0-ecc1**2.0*np.sin(np.pi*lat/180.0)**2.0)
        xs = ns*np.cos(np.pi*lat/180.0)*np.cos(np.pi*lon/180.0)
        ys = ns*np.cos(np.pi*lat/180.0)*np.sin(np.pi*lon/180.0)
        zs = (ns*(1.0-ecc1**2.0))*np.sin(np.pi*lat/180.0)
        points = np.c_[xs,ys,zs]
    else:
        #-- flattened model coordinates
        tree = scipy.spatial.cKDTree(np.c_[gridlon[indy,indx],
            gridlat[indy,indx]])
        #-- output coordinates
        points = np.c_[lon,lat]

    #-- query output data points and find nearest neighbor within cutoff
    dd,ii = tree.query(points,k=1,distance_upper_bound=cutoff)
    ind, = np.nonzero(np.isfinite(dd))
    data.data[ind] = iflat[ii[ind]]
    data.mask[ind] = False
    #-- return extrapolated values
    return data

#-- PURPOSE: calculate Euclidean distances between points
def distance_matrix(c1,c2):
    """
    Calculate Euclidean distances between points

    Arguments
    ---------
    c1: first set of coordinates
    c2: second set of coordinates

    Returns
    -------
    c: Euclidean distance
    """
    #-- decompose Euclidean distance: (x-y)^2 = x^2 - 2xy + y^2
    dx2 = np.sum(c1**2)
    dxy = np.dot(c1[np.newaxis,:], c2.T)
    dy2 = np.sum(c2**2, axis=1)
    #-- calculate Euclidean distance
    D, = np.sqrt(dx2 - 2.0*dxy + dy2)
    return D
