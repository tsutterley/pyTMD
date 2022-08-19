#!/usr/bin/env python
u"""
nearest_extrap.py
Written by Tyler Sutterley (04/2022)
Uses kd-trees for nearest-neighbor extrapolation of tide model data

CALLING SEQUENCE:
    data = nearest_extrap(ilon,ilat,idata,lon,lat)

INPUTS:
    x: x-coordinates of tidal model
    y: y-coordinates of tidal model
    data: tide model data
    XI: output x-coordinates
    YI: output y-coordinates

OPTIONS:
    fill_value: invalid value
    dtype: output data type
    cutoff: return only neighbors within distance [km]
    EPSG: projection of tide model data

OUTPUT:
    DATA: extrapolated data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

PROGRAM DEPENDENCIES:
    spatial.py: utilities for reading and writing spatial data

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
        use valid tide model points when creating kd-trees
    Updated 02/2022: fix equirectangular case for cutoffs near poles
    Updated 05/2021: set ellipsoidal major axis to WGS84 in kilometers
    Updated 03/2021: add checks to prevent runtime exception
        where there are no valid points within the input bounds
        or no points to be extrapolated
        replaced numpy bool/int to prevent deprecation warnings
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Written 12/2020
"""
import numpy as np
import scipy.spatial
import pyTMD.spatial

# PURPOSE: Nearest-neighbor extrapolation of valid data to output data
def nearest_extrap(x, y, data, XI, YI, fill_value=np.nan,
    dtype=np.float64, cutoff=np.inf, EPSG='4326'):
    """
    Nearest-neighbor extrapolation of valid model data

    Parameters
    ----------
    x: float
        x-coordinates of tidal model
    y: float
        y-coordinates of tidal model
    data: float
        Tide model data
    XI: float
        Output x-coordinates
    YI: float
        Output y-coordinates
    fill_value: float, default np.nan
        Invalid value
    dtype: obj, default np.float64
        Output data type
    cutoff: float, default np.inf
        return only neighbors within distance [km]
    EPSG: str, default '4326'
        projection of tide model data

    Returns
    -------
    DATA: float
        interpolated data
    """
    # verify output dimensions
    XI = np.atleast_1d(XI)
    YI = np.atleast_1d(YI)
    # extrapolate valid data values to data
    npts = len(XI)
    # return none if no invalid points
    if (npts == 0):
        return

    # allocate to output extrapolate data array
    DATA = np.ma.zeros((npts), dtype=dtype, fill_value=fill_value)
    DATA.mask = np.ones((npts), dtype=bool)
    # initially set all data to fill value
    DATA.data[:] = data.fill_value

    # create combined valid mask
    valid_mask = (~data.mask) & np.isfinite(data.data)
    # reduce to model points within bounds of input points
    valid_bounds = np.ones_like(data.mask, dtype=bool)

    # calculate coordinates for nearest-neighbors
    if (EPSG == '4326'):
        # global or regional equirectangular model
        # calculate meshgrid of model coordinates
        gridlon,gridlat = np.meshgrid(x,y)
        # ellipsoidal major axis in kilometers
        a_axis = 6378.137
        # calculate Cartesian coordinates of input grid
        gridx,gridy,gridz = pyTMD.spatial.to_cartesian(gridlon,
            gridlat, a_axis=a_axis)
        # calculate Cartesian coordinates of output coordinates
        xs,ys,zs = pyTMD.spatial.to_cartesian(XI, YI, a_axis=a_axis)
        # range of output points in cartesian coordinates
        xmin,xmax = (np.min(xs), np.max(xs))
        ymin,ymax = (np.min(ys), np.max(ys))
        zmin,zmax = (np.min(zs), np.max(zs))
        # reduce to model points within bounds of input points
        valid_bounds = np.ones_like(data.mask, dtype=bool)
        valid_bounds &= (gridx >= (xmin-2.0*cutoff))
        valid_bounds &= (gridx <= (xmax+2.0*cutoff))
        valid_bounds &= (gridy >= (ymin-2.0*cutoff))
        valid_bounds &= (gridy <= (ymax+2.0*cutoff))
        valid_bounds &= (gridz >= (zmin-2.0*cutoff))
        valid_bounds &= (gridz <= (zmax+2.0*cutoff))
        # check if there are any valid points within the input bounds
        if not np.any(valid_mask & valid_bounds):
            # return filled masked array
            return DATA
        # find where input grid is valid and close to output points
        indy,indx = np.nonzero(valid_mask & valid_bounds)
        # create KD-tree of valid points
        tree = scipy.spatial.cKDTree(np.c_[gridx[indy,indx],
            gridy[indy,indx], gridz[indy,indx]])
        # flattened valid data array
        flattened = data.data[indy,indx]
        # output coordinates
        points = np.c_[xs,ys,zs]
    else:
        # projected model
        # calculate meshgrid of model coordinates
        gridx,gridy = np.meshgrid(x,y)
        # range of output points
        xmin,xmax = (np.min(XI),np.max(XI))
        ymin,ymax = (np.min(YI),np.max(YI))
        # reduce to model points within bounds of input points
        valid_bounds = np.ones_like(data.mask, dtype=bool)
        valid_bounds &= (gridx >= (xmin-2.0*cutoff))
        valid_bounds &= (gridx <= (xmax+2.0*cutoff))
        valid_bounds &= (gridy >= (ymin-2.0*cutoff))
        valid_bounds &= (gridy <= (ymax+2.0*cutoff))
        # check if there are any valid points within the input bounds
        if not np.any(valid_mask & valid_bounds):
            # return filled masked array
            return data
        # find where input grid is valid and close to output points
        indy,indx = np.nonzero(valid_mask & valid_bounds)
        # flattened model coordinates
        tree = scipy.spatial.cKDTree(np.c_[gridx[indy,indx],
            gridy[indy,indx]])
        # flattened valid data array
        flattened = data.data[indy,indx]
        # output coordinates
        points = np.c_[XI,YI]

    # query output data points and find nearest neighbor within cutoff
    dd,ii = tree.query(points, k=1, distance_upper_bound=cutoff)
    # spatially extrapolate using nearest neighbors
    if np.any(np.isfinite(dd)):
        ind, = np.nonzero(np.isfinite(dd))
        DATA.data[ind] = flattened[ii[ind]]
        DATA.mask[ind] = False
    # return extrapolated values
    return DATA

# PURPOSE: calculate Euclidean distances between points
def distance_matrix(c1, c2):
    """
    Calculate Euclidean distances between points

    Parameters
    ----------
    c1: float
        first set of coordinates
    c2: float
        second set of coordinates

    Returns
    -------
    c: float
        Euclidean distance
    """
    # decompose Euclidean distance: (x-y)^2 = x^2 - 2xy + y^2
    dx2 = np.sum(c1**2)
    dxy = np.dot(c1[np.newaxis,:], c2.T)
    dy2 = np.sum(c2**2, axis=1)
    # calculate Euclidean distance
    D, = np.sqrt(dx2 - 2.0*dxy + dy2)
    return D
