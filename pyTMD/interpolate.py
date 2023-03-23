#!/usr/bin/env python
u"""
interpolate.py
Written by Tyler Sutterley (12/2022)
Interpolators for spatial data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

UPDATE HISTORY:
    Written 12/2022
"""
from __future__ import annotations

import numpy as np
import scipy.spatial
import scipy.interpolate
import pyTMD.spatial

# PURPOSE: bilinear interpolation of input data to output data
def bilinear(
        ilon: np.ndarray,
        ilat: np.ndarray,
        idata: np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        fill_value: float = np.nan,
        dtype: str | np.dtype = np.float64
    ):
    """
    Bilinear interpolation of input data to output coordinates

    Parameters
    ----------
    ilon: np.ndarray
        longitude of tidal model
    ilat: np.ndarray
        latitude of tidal model
    idata: np.ndarray
        tide model data
    lat: np.ndarray
        output latitude
    lon: np.ndarray
        output longitude
    fill_value: float, default np.nan
        invalid value
    dtype: np.dtype, default np.float64
        output data type

    Returns
    -------
    data: np.ndarray
        interpolated data
    """
    # verify that input data is masked array
    if not isinstance(idata, np.ma.MaskedArray):
        idata = np.ma.array(idata)
        idata.mask = np.zeros_like(idata, dtype=bool)
    # find valid points (within bounds)
    valid, = np.nonzero((lon >= ilon.min()) & (lon <= ilon.max()) &
        (lat > ilat.min()) & (lat < ilat.max()))
    # interpolate gridded data values to data
    npts = len(lon)
    # allocate to output interpolated data array
    data = np.ma.zeros((npts), dtype=dtype, fill_value=fill_value)
    data.mask = np.ones((npts), dtype=bool)
    # initially set all data to fill value
    data.data[:] = data.fill_value
    # for each valid point
    for i in valid:
        # calculating the indices for the original grid
        ix, = np.nonzero((ilon[0:-1] <= lon[i]) & (ilon[1:] > lon[i]))
        iy, = np.nonzero((ilat[0:-1] <= lat[i]) & (ilat[1:] > lat[i]))
        # corner data values for adjacent grid cells
        IM = np.ma.zeros((4), fill_value=fill_value, dtype=dtype)
        IM.mask = np.ones((4), dtype=bool)
        # corner weight values for adjacent grid cells
        WM = np.zeros((4))
        # build data and weight arrays
        for j,XI,YI in zip([0,1,2,3],[ix,ix+1,ix,ix+1],[iy,iy,iy+1,iy+1]):
            IM.data[j], = idata.data[YI,XI].astype(dtype)
            IM.mask[j], = idata.mask[YI,XI]
            WM[3-j], = np.abs(lon[i]-ilon[XI])*np.abs(lat[i]-ilat[YI])
        # if on corner value: use exact
        if (np.isclose(lat[i],ilat[iy]) & np.isclose(lon[i],ilon[ix])):
            data.data[i] = idata.data[iy,ix].astype(dtype)
            data.mask[i] = idata.mask[iy,ix]
        elif (np.isclose(lat[i],ilat[iy+1]) & np.isclose(lon[i],ilon[ix])):
            data.data[i] = idata.data[iy+1,ix].astype(dtype)
            data.mask[i] = idata.mask[iy+1,ix]
        elif (np.isclose(lat[i],ilat[iy]) & np.isclose(lon[i],ilon[ix+1])):
            data.data[i] = idata.data[iy,ix+1].astype(dtype)
            data.mask[i] = idata.mask[iy,ix+1]
        elif (np.isclose(lat[i],ilat[iy+1]) & np.isclose(lon[i],ilon[ix+1])):
            data.data[i] = idata.data[iy+1,ix+1].astype(dtype)
            data.mask[i] = idata.mask[iy+1,ix+1]
        elif np.any(np.isfinite(IM) & (~IM.mask)):
            # find valid indices for data summation and weight matrix
            ii, = np.nonzero(np.isfinite(IM) & (~IM.mask))
            # calculate interpolated value for i
            data.data[i] = np.sum(WM[ii]*IM[ii])/np.sum(WM[ii])
            data.mask[i] = np.all(IM.mask[ii])
    # return interpolated values
    return data


def spline(
        ilon: np.ndarray,
        ilat: np.ndarray,
        idata: np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        fill_value: float = None,
        dtype: str | np.dtype = np.float64,
        reducer=np.ceil,
        **kwargs
    ):
    """
    `Bivariate spline interpolation
    <https://docs.scipy.org/doc/scipy/reference/generated/
    scipy.interpolate.RectBivariateSpline.html>`_
    of input data to output coordinates

    Parameters
    ----------
    ilon: np.ndarray
        longitude of tidal model
    ilat: np.ndarray
        latitude of tidal model
    idata: np.ndarray
        tide model data
    lat: np.ndarray
        output latitude
    lon: np.ndarray
        output longitude
    fill_value: float or NoneType, default None
        invalid value
    dtype: np.dtype, default np.float64
        output data type
    reducer: obj, default np.ceil
        operation for converting mask to boolean
    kx: int, default 1
        degree of the bivariate spline in the x-dimension
    ky: int, default 1
        degree of the bivariate spline in the y-dimension
    kwargs: dict
        additional arguments for ``scipy.interpolate.RectBivariateSpline``

    Returns
    -------
    data: np.ndarray
        interpolated data
    """
    # set default keyword arguments
    kwargs.setdefault('kx', 1)
    kwargs.setdefault('ky', 1)
    # verify that input data is masked array
    if not isinstance(idata, np.ma.MaskedArray):
        idata = np.ma.array(idata)
        idata.mask = np.zeros_like(idata, dtype=bool)
    # interpolate gridded data values to data
    npts = len(lon)
    # allocate to output interpolated data array
    data = np.ma.zeros((npts), dtype=dtype, fill_value=fill_value)
    data.mask = np.ones((npts), dtype=bool)
    # construct splines for input data and mask
    if np.iscomplexobj(idata):
        s1 = scipy.interpolate.RectBivariateSpline(ilon, ilat,
            idata.data.real.T, **kwargs)
        s2 = scipy.interpolate.RectBivariateSpline(ilon, ilat,
            idata.data.imag.T, **kwargs)
        s3 = scipy.interpolate.RectBivariateSpline(ilon, ilat,
            idata.mask.T, **kwargs)
        # evaluate the spline at input coordinates
        data.data.real[:] = s1.ev(lon, lat)
        data.data.imag[:] = s2.ev(lon, lat)
        data.mask[:] = reducer(s3.ev(lon, lat)).astype(bool)
    else:
        s1 = scipy.interpolate.RectBivariateSpline(ilon, ilat,
            idata.data.T, **kwargs)
        s2 = scipy.interpolate.RectBivariateSpline(ilon, ilat,
            idata.mask.T, **kwargs)
        # evaluate the spline at input coordinates
        data.data[:] = s1.ev(lon, lat).astype(dtype)
        data.mask[:] = reducer(s2.ev(lon, lat)).astype(bool)
    # return interpolated values
    return data

def regulargrid(
        ilon: np.ndarray,
        ilat: np.ndarray,
        idata: np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        fill_value: float = None,
        dtype: str | np.dtype = np.float64,
        reducer=np.ceil,
        **kwargs
    ):
    """
    `Regular grid interpolation
    <https://docs.scipy.org/doc/scipy/reference/generated/
    scipy.interpolate.RegularGridInterpolator.html>`_
    of input data to output coordinates

    Parameters
    ----------
    ilon: np.ndarray
        longitude of tidal model
    ilat: np.ndarray
        latitude of tidal model
    idata: np.ndarray
        tide model data
    lat: np.ndarray
        output latitude
    lon: np.ndarray
        output longitude
    fill_value: float or NoneType, default None
        invalid value
    dtype: np.dtype, default np.float64
        output data type
    reducer: obj, default np.ceil
        operation for converting mask to boolean
    bounds_error: bool, default False
        raise Exception when values are requested outside domain
    method: str, default 'linear'
        Method of interpolation

            - ``'linear'``
            - ``'nearest'``
            - ``'slinear'``
            - ``'cubic'``
            - ``'quintic'``
    kwargs: dict
        additional arguments for ``scipy.interpolate.RegularGridInterpolator``

    Returns
    -------
    data: np.ndarray
        interpolated data
    """
    # set default keyword arguments
    kwargs.setdefault('bounds_error', False)
    kwargs.setdefault('method', 'linear')
    # verify that input data is masked array
    if not isinstance(idata, np.ma.MaskedArray):
        idata = np.ma.array(idata)
        idata.mask = np.zeros_like(idata, dtype=bool)
    # interpolate gridded data values to data
    npts = len(lon)
    # allocate to output interpolated data array
    data = np.ma.zeros((npts), dtype=dtype, fill_value=fill_value)
    data.mask = np.ones((npts), dtype=bool)
    # use scipy regular grid to interpolate values for a given method
    r1 = scipy.interpolate.RegularGridInterpolator((ilat, ilon),
        idata.data, fill_value=fill_value, **kwargs)
    r2 = scipy.interpolate.RegularGridInterpolator((ilat, ilon),
        idata.mask, fill_value=1, **kwargs)
    # evaluate the interpolator at input coordinates
    data.data[:] = r1.__call__(np.c_[lat, lon])
    data.mask[:] = reducer(r2.__call__(np.c_[lat, lon])).astype(bool)
    # return interpolated values
    return data

# PURPOSE: Nearest-neighbor extrapolation of valid data to output data
def extrapolate(
        ilon: np.ndarray,
        ilat: np.ndarray,
        idata: np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        fill_value: float = None,
        dtype: str | np.dtype = np.float64,
        cutoff: int | float = np.inf,
        EPSG: str or int = '4326',
        **kwargs
    ):
    """
    Nearest-neighbor (`NN`) extrapolation of valid model data using `kd-trees
    <https://docs.scipy.org/doc/scipy/reference/generated/
    scipy.spatial.cKDTree.html>`_

    Parameters
    ----------
    x: np.ndarray
        x-coordinates of tidal model
    y: np.ndarray
        y-coordinates of tidal model
    data: np.ndarray
        Tide model data
    XI: np.ndarray
        Output x-coordinates
    YI: np.ndarray
        Output y-coordinates
    fill_value: float, default np.nan
        Invalid value
    dtype: np.dtype, default np.float64
        Output data type
    cutoff: float, default np.inf
        return only neighbors within distance [km]

        Set to ``np.inf`` to extrapolate for all points
    EPSG: str, default '4326'
        projection of tide model data

    Returns
    -------
    DATA: np.ndarray
        interpolated data
    """
    # verify output dimensions
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    # extrapolate valid data values to data
    npts = len(lon)
    # return none if no invalid points
    if (npts == 0):
        return

    # allocate to output extrapolate data array
    data = np.ma.zeros((npts), dtype=dtype, fill_value=fill_value)
    data.mask = np.ones((npts), dtype=bool)
    # initially set all data to fill value
    data.data[:] = idata.fill_value

    # create combined valid mask
    valid_mask = (~idata.mask) & np.isfinite(idata.data)
    # reduce to model points within bounds of input points
    valid_bounds = np.ones_like(idata.mask, dtype=bool)

    # calculate coordinates for nearest-neighbors
    if (EPSG == '4326'):
        # global or regional equirectangular model
        # calculate meshgrid of model coordinates
        gridlon, gridlat = np.meshgrid(ilon, ilat)
        # ellipsoidal major axis in kilometers
        a_axis = 6378.137
        # calculate Cartesian coordinates of input grid
        gridx, gridy, gridz = pyTMD.spatial.to_cartesian(
            gridlon, gridlat, a_axis=a_axis)
        # calculate Cartesian coordinates of output coordinates
        xs, ys, zs = pyTMD.spatial.to_cartesian(
            lon, lat, a_axis=a_axis)
        # range of output points in cartesian coordinates
        xmin, xmax = (np.min(xs), np.max(xs))
        ymin, ymax = (np.min(ys), np.max(ys))
        zmin, zmax = (np.min(zs), np.max(zs))
        # reduce to model points within bounds of input points
        valid_bounds = np.ones_like(idata.mask, dtype=bool)
        valid_bounds &= (gridx >= (xmin - 2.0*cutoff))
        valid_bounds &= (gridx <= (xmax + 2.0*cutoff))
        valid_bounds &= (gridy >= (ymin - 2.0*cutoff))
        valid_bounds &= (gridy <= (ymax + 2.0*cutoff))
        valid_bounds &= (gridz >= (zmin - 2.0*cutoff))
        valid_bounds &= (gridz <= (zmax + 2.0*cutoff))
        # check if there are any valid points within the input bounds
        if not np.any(valid_mask & valid_bounds):
            # return filled masked array
            return data
        # find where input grid is valid and close to output points
        indy, indx = np.nonzero(valid_mask & valid_bounds)
        # create KD-tree of valid points
        tree = scipy.spatial.cKDTree(np.c_[gridx[indy, indx],
            gridy[indy, indx], gridz[indy, indx]])
        # flattened valid data array
        flattened = idata.data[indy, indx]
        # output coordinates
        points = np.c_[xs, ys, zs]
    else:
        # projected model
        # calculate meshgrid of model coordinates
        gridx, gridy = np.meshgrid(ilon, ilat)
        # range of output points
        xmin, xmax = (np.min(lon), np.max(lon))
        ymin, ymax = (np.min(lat), np.max(lat))
        # reduce to model points within bounds of input points
        valid_bounds = np.ones_like(idata.mask, dtype=bool)
        valid_bounds &= (gridx >= (xmin - 2.0*cutoff))
        valid_bounds &= (gridx <= (xmax + 2.0*cutoff))
        valid_bounds &= (gridy >= (ymin - 2.0*cutoff))
        valid_bounds &= (gridy <= (ymax + 2.0*cutoff))
        # check if there are any valid points within the input bounds
        if not np.any(valid_mask & valid_bounds):
            # return filled masked array
            return data
        # find where input grid is valid and close to output points
        indy, indx = np.nonzero(valid_mask & valid_bounds)
        # flattened model coordinates
        tree = scipy.spatial.cKDTree(np.c_[gridx[indy, indx],
            gridy[indy, indx]])
        # flattened valid data array
        flattened = idata.data[indy, indx]
        # output coordinates
        points = np.c_[lon, lat]

    # query output data points and find nearest neighbor within cutoff
    dd, ii = tree.query(points, k=1, distance_upper_bound=cutoff)
    # spatially extrapolate using nearest neighbors
    if np.any(np.isfinite(dd)):
        ind, = np.nonzero(np.isfinite(dd))
        data.data[ind] = flattened[ii[ind]]
        data.mask[ind] = False
    # return extrapolated values
    return data

# PURPOSE: calculate Euclidean distances between points
def _distance(c1: np.ndarray, c2: np.ndarray):
    """
    Calculate Euclidean distances between points

    Parameters
    ----------
    c1: np.ndarray
        first set of coordinates
    c2: np.ndarray
        second set of coordinates

    Returns
    -------
    c: np.ndarray
        Euclidean distance
    """
    # decompose Euclidean distance: (x-y)^2 = x^2 - 2xy + y^2
    dx2 = np.sum(c1**2)
    dxy = np.dot(c1[np.newaxis,:], c2.T)
    dy2 = np.sum(c2**2, axis=1)
    # calculate Euclidean distance
    D, = np.sqrt(dx2 - 2.0*dxy + dy2)
    return D
