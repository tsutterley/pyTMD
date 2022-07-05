#!/usr/bin/env python
u"""
bilinear_interp.py
Written by Tyler Sutterley (07/2022)
Bilinear interpolation of input data to output coordinates

CALLING SEQUENCE:
    data = bilinear_interp(ilon,ilat,idata,lon,lat)

INPUTS:
    ilon: longitude of tidal model
    ilat: latitude of tidal model
    idata: tide model data
    lat: output latitude
    lon: output longitude

OPTIONS:
    fill_value: invalid value
    dtype: output data type

OUTPUT:
    data: interpolated data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 07/2022: verify that input data is a masked array
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 12/2020: using numpy isclose to check corner points
    Updated 08/2020: check that output coordinates are within bounds
        allow small extrapolations if individual grid cells are invalid
    Updated 07/2020: split into separate function
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Updated 09/2017: Rewritten in Python
"""
import numpy as np

# PURPOSE: bilinear interpolation of input data to output data
def bilinear_interp(ilon, ilat, idata, lon, lat,
    fill_value=np.nan, dtype=np.float64):
    """
    Bilinear interpolation of input data to output coordinates

    Parameters
    ----------
    ilon: float
        longitude of tidal model
    ilat: float
        latitude of tidal model
    idata: float
        tide model data
    lat: float
        output latitude
    lon: float
        output longitude
    fill_value: float, default np.nan
        invalid value
    dtype: obj, default np.float64
        output data type

    Returns
    -------
    data: float
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
    data = np.ma.zeros((npts),dtype=dtype,fill_value=fill_value)
    data.mask = np.ones((npts),dtype=bool)
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
