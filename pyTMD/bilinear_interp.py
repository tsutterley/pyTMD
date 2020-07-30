#!/usr/bin/env python
u"""
bilinear_interp.py (07/2020)
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
    dtype: output data type

OUTPUT:
    data: interpolated data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 07/2020: split into separate function
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Updated 09/2017: Rewritten in Python
"""
import numpy as np

#-- PURPOSE: bilinear interpolation of input data to output data
def bilinear_interp(ilon,ilat,idata,lon,lat,dtype=np.float):
    """
    Bilinear interpolation of input data to output coordinates

    Arguments
    ---------
    ilon: longitude of tidal model
    ilat: latitude of tidal model
    idata: tide model data
    lat: output latitude
    lon: output longitude

    Keyword arguments
    -----------------
    dtype: output data type

    Returns
    -------
    data: interpolated data
    """
    #-- degrees to radians
    dtr = np.pi/180.0
    #-- grid step size of tide model
    dlon = np.abs(ilon[1] - ilon[0])
    dlat = np.abs(ilat[1] - ilat[0])
    #-- Convert input coordinates to radians
    phi = ilon*dtr
    th = (90.0 - ilat)*dtr
    #-- Convert output data coordinates to radians
    xphi = lon*dtr
    xth = (90.0 - lat)*dtr
    #-- interpolate gridded data values to data
    data = np.zeros_like(lon,dtype=dtype)
    for i,l in enumerate(lon):
        #-- calculating the indices for the original grid
        dx = (ilon - np.floor(lon[i]/dlon)*dlon)**2
        dy = (ilat - np.floor(lat[i]/dlat)*dlat)**2
        iph = np.argmin(dx)
        ith = np.argmin(dy)
        #-- if on corner value: use exact
        if ((lat[i] == ilat[ith]) & (lon[i] == ilon[iph])):
            data[i] = idata[ith,iph]
        elif ((lat[i] == ilat[ith+1]) & (lon[i] == ilon[iph])):
            data[i] = idata[ith+1,iph]
        elif ((lat[i] == ilat[ith]) & (lon[i] == ilon[iph+1])):
            data[i] = idata[ith,iph+1]
        elif ((lat[i] == ilat[ith+1]) & (lon[i] == ilon[iph+1])):
            data[i] = idata[ith+1,iph+1]
        else:
            #-- corner weight values for i,j
            Wa = (xphi[i]-phi[iph])*(xth[i]-th[ith])
            Wb = (phi[iph+1]-xphi[i])*(xth[i]-th[ith])
            Wc = (xphi[i]-phi[iph])*(th[ith+1]-xth[i])
            Wd = (phi[iph+1]-xphi[i])*(th[ith+1]-xth[i])
            #-- divisor weight value
            W = (phi[iph+1]-phi[iph])*(th[ith+1]-th[ith])
            #-- corner data values for i,j
            Ia = idata[ith,iph]#-- (0,0)
            Ib = idata[ith,iph+1]#-- (1,0)
            Ic = idata[ith+1,iph]#-- (0,1)
            Id = idata[ith+1,iph+1]#-- (1,1)
            #-- calculate interpolated value for i
            data[i] = (Ia*Wa + Ib*Wb + Ic*Wc + Id*Wd)/W
    #-- return interpolated values
    return data
