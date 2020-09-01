#!/usr/bin/env python
u"""
read_GOT_model.py (08/2020)
Reads files for Richard Ray's Global Ocean Tide (GOT) models and makes initial
    calculations to run the tide program
Includes functions to extract tidal harmonic constants out of a tidal model for
    given locations

INPUTS:
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    directory: data directory for tide data files
    model_files: list of model files for each constituent

OPTIONS:
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations
    GZIP: input files are compressed
    SCALE: scaling factor for converting to output units

OUTPUTS:
    amplitude: amplitudes of tidal constituents
    phase: phases of tidal constituents

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

PROGRAM DEPENDENCIES:
    bilinear_interp.py: bilinear interpolation of data to specified coordinates

UPDATE HISTORY:
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Updated 07/2020: added function docstrings. separate bilinear interpolation
        update griddata interpolation. add option GZIP for compression
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Updated 11/2019: find invalid mask points for each constituent
    Updated 09/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 07/2019: interpolate fill value mask with bivariate splines
    Updated 12/2018: python3 compatibility updates for division and zip
    Updated 10/2018: added SCALE as load tides are in mm and ocean are in cm
    Updated 08/2018: added multivariate spline interpolation option
    Written 07/2018
"""
from __future__ import division

import os
import gzip
import numpy as np
import scipy.interpolate
from pyTMD.bilinear_interp import bilinear_interp

#-- PURPOSE: extract tidal harmonic constants out of GOT model at coordinates
def extract_GOT_constants(ilon, ilat, directory, model_files,
    METHOD=None, GZIP=True, SCALE=1):
    """
    Reads files for Richard Ray's Global Ocean Tide (GOT) models
    Makes initial calculations to run the tide program
    Spatially interpolates tidal constituents to input coordinates

    Arguments
    ---------
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    directory: data directory for tide data files
    model_files: list of model files for each constituent

    Keyword arguments
    -----------------
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations
    GZIP: input files are compressed
    SCALE: scaling factor for converting to output units

    Returns
    -------
    amplitude: amplitudes of tidal constituents
    phase: phases of tidal constituents
    """
    #-- adjust longitudinal convention of input latitude and longitude
    #-- to fit tide model convention
    if (np.min(ilon) < 0.0):
        lt0, = np.nonzero(ilon < 0)
        ilon[lt0] += 360.0

    #-- number of points
    npts = len(ilon)
    #-- amplitude and phase
    nc = len(model_files)
    amplitude = np.ma.zeros((npts,nc))
    amplitude.mask = np.zeros((npts,nc),dtype=np.bool)
    phase = np.ma.zeros((npts,nc))
    phase.mask = np.zeros((npts,nc),dtype=np.bool)
    #-- read and interpolate each constituent
    for i,model_file in enumerate(model_files):
        #-- read constituent from elevation file
        hc,lon,lat = read_GOT_grid(os.path.join(directory,model_file),
            GZIP=GZIP)
        #-- grid step size of tide model
        dlon = np.abs(lon[1] - lon[0])
        dlat = np.abs(lat[1] - lat[0])
        #-- replace original values with extend matrices
        lon = extend_array(lon,dlon)
        hc = extend_matrix(hc)
        #-- interpolated complex form of constituent oscillation
        hci = np.ma.zeros((npts),dtype=hc.dtype,fill_value=hc.fill_value)
        hci.mask = np.zeros((npts),dtype=np.bool)
        #-- interpolate amplitude and phase of the constituent
        if (METHOD == 'bilinear'):
            #-- replace invalid values with nan
            hc[hc.mask] = np.nan
            #-- use quick bilinear to interpolate values
            hci.data[:] = bilinear_interp(lon,lat,hc,ilon,ilat,dtype=hc.dtype)
            #-- replace nan values with fill_value
            hci.mask[:] |= np.isnan(hci.data)
            hci.data[hci.mask] = hci.fill_value
        elif (METHOD == 'spline'):
            #-- interpolate complex form of the constituent with scipy
            f1=scipy.interpolate.RectBivariateSpline(lon,lat,
                hc.data.real.T,kx=1,ky=1)
            f2=scipy.interpolate.RectBivariateSpline(lon,lat,
                hc.data.imag.T,kx=1,ky=1)
            f3=scipy.interpolate.RectBivariateSpline(lon,lat,
                hc.mask.T,kx=1,ky=1)
            hci.data.real[:] = f1.ev(ilon,ilat)
            hci.data.imag[:] = f2.ev(ilon,ilat)
            hci.mask[:] = f3.ev(ilon,ilat).astype(np.bool)
        else:
            #-- use scipy regular grid to interpolate values for a given method
            r1 = scipy.interpolate.RegularGridInterpolator((lat,lon),
                hc.data, method=METHOD)
            r2 = scipy.interpolate.RegularGridInterpolator((lat,lon),
                hc.mask, method=METHOD)
            hci.data[:] = r1.__call__(np.c_[ilat,ilon])
            hci.mask[:] = np.ceil(r2.__call__(np.c_[ilat,ilon])).astype(np.bool)
        #-- convert amplitude from input units to meters
        amplitude.data[:,i] = np.abs(hci)*SCALE
        amplitude.mask[:,i] = np.copy(hci.mask)
        #-- convert phase to degrees
        phase.data[:,i] = np.arctan2(-np.imag(hci),np.real(hci))*180.0/np.pi
        phase.mask[:,i] = np.copy(hci.mask)
        phase.data[phase.data < 0] += 360.0

    #-- replace data for invalid mask values
    amplitude.data[amplitude.mask] = amplitude.fill_value
    phase.data[phase.mask] = phase.fill_value
    #-- return the interpolated values
    return (amplitude,phase)

#-- PURPOSE: wrapper function to extend an array
def extend_array(input_array,step_size):
    """
    Wrapper function to extend an array

    Arguments
    ---------
    input_array: array to extend
    step_size: step size between elements of array

    Returns
    -------
    temp: extended array
    """
    n = len(input_array)
    temp = np.zeros((n+3),dtype=input_array.dtype)
    #-- extended array [x-1,x0,...,xN,xN+1,xN+2]
    temp[0] = input_array[0] - step_size
    temp[1:-2] = input_array[:]
    temp[-2] = input_array[-1] + step_size
    temp[-1] = input_array[-1] + 2.0*step_size
    return temp

#-- PURPOSE: wrapper function to extend a matrix
def extend_matrix(input_matrix):
    """
    Wrapper function to extend a matrix

    Arguments
    ---------
    input_matrix: matrix to extend

    Returns
    -------
    temp: extended matrix
    """
    ny,nx = np.shape(input_matrix)
    temp = np.ma.zeros((ny,nx+3),dtype=input_matrix.dtype)
    temp[:,0] = input_matrix[:,-1]
    temp[:,1:-2] = input_matrix[:,:]
    temp[:,-2] = input_matrix[:,0]
    temp[:,-1] = input_matrix[:,1]
    return temp

#-- PURPOSE: read GOT model grid files
def read_GOT_grid(input_file, GZIP=False):
    """
    Read Richard Ray's Global Ocean Tide (GOT) model file

    Arguments
    ---------
    input_file: model file

    Returns
    -------
    hc: complex form of tidal constituent oscillation
    lon: longitude of tidal model
    lat: latitude of tidal model
    """
    #-- read input tide model file
    if GZIP:
        with gzip.open(os.path.expanduser(input_file),'rb') as f:
            file_contents = f.read().splitlines()
    else:
        with open(os.path.expanduser(input_file),'r') as f:
            file_contents = f.read().splitlines()
    #-- parse header text
    nlat,nlon = np.array(file_contents[2].split(), dtype=np.int)
    #-- longitude range
    ilat = np.array(file_contents[3].split(), dtype=np.float)
    #-- latitude range
    ilon = np.array(file_contents[4].split(), dtype=np.float)
    #-- mask fill value
    fill_value = np.array(file_contents[5].split(), dtype=np.float)
    #-- create output variables
    lat = np.linspace(ilat[0],ilat[1],nlat)
    lon = np.linspace(ilon[0],ilon[1],nlon)
    amp = np.ma.zeros((nlat,nlon),fill_value=fill_value[0],dtype=np.float32)
    ph = np.ma.zeros((nlat,nlon),fill_value=fill_value[0],dtype=np.float32)
    #-- create masks for output variables (0=valid)
    amp.mask = np.zeros((nlat,nlon),dtype=np.bool)
    ph.mask = np.zeros((nlat,nlon),dtype=np.bool)
    #-- starting lines to fill amplitude and phase variables
    l1 = 7
    l2 = 14 + np.int(nlon//11)*nlat + nlat
    #-- for each latitude
    for i in range(nlat):
        for j in range(nlon//11):
            j1 = j*11
            amp.data[i,j1:j1+11] = np.array(file_contents[l1].split(),dtype='f')
            ph.data[i,j1:j1+11] = np.array(file_contents[l2].split(),dtype='f')
            l1 += 1
            l2 += 1
        #-- add last tidal variables
        j1 = (j+1)*11; j2 = nlon % 11
        amp.data[i,j1:j1+j2] = np.array(file_contents[l1].split(),dtype='f')
        ph.data[i,j1:j1+j2] = np.array(file_contents[l2].split(),dtype='f')
        l1 += 1
        l2 += 1
    #-- calculate complex form of constituent oscillation
    hc = amp*np.exp(-1j*ph*np.pi/180.0)
    #-- set masks
    hc.mask = (amp.data == amp.fill_value) | (ph.data == ph.fill_value)
    #-- return output variables
    return (hc,lon,lat)
