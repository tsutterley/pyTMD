#!/usr/bin/env python
u"""
read_GOT_model.py (11/2022)
Reads files for Richard Ray's Global Ocean Tide (GOT) models and makes initial
    calculations to run the tide program
Includes functions to extract tidal harmonic constants out of a tidal model for
    given locations

INPUTS:
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    model_files: list of model files for each constituent

OPTIONS:
    method: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations
    extrapolate: extrapolate model using nearest-neighbors
    cutoff: extrapolation cutoff in kilometers
        set to np.inf to extrapolate for all points
    compressed: input files are gzip compressed
    scale: scaling factor for converting to output units

OUTPUTS:
    amplitude: amplitudes of tidal constituents
    phase: phases of tidal constituents
    constituents: list of model constituents

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

PROGRAM DEPENDENCIES:
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates

UPDATE HISTORY:
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: reformat arguments to extract_GOT_constants definition
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 12/2021: adjust longitude convention based on model longitude
    Updated 07/2021: added check that tide model files are accessible
    Updated 06/2021: add warning for tide models being entered as string
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: add extrapolation check where there are no invalid points
        prevent ComplexWarning for fill values when calculating amplitudes
        simplified inputs to be similar to binary OTIS read program
        replaced numpy bool/int to prevent deprecation warnings
    Updated 02/2021: set invalid values to nan in extrapolation
        replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added valid data extrapolation with nearest_extrap
    Updated 09/2020: set bounds error to false for regular grid interpolations
        adjust dimensions of input coordinates to be iterable
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Updated 07/2020: added function docstrings. separate bilinear interpolation
        update griddata interpolation. add option for compression
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Updated 11/2019: find invalid mask points for each constituent
    Updated 09/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 07/2019: interpolate fill value mask with bivariate splines
    Updated 12/2018: python3 compatibility updates for division and zip
    Updated 10/2018: added scale as load tides are in mm and ocean are in cm
    Updated 08/2018: added multivariate spline interpolation option
    Written 07/2018
"""
from __future__ import division

import os
import re
import copy
import gzip
import warnings
import numpy as np
import scipy.interpolate
from pyTMD.bilinear_interp import bilinear_interp
from pyTMD.nearest_extrap import nearest_extrap

# PURPOSE: extract tidal harmonic constants out of GOT model at coordinates
def extract_GOT_constants(ilon, ilat, model_files=None, **kwargs):
    """
    Reads files for Richard Ray's Global Ocean Tide (GOT) models

    Makes initial calculations to run the tide program

    Spatially interpolates tidal constituents to input coordinates

    Parameters
    ----------
    ilon: float
        longitude to interpolate
    ilat: float
        latitude to interpolate
    model_files: list or NoneType, default None
        list of model files for each constituent
    method: str, default 'spline'
        Interpolation method

            - ``'bilinear'``: quick bilinear interpolation
            - ``'spline'``: scipy bivariate spline interpolation
            - ``'linear'``, ``'nearest'``: scipy regular grid interpolations
    extrapolate: bool, default False
        Extrapolate model using nearest-neighbors
    cutoff: float, default 10.0
        Extrapolation cutoff in kilometers

        Set to np.inf to extrapolate for all points
    compressed: bool, default False
        Input files are gzip compressed
    scale: float, default 1.0
        Scaling factor for converting to output units

    Returns
    -------
    amplitude: float
        amplitudes of tidal constituents
    phase: float
        phases of tidal constituents
    constituents: list
        list of model constituents
    """
    # set default keyword arguments
    kwargs.setdefault('method', 'spline')
    kwargs.setdefault('extrapolate', False)
    kwargs.setdefault('cutoff', 10.0)
    kwargs.setdefault('compressed', False)
    kwargs.setdefault('scale', 1.0)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(METHOD='method',
        EXTRAPOLATE='extrapolate',CUTOFF='cutoff',
        GZIP='compressed',SCALE='scale')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # raise warning if model files are entered as a string
    if isinstance(model_files,str):
        warnings.warn("Tide model is entered as a string")
        model_files = [model_files]

    # adjust dimensions of input coordinates to be iterable
    ilon = np.atleast_1d(ilon)
    ilat = np.atleast_1d(ilat)
    # number of points
    npts = len(ilon)
    # number of constituents
    nc = len(model_files)
    constituents = []

    # amplitude and phase
    amplitude = np.ma.zeros((npts,nc))
    amplitude.mask = np.zeros((npts,nc),dtype=bool)
    ph = np.ma.zeros((npts,nc))
    ph.mask = np.zeros((npts,nc),dtype=bool)
    # read and interpolate each constituent
    for i,model_file in enumerate(model_files):
        # check that model file is accessible
        if not os.access(os.path.expanduser(model_file), os.F_OK):
            raise FileNotFoundError(os.path.expanduser(model_file))
        # read constituent from elevation file
        hc,lon,lat,cons = read_GOT_grid(os.path.expanduser(model_file),
            compressed=kwargs['compressed'])
        # append to the list of constituents
        constituents.append(cons)
        # adjust longitudinal convention of input latitude and longitude
        # to fit tide model convention
        if (np.min(ilon) < 0.0) & (np.max(lon) > 180.0):
            # input points convention (-180:180)
            # tide model convention (0:360)
            ilon[ilon<0.0] += 360.0
        elif (np.max(ilon) > 180.0) & (np.min(lon) < 0.0):
            # input points convention (0:360)
            # tide model convention (-180:180)
            ilon[ilon>180.0] -= 360.0
        # grid step size of tide model
        dlon = np.abs(lon[1] - lon[0])
        dlat = np.abs(lat[1] - lat[0])
        # replace original values with extend matrices
        lon = extend_array(lon,dlon)
        hc = extend_matrix(hc)
        # interpolated complex form of constituent oscillation
        hci = np.ma.zeros((npts), dtype=hc.dtype, fill_value=hc.fill_value)
        hci.mask = np.zeros((npts),dtype=bool)
        # interpolate amplitude and phase of the constituent
        if (kwargs['method'] == 'bilinear'):
            # replace invalid values with nan
            hc[hc.mask] = np.nan
            # use quick bilinear to interpolate values
            hci.data[:] = bilinear_interp(lon, lat, hc, ilon, ilat,
                dtype=hc.dtype)
            # replace nan values with fill_value
            hci.mask[:] |= np.isnan(hci.data)
            hci.data[hci.mask] = hci.fill_value
        elif (kwargs['method'] == 'spline'):
            # interpolate complex form of the constituent with scipy
            f1=scipy.interpolate.RectBivariateSpline(lon, lat,
                hc.data.real.T, kx=1, ky=1)
            f2=scipy.interpolate.RectBivariateSpline(lon, lat,
                hc.data.imag.T, kx=1, ky=1)
            f3=scipy.interpolate.RectBivariateSpline(lon, lat,
                hc.mask.T, kx=1, ky=1)
            hci.data.real[:] = f1.ev(ilon,ilat)
            hci.data.imag[:] = f2.ev(ilon,ilat)
            hci.mask[:] = f3.ev(ilon,ilat).astype(bool)
            # replace invalid values with fill_value
            hci.data[hci.mask] = hci.fill_value
        else:
            # use scipy regular grid to interpolate values for a given method
            r1 = scipy.interpolate.RegularGridInterpolator((lat,lon),
                hc.data, method=kwargs['method'], bounds_error=False,
                fill_value=hci.fill_value)
            r2 = scipy.interpolate.RegularGridInterpolator((lat,lon),
                hc.mask, method=kwargs['method'], bounds_error=False,
                fill_value=1)
            hci.data[:] = r1.__call__(np.c_[ilat,ilon])
            hci.mask[:] = np.ceil(r2.__call__(np.c_[ilat,ilon])).astype(bool)
            # replace invalid values with fill_value
            hci.mask[:] |= (hci.data == hci.fill_value)
            hci.data[hci.mask] = hci.fill_value
        # extrapolate data using nearest-neighbors
        if kwargs['extrapolate'] and np.any(hci.mask):
            # find invalid data points
            inv, = np.nonzero(hci.mask)
            # replace invalid values with nan
            hc[hc.mask] = np.nan
            # extrapolate points within cutoff of valid model points
            hci[inv] = nearest_extrap(lon, lat, hc, ilon[inv], ilat[inv],
                dtype=hc.dtype, cutoff=kwargs['cutoff'])
        # convert amplitude from input units to meters
        amplitude.data[:,i] = np.abs(hci.data)*kwargs['scale']
        amplitude.mask[:,i] = np.copy(hci.mask)
        # phase of the constituent in radians
        ph.data[:,i] = np.arctan2(-np.imag(hci.data),np.real(hci.data))
        ph.mask[:,i] = np.copy(hci.mask)

    # convert phase to degrees
    phase = ph*180.0/np.pi
    phase.data[phase.data < 0] += 360.0
    # replace data for invalid mask values
    amplitude.data[amplitude.mask] = amplitude.fill_value
    phase.data[phase.mask] = phase.fill_value
    # return the interpolated values
    return (amplitude, phase, constituents)

# PURPOSE: wrapper function to extend an array
def extend_array(input_array, step_size):
    """
    Wrapper function to extend an array

    Parameters
    ----------
    input_array: float
        array to extend
    step_size: float
        step size between elements of array

    Returns
    -------
    temp: float
        extended array
    """
    n = len(input_array)
    temp = np.zeros((n+3),dtype=input_array.dtype)
    # extended array [x-1,x0,...,xN,xN+1,xN+2]
    temp[0] = input_array[0] - step_size
    temp[1:-2] = input_array[:]
    temp[-2] = input_array[-1] + step_size
    temp[-1] = input_array[-1] + 2.0*step_size
    return temp

# PURPOSE: wrapper function to extend a matrix
def extend_matrix(input_matrix):
    """
    Wrapper function to extend a matrix

    Parameters
    ----------
    input_matrix: float
        matrix to extend

    Returns
    -------
    temp: float
        extended matrix
    """
    ny,nx = np.shape(input_matrix)
    temp = np.ma.zeros((ny,nx+3),dtype=input_matrix.dtype)
    temp[:,0] = input_matrix[:,-1]
    temp[:,1:-2] = input_matrix[:,:]
    temp[:,-2] = input_matrix[:,0]
    temp[:,-1] = input_matrix[:,1]
    return temp

# PURPOSE: read GOT model grid files
def read_GOT_grid(input_file, **kwargs):
    """
    Read Richard Ray's Global Ocean Tide (GOT) model file

    Parameters
    ----------
    input_file: str
        Model file
    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    hc: complex
        complex form of tidal constituent oscillation
    lon: float
        longitude of tidal model
    lat: float
        latitude of tidal model
    cons: str
        tidal constituent ID
    """
    # set default keyword arguments
    kwargs.setdefault('compressed', False)
    # tilde-expand input file
    input_file = os.path.expanduser(input_file)
    # read input tide model file
    if kwargs['compressed']:
        # read gzipped ascii file
        with gzip.open(input_file, 'rb') as f:
            file_contents = f.read().decode('utf8').splitlines()
    else:
        with open(input_file, mode="r", encoding='utf8') as f:
            file_contents = f.read().splitlines()
    # parse header text
    constituent_list = ['Q1','O1','P1','K1','N2','M2','S2','K2','S1','M4']
    regex = re.compile(r'|'.join(constituent_list), re.IGNORECASE)
    cons = regex.findall(file_contents[0]).pop().lower()
    nlat,nlon = np.array(file_contents[2].split(), dtype=int)
    # longitude range
    ilat = np.array(file_contents[3].split(), dtype=np.float64)
    # latitude range
    ilon = np.array(file_contents[4].split(), dtype=np.float64)
    # mask fill value
    fill_value = np.array(file_contents[5].split(), dtype=np.float64)
    # create output variables
    lat = np.linspace(ilat[0],ilat[1],nlat)
    lon = np.linspace(ilon[0],ilon[1],nlon)
    amp = np.ma.zeros((nlat,nlon),fill_value=fill_value[0],dtype=np.float32)
    ph = np.ma.zeros((nlat,nlon),fill_value=fill_value[0],dtype=np.float32)
    # create masks for output variables (0=valid)
    amp.mask = np.zeros((nlat,nlon),dtype=bool)
    ph.mask = np.zeros((nlat,nlon),dtype=bool)
    # starting lines to fill amplitude and phase variables
    l1 = 7
    l2 = 14 + int(nlon//11)*nlat + nlat
    # for each latitude
    for i in range(nlat):
        for j in range(nlon//11):
            j1 = j*11
            amp.data[i,j1:j1+11] = np.array(file_contents[l1].split(),dtype='f')
            ph.data[i,j1:j1+11] = np.array(file_contents[l2].split(),dtype='f')
            l1 += 1
            l2 += 1
        # add last tidal variables
        j1 = (j+1)*11; j2 = nlon % 11
        amp.data[i,j1:j1+j2] = np.array(file_contents[l1].split(),dtype='f')
        ph.data[i,j1:j1+j2] = np.array(file_contents[l2].split(),dtype='f')
        l1 += 1
        l2 += 1
    # calculate complex form of constituent oscillation
    hc = amp*np.exp(-1j*ph*np.pi/180.0)
    # set masks
    hc.mask = (amp.data == amp.fill_value) | (ph.data == ph.fill_value)
    # return output variables
    return (hc,lon,lat,cons)
