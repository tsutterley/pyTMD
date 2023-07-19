#!/usr/bin/env python
u"""
GOT.py
Written by Tyler Sutterley (04/2023)

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
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    interpolate.py: interpolation routines for spatial data

UPDATE HISTORY:
    Updated 04/2023: fix repeated longitudinal convention adjustment
        using pathlib to define and expand tide model paths
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: refactor tide read programs under io
        new functions to read and interpolate from constituents class
        new functions to read and write GOT netCDF4 files
        refactored interpolation routines into new module
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
from __future__ import division, annotations

import re
import copy
import gzip
import uuid
import logging
import pathlib
import datetime
import warnings
import numpy as np
import pyTMD.version
import pyTMD.interpolate
import pyTMD.io.constituents

# attempt imports
try:
    import netCDF4
except (ImportError, ModuleNotFoundError) as exc:
    logging.debug("netCDF4 not available")

# PURPOSE: extract harmonic constants from tide models at coordinates
def extract_constants(
        ilon: np.ndarray,
        ilat: np.ndarray,
        model_files: str | pathlib.Path | list | None = None,
        **kwargs
    ):
    """
    Reads files for Richard Ray's Global Ocean Tide (GOT) models

    Makes initial calculations to run the tide program

    Spatially interpolates tidal constituents to input coordinates

    Parameters
    ----------
    ilon: np.ndarray
        longitude to interpolate
    ilat: np.ndarray
        latitude to interpolate
    model_files: str, list, pathlib.Path or NoneType, default None
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

        Set to ``np.inf`` to extrapolate for all points
    grid: str, default 'ascii'
        Tide model file type to read

            - ``'ascii'``: traditional GOT ascii format
            - ``'netcdf'``: GOT netCDF4 format
    compressed: bool, default False
        Input files are gzip compressed
    scale: float, default 1.0
        Scaling factor for converting to output units

    Returns
    -------
    amplitude: np.ndarray
        amplitudes of tidal constituents
    phase: np.ndarray
        phases of tidal constituents
    constituents: np.ndarray
        list of model constituents
    """
    # set default keyword arguments
    kwargs.setdefault('method', 'spline')
    kwargs.setdefault('extrapolate', False)
    kwargs.setdefault('cutoff', 10.0)
    kwargs.setdefault('grid', 'ascii')
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

    # raise warning if model files are entered as a string or path
    if isinstance(model_files, (str, pathlib.Path)):
        warnings.warn("Tide model is entered as a string")
        model_files = [model_files]

    # adjust dimensions of input coordinates to be iterable
    ilon = np.atleast_1d(np.copy(ilon))
    ilat = np.atleast_1d(np.copy(ilat))
    # number of points
    npts = len(ilon)
    # number of constituents
    nc = len(model_files)
    # list of constituents
    constituents = []

    # amplitude and phase
    amplitude = np.ma.zeros((npts,nc))
    amplitude.mask = np.zeros((npts,nc),dtype=bool)
    ph = np.ma.zeros((npts,nc))
    ph.mask = np.zeros((npts,nc),dtype=bool)
    # read and interpolate each constituent
    for i,model_file in enumerate(model_files):
        # check that model file is accessible
        model_file = pathlib.Path(model_file).expanduser()
        if not model_file.exists():
            raise FileNotFoundError(str(model_file))
        # read constituent from elevation file
        if (kwargs['grid'] == 'ascii'):
            hc, lon, lat, cons = read_ascii_file(model_file,
                compressed=kwargs['compressed'])
        elif (kwargs['grid'] == 'netcdf'):
            hc, lon, lat, cons = read_netcdf_file(model_file,
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
        # replace original values with extend matrices
        lon = extend_array(lon,dlon)
        hc = extend_matrix(hc)
        # interpolate amplitude and phase of the constituent
        if (kwargs['method'] == 'bilinear'):
            # replace invalid values with nan
            hc[hc.mask] = np.nan
            # use quick bilinear to interpolate values
            hci = pyTMD.interpolate.bilinear(lon, lat, hc, ilon, ilat,
                dtype=hc.dtype)
            # replace nan values with fill_value
            hci.mask[:] |= np.isnan(hci.data)
            hci.data[hci.mask] = hci.fill_value
        elif (kwargs['method'] == 'spline'):
            # interpolate complex form of the constituent
            # use scipy splines to interpolate values
            hci = pyTMD.interpolate.spline(lon, lat, hc, ilon, ilat,
                dtype=hc.dtype, reducer=np.ceil, kx=1, ky=1)
            # replace invalid values with fill_value
            hci.data[hci.mask] = hci.fill_value
        else:
            # interpolate complex form of the constituent
            # use scipy regular grid to interpolate values
            hci = pyTMD.interpolate.regulargrid(lon, lat, hc, ilon, ilat,
                dtype=hc.dtype, method=kwargs['method'], reducer=np.ceil,
                bounds_error=False)
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
            hci[inv] = pyTMD.interpolate.extrapolate(lon, lat, hc,
                ilon[inv], ilat[inv], dtype=hc.dtype,
                cutoff=kwargs['cutoff'])
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

# PURPOSE: read harmonic constants from tide models
def read_constants(
        model_files: str | pathlib.Path | list | None = None,
        **kwargs
    ):
    """
    Reads files for Richard Ray's Global Ocean Tide (GOT) models

    Parameters
    ----------
    model_files: str, list, pathlib.Path or NoneType, default None
        list of model files for each constituent
    grid: str, default 'ascii'
        Tide model file type to read

            - ``'ascii'``: traditional GOT ascii format
            - ``'netcdf'``: GOT netCDF4 format
    compressed: bool, default False
        Input files are gzip compressed

    Returns
    -------
    constituents: obj
        complex form of tide model constituents
    """
    # set default keyword arguments
    kwargs.setdefault('grid', 'ascii')
    kwargs.setdefault('compressed', False)

    # raise warning if model files are entered as a string
    if isinstance(model_files, (str, pathlib.Path)):
        warnings.warn("Tide model is entered as a string")
        model_files = [model_files]

    # save output constituents
    constituents = pyTMD.io.constituents()
    # read each model constituent
    for i, model_file in enumerate(model_files):
        # check that model file is accessible
        model_file = pathlib.Path(model_file).expanduser()
        if not model_file.exists():
            raise FileNotFoundError(str(model_file))
        # read constituent from elevation file
        if (kwargs['grid'] == 'ascii'):
            hc, lon, lat, cons = read_ascii_file(model_file,
                compressed=kwargs['compressed'])
        elif (kwargs['grid'] == 'netcdf'):
            hc, lon, lat, cons = read_netcdf_file(model_file,
                compressed=kwargs['compressed'])
        # grid step size of tide model
        dlon = np.abs(lon[1] - lon[0])
        # replace original values with extend matrices
        lon = extend_array(lon, dlon)
        hc = extend_matrix(hc)
        # append extended constituent
        constituents.append(cons, hc)
        # set model coordinates
        setattr(constituents, 'longitude', lon)
        setattr(constituents, 'latitude', lat)

    # return the complex form of the model constituents
    return constituents

# PURPOSE: interpolate constants from tide models to input coordinates
def interpolate_constants(
        ilon: np.ndarray,
        ilat: np.ndarray,
        constituents,
        **kwargs
    ):
    """
    Interpolate constants from GOT tidal models to input coordinates

    Makes initial calculations to run the tide program

    Parameters
    ----------
    ilon: np.ndarray
        longitude to interpolate
    ilat: np.ndarray
        latitude to interpolate
    constituents: obj
        Tide model constituents (complex form)
    method: str, default 'spline'
        Interpolation method

            - ``'bilinear'``: quick bilinear interpolation
            - ``'spline'``: scipy bivariate spline interpolation
            - ``'linear'``, ``'nearest'``: scipy regular grid interpolations
    extrapolate: bool, default False
        Extrapolate model using nearest-neighbors
    cutoff: float, default 10.0
        Extrapolation cutoff in kilometers

        Set to ``np.inf`` to extrapolate for all points
    scale: float, default 1.0
        Scaling factor for converting to output units

    Returns
    -------
    amplitude: np.ndarray
        amplitudes of tidal constituents
    phase: np.ndarray
        phases of tidal constituents
    """
    # set default keyword arguments
    kwargs.setdefault('method', 'spline')
    kwargs.setdefault('extrapolate', False)
    kwargs.setdefault('cutoff', 10.0)
    kwargs.setdefault('scale', 1.0)

    # verify that constituents are valid class instance
    assert isinstance(constituents, pyTMD.io.constituents)
    # extract model coordinates
    lon = np.copy(constituents.longitude)
    lat = np.copy(constituents.latitude)

    # adjust dimensions of input coordinates to be iterable
    ilon = np.atleast_1d(np.copy(ilon))
    ilat = np.atleast_1d(np.copy(ilat))
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
    # number of points
    npts = len(ilon)
    # number of constituents
    nc = len(constituents)

    # amplitude and phase
    amplitude = np.ma.zeros((npts,nc))
    amplitude.mask = np.zeros((npts,nc),dtype=bool)
    ph = np.ma.zeros((npts,nc))
    ph.mask = np.zeros((npts,nc),dtype=bool)
    # default complex fill value
    fill_value = np.ma.default_fill_value(np.dtype(complex))
    # interpolate each constituent
    for i, c in enumerate(constituents.fields):
        # get model constituent
        hc = constituents.get(c)
        # interpolate amplitude and phase of the constituent
        if (kwargs['method'] == 'bilinear'):
            # replace invalid values with nan
            hc.data[hc.mask] = np.nan
            # use quick bilinear to interpolate values
            hci = pyTMD.interpolate.bilinear(lon, lat, hc, ilon, ilat,
                fill_value=fill_value,
                dtype=hc.dtype)
            # replace nan values with fill_value
            hci.mask[:] |= np.isnan(hci.data)
            hci.data[hci.mask] = hci.fill_value
        elif (kwargs['method'] == 'spline'):
            # replace invalid values with fill value
            hc.data[hc.mask] = fill_value
            # interpolate complex form of the constituent
            # use scipy splines to interpolate values
            hci = pyTMD.interpolate.spline(lon, lat, hc, ilon, ilat,
                fill_value=fill_value,
                dtype=hc.dtype,
                reducer=np.ceil,
                kx=1, ky=1)
            # replace invalid values with fill_value
            hci.data[hci.mask] = hci.fill_value
        else:
            # replace invalid values with fill value
            hc.data[hc.mask] = fill_value
            # interpolate complex form of the constituent
            # use scipy regular grid to interpolate values
            hci = pyTMD.interpolate.regulargrid(lon, lat, hc, ilon, ilat,
                fill_value=fill_value,
                dtype=hc.dtype,
                method=kwargs['method'],
                reducer=np.ceil,
                bounds_error=False)
            # replace invalid values with fill_value
            hci.mask[:] |= (hci.data == hci.fill_value)
            hci.data[hci.mask] = hci.fill_value
        # extrapolate data using nearest-neighbors
        if kwargs['extrapolate'] and np.any(hci.mask):
            # find invalid data points
            inv, = np.nonzero(hci.mask)
            # replace invalid values with nan
            hc.data[hc.mask] = np.nan
            # extrapolate points within cutoff of valid model points
            hci[inv] = pyTMD.interpolate.extrapolate(lon, lat, hc,
                ilon[inv], ilat[inv], dtype=hc.dtype,
                cutoff=kwargs['cutoff'])
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
    return (amplitude, phase)

# PURPOSE: Extend a longitude array
def extend_array(input_array: np.ndarray, step_size: float):
    """
    Extends a longitude array

    Parameters
    ----------
    input_array: np.ndarray
        array to extend
    step_size: float
        step size between elements of array

    Returns
    -------
    temp: np.ndarray
        extended array
    """
    n = len(input_array)
    temp = np.zeros((n+3), dtype=input_array.dtype)
    # extended array [x-1,x0,...,xN,xN+1,xN+2]
    temp[0] = input_array[0] - step_size
    temp[1:-2] = input_array[:]
    temp[-2] = input_array[-1] + step_size
    temp[-1] = input_array[-1] + 2.0*step_size
    return temp

# PURPOSE: Extend a global matrix
def extend_matrix(input_matrix: np.ndarray):
    """
    Extends a global matrix

    Parameters
    ----------
    input_matrix: np.ndarray
        matrix to extend

    Returns
    -------
    temp: np.ndarray
        extended matrix
    """
    ny, nx = np.shape(input_matrix)
    temp = np.ma.zeros((ny,nx+3), dtype=input_matrix.dtype)
    temp[:,0] = input_matrix[:,-1]
    temp[:,1:-2] = input_matrix[:,:]
    temp[:,-2] = input_matrix[:,0]
    temp[:,-1] = input_matrix[:,1]
    return temp

# PURPOSE: read GOT model grid files
def read_ascii_file(
        input_file: str | pathlib.Path,
        **kwargs
    ):
    """
    Read Richard Ray's Global Ocean Tide (GOT) model file

    Parameters
    ----------
    input_file: str or pathlib.Path
        Model file
    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    hc: np.ndarray
        complex form of tidal constituent oscillation
    lon: np.ndarray
        longitude of tidal model
    lat: np.ndarray
        latitude of tidal model
    cons: str
        tidal constituent ID
    """
    # set default keyword arguments
    kwargs.setdefault('compressed', False)
    # tilde-expand input file
    input_file = pathlib.Path(input_file).expanduser()
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
    return (hc, lon, lat, cons)

# PURPOSE: read GOT netCDF4 tide model files
def read_netcdf_file(
        input_file: str | pathlib.Path,
        **kwargs
    ):
    """
    Read Richard Ray's Global Ocean Tide (GOT) netCDF4 model file

    Parameters
    ----------
    input_file: str or pathlib.Path
        model file
    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    hc: np.ndarray
        complex form of tidal constituent oscillation
    lon: np.ndarray
        longitude of tidal model
    lat: np.ndarray
        latitude of tidal model
    cons: str
        tidal constituent ID
    """
    # set default keyword arguments
    kwargs.setdefault('compressed', False)
    # tilde-expand input file
    input_file = pathlib.Path(input_file).expanduser()
    # read the netcdf format tide elevation file
    if kwargs['compressed']:
        # read gzipped netCDF4 file
        f = gzip.open(input_file, 'rb')
        fileID = netCDF4.Dataset(uuid.uuid4().hex, 'r', memory=f.read())
    else:
        fileID = netCDF4.Dataset(input_file, 'r')
    # variable dimensions
    lon = fileID.variables['longitude'][:]
    lat = fileID.variables['latitude'][:]
    # get amplitude and phase components
    amp = fileID.variables['amplitude'][:]
    ph = fileID.variables['phase'][:]
    # extract constituent from attribute
    cons = fileID.Constituent.lower()
    # close the file
    fileID.close()
    f.close() if kwargs['compressed'] else None
    # calculate complex form of constituent oscillation
    mask = (amp.data == amp.fill_value) | \
        (ph.data == ph.fill_value) | \
        np.isnan(amp.data) | np.isnan(ph.data)
    hc = np.ma.array(amp*np.exp(-1j*ph*np.pi/180.0), mask=mask,
        fill_value=np.ma.default_fill_value(np.dtype(complex)))
    # return output variables
    return (hc, lon, lat, cons)

# PURPOSE: output tidal constituent file in GOT netCDF format
def output_netcdf_file(
        FILE: str | pathlib.Path,
        hc: np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        constituent: str
    ):
    """
    Writes tidal constituent files in GOT netCDF format

    Parameters
    ----------
    FILE: str or pathlib.Path
        output GOT model file name
    hc: np.ndarray
        Eulerian form of tidal constituent
    lon: np.ndarray
        longitude coordinates
    lat: np.ndarray
        latitude coordinates
    constituent: str
        tidal constituent ID
    """
    # tilde-expand output file
    FILE = pathlib.Path(FILE).expanduser()
    # opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILE, 'w', format="NETCDF4")
    # define the NetCDF dimensions
    fileID.createDimension('longitude', len(lon))
    fileID.createDimension('latitude', len(lat))
    # calculate amplitude and phase
    amp = np.abs(hc)
    ph = 180.0*np.arctan2(-np.imag(hc), np.real(hc))/np.pi
    ph.data[ph.data < 0] += 360.0
    # update masks and fill values
    amp.mask = np.copy(hc.mask)
    amp.data[amp.mask] = amp.fill_value
    ph.mask = np.copy(hc.mask)
    ph.data[ph.mask] = ph.fill_value
    # defining the NetCDF variables
    nc = {}
    nc['longitude'] = fileID.createVariable('longitude', lon.dtype,
        ('longitude',))
    nc['latitude'] = fileID.createVariable('latitude', lat.dtype,
        ('latitude',))
    nc['amplitude'] = fileID.createVariable('amplitude', amp.dtype,
        ('latitude','longitude',), fill_value=amp.fill_value, zlib=True)
    nc['phase'] = fileID.createVariable('phase', ph.dtype,
        ('latitude','longitude',), fill_value=ph.fill_value, zlib=True)
    # filling the NetCDF variables
    nc['longitude'][:] = lon[:]
    nc['latitude'][:] = lat[:]
    nc['amplitude'][:] = amp[:]
    nc['phase'][:] = ph[:]
    # set variable attributes for coordinates
    nc['longitude'].setncattr('units', 'degrees_east')
    nc['longitude'].setncattr('long_name', 'longitude')
    nc['latitude'].setncattr('units', 'degrees_north')
    nc['latitude'].setncattr('long_name', 'latitude')
    # set variable attributes
    nc['amplitude'].setncattr('units', 'cm')
    nc['amplitude'].setncattr('long_name', 'Tide amplitude')
    nc['phase'].setncattr('units', 'degrees')
    nc['phase'].setncattr('long_name', 'Greenwich tide phase lag')
    # add global attributes
    fileID.title = 'GOT tide file'
    fileID.authors = 'Richard Ray'
    fileID.institution = 'NASA Goddard Space Flight Center'
    # add attribute for tidal constituent ID
    fileID.Constituent = constituent.upper()
    # add attribute for date created
    fileID.date_created = datetime.datetime.now().isoformat()
    # add attributes for software information
    fileID.software_reference = pyTMD.version.project_name
    fileID.software_version = pyTMD.version.full_version
    # Output NetCDF structure information
    logging.info(str(FILE))
    logging.info(list(fileID.variables.keys()))
    # Closing the NetCDF file
    fileID.close()
