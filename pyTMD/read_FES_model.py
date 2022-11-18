#!/usr/bin/env python
u"""
read_FES_model.py (11/2022)
Reads files for a tidal model and makes initial calculations to run tide program
Includes functions to extract tidal harmonic constants from the
    FES (Finite Element Solution) tide models for given locations
ascii and netCDF4 files can be been compressed using gzip

Reads ascii and netCDF4 FES tidal solutions provided by AVISO
    https://www.aviso.altimetry.fr/data/products/auxiliary-products/
        global-tide-fes.html

INPUTS:
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    model_files: list of model files for each constituent

OPTIONS:
    type: tidal variable to run
        z: heights
        u: horizontal transport velocities
        v: vertical transport velocities
    version: model version to run
        FES1999
        FES2004
        FES2012
        FES2014
        EOT20
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

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates

UPDATE HISTORY:
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 05/2022: reformat arguments to extract_FES_constants definition
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
        fix netCDF4 masks for nan values
    Updated 01/2022: added global Empirical Ocean Tide model (EOT20)
    Updated 12/2021: adjust longitude convention based on model longitude
    Updated 07/2021: added check that tide model files are accessible
    Updated 06/2021: add warning for tide models being entered as string
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: add extrapolation check where there are no invalid points
        prevent ComplexWarning for fill values when calculating amplitudes
        simplified inputs to be similar to binary OTIS read program
        replaced numpy bool/int to prevent deprecation warnings
        use uuid for reading from gzipped netCDF4 files
    Updated 02/2021: set invalid values to nan in extrapolation
        replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added nearest-neighbor data extrapolation
    Updated 09/2020: set bounds error to false for regular grid interpolations
        adjust dimensions of input coordinates to be iterable
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Written 07/2020
"""
import os
import copy
import gzip
import uuid
import warnings
import numpy as np
import scipy.interpolate
from pyTMD.bilinear_interp import bilinear_interp
from pyTMD.nearest_extrap import nearest_extrap

# attempt imports
try:
    import netCDF4
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("netCDF4 not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: extract tidal harmonic constants from tide models at coordinates
def extract_FES_constants(ilon, ilat, model_files=None, **kwargs):
    """
    Reads files for an ascii or netCDF4 tidal model

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
    type: str, default 'z'
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'v'``: vertical transport velocities
    version: str or NoneType, default None
        Model version to read

            - ``'FES1999'``
            - ``'FES2004'``
            - ``'FES2012'``
            - ``'FES2014'``
            - ``'EOT20'``
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
    """
    # set default keyword arguments
    kwargs.setdefault('type', 'z')
    kwargs.setdefault('version', None)
    kwargs.setdefault('method', 'spline')
    kwargs.setdefault('extrapolate', False)
    kwargs.setdefault('cutoff', 10.0)
    kwargs.setdefault('compressed', False)
    kwargs.setdefault('scale', 1.0)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(TYPE='type',VERSION='version',
        METHOD='method',EXTRAPOLATE='extrapolate',CUTOFF='cutoff',
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

    # amplitude and phase
    amplitude = np.ma.zeros((npts,nc))
    amplitude.mask = np.zeros((npts,nc),dtype=bool)
    ph = np.ma.zeros((npts,nc))
    ph.mask = np.zeros((npts,nc),dtype=bool)
    # read and interpolate each constituent
    for i,fi in enumerate(model_files):
        # check that model file is accessible
        if not os.access(os.path.expanduser(fi), os.F_OK):
            raise FileNotFoundError(os.path.expanduser(fi))
        # read constituent from elevation file
        if kwargs['version'] in ('FES1999','FES2004'):
            # FES ascii constituent files
            hc,lon,lat = read_ascii_file(os.path.expanduser(fi), **kwargs)
        elif kwargs['version'] in ('FES2012','FES2014','EOT20'):
            # FES netCDF4 constituent files
            hc,lon,lat = read_netcdf_file(os.path.expanduser(fi), **kwargs)
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
    return (amplitude,phase)

# PURPOSE: read FES ascii tide model grid files
def read_ascii_file(input_file, **kwargs):
    """
    Read FES (Finite Element Solution) tide model file

    Parameters
    ----------
    input_file: str
        model file
    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    hc: complex form of tidal constituent oscillation
    lon: longitude of tidal model
    lat: latitude of tidal model
    """
    # set default keyword arguments
    kwargs.setdefault('compressed', False)
    # tilde-expand input file
    input_file = os.path.expanduser(input_file)
    # read input tide model file
    if kwargs['compressed']:
        # read gzipped ascii file
        with gzip.open(input_file, 'rb') as f:
            file_contents = f.read(input_file).splitlines()
    else:
        with open(input_file, mode="r", encoding='utf8') as f:
            file_contents = f.read().splitlines()
    # parse header text
    # longitude range (lonmin, lonmax)
    lonmin,lonmax = np.array(file_contents[0].split(), dtype=np.float64)
    # latitude range (latmin, latmax)
    latmin,latmax = np.array(file_contents[1].split(), dtype=np.float64)
    # grid step size (dlon, dlat)
    dlon,dlat = np.array(file_contents[2].split(), dtype=np.float64)
    # grid dimensions (nlon, nlat)
    nlon,nlat = np.array(file_contents[3].split(), dtype=int)
    # mask fill value
    masked_values = file_contents[4].split()
    fill_value = np.float64(masked_values[0])
    # create output variables
    lat = np.linspace(latmin, latmax, nlat)
    lon = np.linspace(lonmin,lonmax,nlon)
    amp = np.ma.zeros((nlat,nlon),fill_value=fill_value,dtype=np.float32)
    ph = np.ma.zeros((nlat,nlon),fill_value=fill_value,dtype=np.float32)
    # create masks for output variables (0=valid)
    amp.mask = np.zeros((nlat,nlon),dtype=bool)
    ph.mask = np.zeros((nlat,nlon),dtype=bool)
    # starting line to fill amplitude and phase variables
    i1 = 5
    # for each latitude
    for i in range(nlat):
        for j in range(nlon//30):
            j1 = j*30
            amp.data[i,j1:j1+30]=np.array(file_contents[i1].split(),dtype='f')
            ph.data[i,j1:j1+30]=np.array(file_contents[i1+1].split(),dtype='f')
            i1 += 2
        # add last tidal variables
        j1 = (j+1)*30
        j2 = nlon % 30
        amp.data[i,j1:j1+j2] = np.array(file_contents[i1].split(),dtype='f')
        ph.data[i,j1:j1+j2] = np.array(file_contents[i1+1].split(),dtype='f')
        i1 += 2
    # calculate complex form of constituent oscillation
    hc = amp*np.exp(-1j*ph*np.pi/180.0)
    # set masks
    hc.mask = (amp.data == amp.fill_value) | (ph.data == ph.fill_value)
    # return output variables
    return (hc,lon,lat)

# PURPOSE: read FES netCDF4 tide model files
def read_netcdf_file(input_file, **kwargs):
    """
    Read FES (Finite Element Solution) tide model netCDF4 file

    Parameters
    ----------
    input_file: str
        model file
    type: str or NoneType, default None
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'v'``: vertical transport velocities
    version: str or NoneType, default None
        FES model version
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
    """
    # set default keyword arguments
    kwargs.setdefault('type', None)
    kwargs.setdefault('version', None)
    kwargs.setdefault('compressed', False)
    # read the netcdf format tide elevation file
    if kwargs['compressed']:
        # read gzipped netCDF4 file
        f = gzip.open(os.path.expanduser(input_file),'rb')
        fileID = netCDF4.Dataset(uuid.uuid4().hex,'r',memory=f.read())
    else:
        fileID = netCDF4.Dataset(os.path.expanduser(input_file), 'r')
    # variable dimensions for each model
    if kwargs['version'] in ('FES2012',):
        lon = fileID.variables['longitude'][:]
        lat = fileID.variables['latitude'][:]
    elif kwargs['version'] in ('FES2014','EOT20'):
        lon = fileID.variables['lon'][:]
        lat = fileID.variables['lat'][:]
    # amplitude and phase components for each type
    if (kwargs['type'] == 'z'):
        amp = fileID.variables['amplitude'][:]
        ph = fileID.variables['phase'][:]
    elif (kwargs['type'] == 'u'):
        amp = fileID.variables['Ua'][:]
        ph = fileID.variables['Ug'][:]
    elif (kwargs['type'] == 'v'):
        amp = fileID.variables['Va'][:]
        ph = fileID.variables['Vg'][:]
    # close the file
    fileID.close()
    f.close() if kwargs['compressed'] else None
    # calculate complex form of constituent oscillation
    hc = amp*np.exp(-1j*ph*np.pi/180.0)
    # set masks
    hc.mask = (amp.data == amp.fill_value) | \
        (ph.data == ph.fill_value) | \
        np.isnan(amp.data) | np.isnan(ph.data)
    # return output variables
    return (hc,lon,lat)
