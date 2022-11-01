#!/usr/bin/env python
u"""
read_netcdf_model.py (11/2022)
Reads files for a tidal model and makes initial calculations to run tide program
Includes functions to extract tidal harmonic constants from OTIS tide models for
    given locations
netCDF4 files can be been compressed using gzip

Reads netCDF4 ATLAS tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/

INPUTS:
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    grid_file: grid file for model (can be gzipped)
    model_files: list of model files for each constituent (can be gzipped)

OPTIONS:
    type: tidal variable to run
        z: heights
        u: horizontal transport velocities
        U: horizontal depth-averaged transport
        v: vertical transport velocities
        V: vertical depth-averaged transport
    method: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations
    extrapoalte: extrapolate model using nearest-neighbors
    cutoff: extrapolation cutoff in kilometers
        set to np.inf to extrapolate for all points
    compressed: input netCDF4 files are gzip compressed
    scale: scaling factor for converting to output units

OUTPUTS:
    amplitude: amplitudes of tidal constituents
    phase: phases of tidal constituents
    D: bathymetry of tide model
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
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates

UPDATE HISTORY:
    Updated 11/2022: place some imports within try/except statements
    Updated 07/2022: fix setting of masked array data to NaN
    Updated 05/2022: reformat arguments to extract_netcdf_constants definition
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 12/2021: adjust longitude convention based on model longitude
    Updated 09/2021: fix cases where there is no mask on constituent files
    Updated 07/2021: added check that tide model files are accessible
    Updated 06/2021: add warning for tide models being entered as string
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: add extrapolation check where there are no invalid points
        prevent ComplexWarning for fill values when calculating amplitudes
        simplified inputs to be similar to binary OTIS read program
    Updated 02/2021: set invalid values to nan in extrapolation
        replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added valid data extrapolation with nearest_extrap
        replace tostring with tobytes to fix DeprecationWarning
    Updated 11/2020: create function to read bathymetry and spatial coordinates
    Updated 09/2020: set bounds error to false for regular grid interpolations
        adjust dimensions of input coordinates to be iterable
        reduce number of interpolations by copying bathymetry mask to variables
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Updated 07/2020: added function docstrings. separate bilinear interpolation
        changed TYPE variable to keyword argument. update griddata interpolation
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Written 09/2019
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
def extract_netcdf_constants(ilon, ilat,
    grid_file=None,
    model_files=None,
    **kwargs):
    """
    Reads files for ATLAS netCDF4 tidal models

    Makes initial calculations to run the tide program

    Spatially interpolates tidal constituents to input coordinates

    Parameters
    ----------
    ilon: float
        longitude to interpolate
    ilat: float
        latitude to interpolate
    grid_file: str or NoneType, default None
        grid file for model
    model_files: list or NoneType, default None
        list of model files for each constituent
    type: str, default 'z'
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'U'``: horizontal depth-averaged transport
            - ``'v'``: vertical transport velocities
            - ``'V'``: vertical depth-averaged transport
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
    D: float
        bathymetry of tide model
    constituents: list
        list of model constituents
    """
    # set default keyword arguments
    kwargs.setdefault('type', 'z')
    kwargs.setdefault('method', 'spline')
    kwargs.setdefault('extrapolate', False)
    kwargs.setdefault('cutoff', 10.0)
    kwargs.setdefault('compressed', True)
    kwargs.setdefault('scale', 1.0)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(TYPE='type', METHOD='method',
        EXTRAPOLATE='extrapolate', CUTOFF='cutoff',
        GZIP='compressed', SCALE='scale')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn("""Deprecated keyword argument {0}.
                Changed to '{1}'""".format(old,new),
                DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # raise warning if model files are entered as a string
    if isinstance(model_files,str):
        warnings.warn("Tide model is entered as a string")
        model_files = [model_files]

    # check that grid file is accessible
    if not os.access(os.path.expanduser(grid_file), os.F_OK):
        raise FileNotFoundError(os.path.expanduser(grid_file))

    # read the tide grid file for bathymetry and spatial coordinates
    lon, lat, bathymetry = read_netcdf_grid(grid_file, kwargs['type'],
        compressed=kwargs['compressed'])

    # adjust dimensions of input coordinates to be iterable
    ilon = np.atleast_1d(ilon)
    ilat = np.atleast_1d(ilat)
    # adjust longitudinal convention of input latitude and longitude
    # to fit tide model convention
    if (np.min(ilon) < 0.0) & (np.max(lon) > 180.0):
        # input points convention (-180:180)
        # tide model convention (0:360)
        ilon[ilon < 0.0] += 360.0
    elif (np.max(ilon) > 180.0) & (np.min(lon) < 0.0):
        # input points convention (0:360)
        # tide model convention (-180:180)
        ilon[ilon > 180.0] -= 360.0

    # grid step size of tide model
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]
    # replace original values with extend arrays/matrices
    lon = extend_array(lon, dlon)
    bathymetry = extend_matrix(bathymetry)
    # create masks
    bathymetry.mask = (bathymetry.data == 0)

    # number of points
    npts = len(ilon)
    # interpolate bathymetry and mask to output points
    D = np.ma.zeros((npts))
    D.mask = np.zeros((npts), dtype=bool)
    if (kwargs['method'] == 'bilinear'):
        # replace invalid values with nan
        bathymetry.data[bathymetry.mask] = np.nan
        # use quick bilinear to interpolate values
        D.data[:] = bilinear_interp(lon, lat, bathymetry, ilon, ilat)
        # replace nan values with fill_value
        D.mask[:] = np.isnan(D.data)
        D.data[D.mask] = D.fill_value
    elif (kwargs['method'] == 'spline'):
        # use scipy bivariate splines to interpolate values
        f1 = scipy.interpolate.RectBivariateSpline(lon, lat,
            bathymetry.data.T, kx=1, ky=1)
        f2 = scipy.interpolate.RectBivariateSpline(lon, lat,
            bathymetry.mask.T, kx=1, ky=1)
        D.data[:] = f1.ev(ilon,ilat)
        D.mask[:] = np.ceil(f2.ev(ilon,ilat)).astype(bool)
    else:
        # use scipy regular grid to interpolate values for a given method
        r1 = scipy.interpolate.RegularGridInterpolator((lat, lon),
            bathymetry.data, method=kwargs['method'], bounds_error=False)
        r2 = scipy.interpolate.RegularGridInterpolator((lat, lon),
            bathymetry.mask, method=kwargs['method'], bounds_error=False,
            fill_value=1)
        D.data[:] = r1.__call__(np.c_[ilat, ilon])
        D.mask[:] = np.ceil(r2.__call__(np.c_[ilat, ilon])).astype(bool)

    # u and v are velocities in cm/s
    if kwargs['type'] in ('v','u'):
        unit_conv = (D.data/100.0)
    # U and V are transports in m^2/s
    elif kwargs['type'] in ('V','U'):
        unit_conv = 1.0

    # number of constituents
    nc = len(model_files)
    # list of constituents
    constituents = []
    # amplitude and phase
    ampl = np.ma.zeros((npts, nc))
    ampl.mask = np.zeros((npts, nc), dtype=bool)
    ph = np.ma.zeros((npts, nc))
    ph.mask = np.zeros((npts, nc), dtype=bool)
    # read and interpolate each constituent
    for i,model_file in enumerate(model_files):
        # check that model file is accessible
        if not os.access(os.path.expanduser(model_file), os.F_OK):
            raise FileNotFoundError(os.path.expanduser(model_file))
        if (kwargs['type'] == 'z'):
            # read constituent from elevation file
            z,con = read_elevation_file(model_file,
                compressed=kwargs['compressed'])
            # append constituent to list
            constituents.append(con)
            # replace original values with extend matrices
            z = extend_matrix(z)
            # update constituent mask with bathymetry mask
            z.mask[:] |= bathymetry.mask[:]
            # interpolate amplitude and phase of the constituent
            z1 = np.ma.zeros((npts), dtype=z.dtype)
            z1.mask = np.zeros((npts), dtype=bool)
            if (kwargs['method'] == 'bilinear'):
                # replace invalid values with nan
                z.data[z.mask] = np.nan
                z1.data[:] = bilinear_interp(lon, lat, z, ilon, ilat,
                    dtype=z.dtype)
                # mask invalid values
                z1.mask[:] |= np.copy(D.mask)
                z1.data[z1.mask] = z1.fill_value
            elif (kwargs['method'] == 'spline'):
                f1 = scipy.interpolate.RectBivariateSpline(lon, lat,
                    z.data.real.T, kx=1, ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(lon, lat,
                    z.data.imag.T, kx=1, ky=1)
                z1.data.real = f1.ev(ilon,ilat)
                z1.data.imag = f2.ev(ilon,ilat)
                # mask invalid values
                z1.mask[:] |= np.copy(D.mask)
                z1.data[z1.mask] = z1.fill_value
            else:
                # use scipy regular grid to interpolate values
                r1 = scipy.interpolate.RegularGridInterpolator((lat, lon),
                    z.data, method=kwargs['method'], bounds_error=False,
                    fill_value=z1.fill_value)
                z1.data[:] = r1.__call__(np.c_[ilat, ilon])
                # mask invalid values
                z1.mask[:] |= np.copy(D.mask)
                z1.data[z1.mask] = z1.fill_value
            # extrapolate data using nearest-neighbors
            if kwargs['extrapolate'] and np.any(z1.mask):
                # find invalid data points
                inv, = np.nonzero(z1.mask)
                # replace invalid values with nan
                z.data[z.mask] = np.nan
                # extrapolate points within cutoff of valid model points
                z1[inv] = nearest_extrap(lon, lat, z, ilon[inv], ilat[inv],
                    dtype=z.dtype, cutoff=kwargs['cutoff'])
            # amplitude and phase of the constituent
            ampl.data[:,i] = np.abs(z1.data)
            ampl.mask[:,i] = np.copy(z1.mask)
            ph.data[:,i] = np.arctan2(-np.imag(z1.data), np.real(z1.data))
            ph.mask[:,i] = np.copy(z1.mask)
        elif kwargs['type'] in ('U','u','V','v'):
            # read constituent from transport file
            tr,con = read_transport_file(model_file, kwargs['type'],
                compressed=kwargs['compressed'])
            # append constituent to list
            constituents.append(con)
            # replace original values with extend matrices
            tr = extend_matrix(tr)
            # update constituent mask with bathymetry mask
            tr.mask[:] |= bathymetry.mask[:]
            # interpolate amplitude and phase of the constituent
            tr1 = np.ma.zeros((npts), dtype=tr.dtype)
            tr1.mask = np.zeros((npts), dtype=bool)
            if (kwargs['method'] == 'bilinear'):
                # replace invalid values with nan
                tr.data[tr.mask] = np.nan
                tr1.data[:]=bilinear_interp(lon, lat, tr, ilon, ilat,
                    dtype=tr.dtype)
                # mask invalid values
                tr1.mask[:] |= np.copy(D.mask)
                tr1.data[tr1.mask] = tr1.fill_value
            elif (kwargs['method'] == 'spline'):
                f1 = scipy.interpolate.RectBivariateSpline(lon, lat,
                    tr.data.real.T, kx=1, ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(lon, lat,
                    tr.data.imag.T, kx=1, ky=1)
                tr1.data.real = f1.ev(ilon,ilat)
                tr1.data.imag = f2.ev(ilon,ilat)
                # mask invalid values
                tr1.mask[:] |= np.copy(D.mask)
                tr1.data[tr1.mask] = z1.fill_value
            else:
                # use scipy regular grid to interpolate values
                r1 = scipy.interpolate.RegularGridInterpolator((lat, lon),
                    tr.data, method=kwargs['method'], bounds_error=False,
                    fill_value=tr1.fill_value)
                tr1.data[:] = r1.__call__(np.c_[ilat, ilon])
                # mask invalid values
                tr1.mask[:] |= np.copy(D.mask)
                tr1.data[tr1.mask] = tr1.fill_value
            # extrapolate data using nearest-neighbors
            if kwargs['extrapolate'] and np.any(tr1.mask):
                # find invalid data points
                inv, = np.nonzero(tr1.mask)
                # replace invalid values with nan
                tr.data[tr.mask] = np.nan
                # extrapolate points within cutoff of valid model points
                tr1[inv] = nearest_extrap(lon, lat, tr, ilon[inv], ilat[inv],
                    dtype=tr.dtype, cutoff=kwargs['cutoff'])
            # convert units
            # amplitude and phase of the constituent
            ampl.data[:,i] = np.abs(tr1.data)/unit_conv
            ampl.mask[:,i] = np.copy(tr1.mask)
            ph.data[:,i] = np.arctan2(-np.imag(tr1.data), np.real(tr1.data))
            ph.mask[:,i] = np.copy(tr1.mask)

    # convert amplitude from input units to meters
    amplitude = ampl*kwargs['scale']
    # convert phase to degrees
    phase = ph*180.0/np.pi
    phase[phase < 0] += 360.0
    # return the interpolated values
    return (amplitude, phase, D, constituents)

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
    temp = np.zeros((n+2), dtype=input_array.dtype)
    # extended array [x-1,x0,...,xN,xN+1]
    temp[0] = input_array[0] - step_size
    temp[1:-1] = input_array[:]
    temp[-1] = input_array[-1] + step_size
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
    temp = np.ma.zeros((ny,nx+2), dtype=input_matrix.dtype)
    temp[:,0] = input_matrix[:,-1]
    temp[:,1:-1] = input_matrix[:,:]
    temp[:,-1] = input_matrix[:,0]
    return temp

# PURPOSE: read grid file
def read_netcdf_grid(input_file, variable, **kwargs):
    """
    Read grid file to extract model coordinates and bathymetry

    Parameters
    ----------
    input_file: str
        input grid file
    variable: str
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'U'``: horizontal depth-averaged transport
            - ``'v'``: vertical transport velocities
            - ``'V'``: vertical depth-averaged transport

    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    lon: float
        longitudinal coordinates of input grid
    lat: float
        latitudinal coordinates of input grid
    bathymetry: float
        model bathymetry
    """
    # set default keyword arguments
    kwargs.setdefault('compressed', False)
    # read the netcdf format tide grid file
    # reading a combined global solution with localized solutions
    if kwargs['compressed']:
        # read gzipped netCDF4 file
        f = gzip.open(os.path.expanduser(input_file), 'rb')
        fileID=netCDF4.Dataset(uuid.uuid4().hex, 'r', memory=f.read())
    else:
        fileID=netCDF4.Dataset(os.path.expanduser(input_file), 'r')
    # variable dimensions
    nx = fileID.dimensions['nx'].size
    ny = fileID.dimensions['ny'].size
    # allocate numpy masked array for bathymetry
    bathymetry = np.ma.zeros((ny,nx))
    # read bathymetry and coordinates for variable type
    if (variable == 'z'):
        # get bathymetry at nodes
        bathymetry.data[:,:] = fileID.variables['hz'][:,:].T
        # read latitude and longitude at z-nodes
        lon = fileID.variables['lon_z'][:].copy()
        lat = fileID.variables['lat_z'][:].copy()
    elif variable in ('U','u'):
        # get bathymetry at u-nodes
        bathymetry.data[:,:] = fileID.variables['hu'][:,:].T
        # read latitude and longitude at u-nodes
        lon = fileID.variables['lon_u'][:].copy()
        lat = fileID.variables['lat_u'][:].copy()
    elif variable in ('V','v'):
        # get bathymetry at v-nodes
        bathymetry.data[:,:] = fileID.variables['hv'][:,:].T
        # read latitude and longitude at v-nodes
        lon = fileID.variables['lon_v'][:].copy()
        lat = fileID.variables['lat_v'][:].copy()
    # set bathymetry mask
    bathymetry.mask = (bathymetry.data == 0.0)
    # close the grid file
    fileID.close()
    f.close() if kwargs['compressed'] else None
    return (lon, lat, bathymetry)

# PURPOSE: read elevation file to extract real and imaginary components for
# constituent
def read_elevation_file(input_file, **kwargs):
    """
    Read elevation file to extract real and imaginary components for constituent

    Parameters
    ----------
    input_file: str
        input elevation file
    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    h: float
        tidal elevation
    con: str
        tidal constituent ID
    """
    # set default keyword arguments
    kwargs.setdefault('compressed', False)
    # read the netcdf format tide elevation file
    # reading a combined global solution with localized solutions
    if kwargs['compressed']:
        # read gzipped netCDF4 file
        f = gzip.open(os.path.expanduser(input_file), 'rb')
        fileID = netCDF4.Dataset(uuid.uuid4().hex, 'r', memory=f.read())
    else:
        fileID = netCDF4.Dataset(os.path.expanduser(input_file), 'r')
    # constituent name
    con = fileID.variables['con'][:].tobytes().decode('utf8')
    # variable dimensions
    nx = fileID.dimensions['nx'].size
    ny = fileID.dimensions['ny'].size
    # real and imaginary components of elevation
    h = np.ma.zeros((ny,nx), dtype=np.complex64)
    h.mask = np.zeros((ny,nx), dtype=bool)
    h.data.real[:,:] = fileID.variables['hRe'][:,:].T
    h.data.imag[:,:] = fileID.variables['hIm'][:,:].T
    # close the file
    fileID.close()
    f.close() if kwargs['compressed'] else None
    # return the elevation and constituent
    return (h, con.strip())

# PURPOSE: read transport file to extract real and imaginary components for
# constituent
def read_transport_file(input_file, variable, **kwargs):
    """
    Read transport file to extract real and imaginary components for constituent

    Parameters
    ----------
    input_file: str
        input transport file
    variable: str
        Tidal variable to read

            - ``'u'``: horizontal transport velocities
            - ``'U'``: horizontal depth-averaged transport
            - ``'v'``: vertical transport velocities
            - ``'V'``: vertical depth-averaged transport

    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    tr: float
        tidal transport
    con: str
        tidal constituent ID
    """
    # set default keyword arguments
    kwargs.setdefault('compressed', False)
    # read the netcdf format tide transport file
    # reading a combined global solution with localized solutions
    if kwargs['compressed']:
        # read gzipped netCDF4 file
        f = gzip.open(os.path.expanduser(input_file), 'rb')
        fileID = netCDF4.Dataset(uuid.uuid4().hex, 'r', memory=f.read())
    else:
        fileID = netCDF4.Dataset(os.path.expanduser(input_file), 'r')
    # constituent name
    con = fileID.variables['con'][:].tobytes().decode('utf8')
    # variable dimensions
    nx = fileID.dimensions['nx'].size
    ny = fileID.dimensions['ny'].size
    # real and imaginary components of transport
    tr = np.ma.zeros((ny,nx), dtype=np.complex64)
    tr.mask = np.zeros((ny,nx), dtype=bool)
    if variable in ('U','u'):
        tr.data.real[:,:] = fileID.variables['uRe'][:,:].T
        tr.data.imag[:,:] = fileID.variables['uIm'][:,:].T
    elif variable in ('V','v'):
        tr.data.real[:,:] = fileID.variables['vRe'][:,:].T
        tr.data.imag[:,:] = fileID.variables['vIm'][:,:].T
    # close the file
    fileID.close()
    f.close() if kwargs['compressed'] else None
    # return the transport components and constituent
    return (tr, con.strip())
