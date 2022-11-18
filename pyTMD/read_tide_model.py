#!/usr/bin/env python
u"""
read_tide_model.py (11/2022)
Reads files for a tidal model and makes initial calculations to run tide program
Includes functions to extract tidal harmonic constants from OTIS tide models for
    given locations

Reads OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/

INPUTS:
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    grid_file: grid file for model
    model_file: model file containing each constituent
    EPSG: projection of tide model data

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
    extrapolate: extrapolate model using nearest-neighbors
    cutoff: extrapolation cutoff in kilometers
        set to np.inf to extrapolate for all points
    grid: binary file type to read
        ATLAS: reading a global solution with localized solutions
        ESR: combined global or local netCDF4 solution
        OTIS: combined global or local solution
    apply_flexure: apply ice flexure scaling factor to constituents

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
    convert_ll_xy.py: converts lat/lon points to and from projected coordinates
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates

UPDATE HISTORY:
    Updated 11/2022: place some imports within try/except statements
        fix variable reads for ATLAS compact data formats
        use f-strings for formatting verbose or ascii output
    Updated 10/2022: invert current tide masks to be True for invalid points
    Updated 06/2022: unit updates in the ESR netCDF4 format
    Updated 05/2022: add functions for using ESR netCDF4 format models
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
        use longcomplex data format to be windows compliant
    Updated 03/2022: invert tide mask to be True for invalid points
        add separate function for resampling ATLAS compact global model
        decode ATLAS compact constituents for Python3 compatibility
        reduce iterative steps when combining ATLAS local models
    Updated 02/2022: use ceiling of masks for interpolation
    Updated 07/2021: added checks that tide model files are accessible
    Updated 06/2021: fix tidal currents for bilinear interpolation
        check for nan points when reading elevation and transport files
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: add extrapolation check where there are no invalid points
        prevent ComplexWarning for fill values when calculating amplitudes
        can read from single constituent TPXO9 ATLAS binary files
        replaced numpy bool/int to prevent deprecation warnings
    Updated 02/2021: set invalid values to nan in extrapolation
        replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added valid data extrapolation with nearest_extrap
    Updated 09/2020: set bounds error to false for regular grid interpolations
        adjust dimensions of input coordinates to be iterable
        use masked arrays with atlas models and grids. make 2' grid with nearest
    Updated 08/2020: check that interpolated points are within range of model
        replaced griddata interpolation with scipy regular grid interpolators
    Updated 07/2020: added function docstrings. separate bilinear interpolation
        update griddata interpolation. changed type variable to keyword argument
    Updated 06/2020: output currents as numpy masked arrays
        use argmin and argmax in bilinear interpolation
    Updated 11/2019: interpolate heights and fluxes to numpy masked arrays
    Updated 09/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 01/2019: decode constituents for Python3 compatibility
    Updated 08/2018: added option grid for using ATLAS outputs that
        combine both global and localized tidal solutions
        added multivariate spline interpolation option
    Updated 07/2018: added different interpolation methods
    Updated 09/2017: Adapted for Python
"""
import os
import copy
import warnings
import numpy as np
import scipy.interpolate
from pyTMD.convert_ll_xy import convert_ll_xy
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
def extract_tidal_constants(ilon, ilat,
    grid_file=None,
    model_file=None,
    EPSG=None,
    **kwargs):
    """
    Reads files for an OTIS-formatted tidal model

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
    model_file: str, list or NoneType, default None
        model file containing each constituent
    EPSG: str or NoneType, default None,
        projection of tide model data
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
    grid: str, default 'OTIS'
        Tide model file type to read

            - ``'ATLAS'``: reading a global solution with localized solutions
            - ``'ESR'``: combined global or local netCDF4 solution
            - ``'OTIS'``: combined global or local solution
    apply_flexure: bool, default False
        Apply ice flexure scaling factor to height constituents

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
    kwargs.setdefault('grid', 'OTIS')
    kwargs.setdefault('apply_flexure', False)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(TYPE='type',METHOD='method',
        EXTRAPOLATE='extrapolate',CUTOFF='cutoff',GRID='grid')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # check that grid file is accessible
    if not os.access(os.path.expanduser(grid_file), os.F_OK):
        raise FileNotFoundError(os.path.expanduser(grid_file))
    # read the OTIS-format tide grid file
    if (kwargs['grid'] == 'ATLAS'):
        # if reading a global solution with localized solutions
        x0,y0,hz0,mz0,iob,dt,pmask,local = read_atlas_grid(grid_file)
        xi,yi,hz = combine_atlas_model(x0,y0,hz0,pmask,local,variable='depth')
        mz = create_atlas_mask(x0,y0,mz0,local,variable='depth')
    elif (kwargs['grid'] == 'ESR'):
        # if reading a single ESR netCDF4 solution
        xi,yi,hz,mz,sf = read_netcdf_grid(grid_file)
    else:
        # if reading a single OTIS solution
        xi,yi,hz,mz,iob,dt = read_tide_grid(grid_file)
    # invert tide mask to be True for invalid points
    mz = np.logical_not(mz).astype(mz.dtype)
    # adjust dimensions of input coordinates to be iterable
    # run wrapper function to convert coordinate systems of input lat/lon
    x,y = convert_ll_xy(np.atleast_1d(ilon), np.atleast_1d(ilat), EPSG, 'F')
    # grid step size of tide model
    dx = xi[1] - xi[0]
    dy = yi[1] - yi[0]

    # create current masks and bathymetry estimates
    if (kwargs['type'] != 'z'):
        mz,mu,mv = Muv(hz)
        hu,hv = Huv(hz)
        # invert current masks to be True for invalid points
        mu = np.logical_not(mu).astype(mu.dtype)
        mv = np.logical_not(mv).astype(mv.dtype)

    # if global: extend limits
    global_grid = False
    # replace original values with extend arrays/matrices
    if ((xi[-1] - xi[0]) == (360.0 - dx)) & (EPSG == '4326'):
        xi = extend_array(xi, dx)
        hz = extend_matrix(hz)
        mz = extend_matrix(mz)
        # set global grid flag
        global_grid = True

    # adjust longitudinal convention of input latitude and longitude
    # to fit tide model convention
    if (np.min(x) < np.min(xi)) & (EPSG == '4326'):
        x[x < 0] += 360.0
    if (np.max(x) > np.max(xi)) & (EPSG == '4326'):
        x[x > 180] -= 360.0
    # determine if any input points are outside of the model bounds
    invalid = (x < xi.min()) | (x > xi.max()) | (y < yi.min()) | (y > yi.max())

    # masks zero values
    hz = np.ma.array(hz,mask=(hz==0))
    if (kwargs['type'] != 'z'):
        # replace original values with extend matrices
        if global_grid:
            hu = extend_matrix(hu)
            hv = extend_matrix(hv)
            mu = extend_matrix(mu)
            mv = extend_matrix(mv)
        # masks zero values
        hu = np.ma.array(hu,mask=(hu==0))
        hv = np.ma.array(hv,mask=(hv==0))

    # interpolate depth and mask to output points
    if (kwargs['method'] == 'bilinear'):
        # use quick bilinear to interpolate values
        D = bilinear_interp(xi, yi, hz, x, y)
        mz1 = bilinear_interp(xi, yi, mz, x, y)
        mz1 = np.ceil(mz1).astype(mz.dtype)
        if (kwargs['type'] != 'z'):
            mu1 = bilinear_interp(xi, yi, mu, x, y)
            mu1 = np.ceil(mu1).astype(mu.dtype)
            mv1 = bilinear_interp(xi, yi, mv, x, y)
            mv1 = np.ceil(mv1).astype(mz.dtype)
    elif (kwargs['method'] == 'spline'):
        # use scipy bivariate splines to interpolate values
        f1=scipy.interpolate.RectBivariateSpline(xi, yi, hz.T, kx=1, ky=1)
        f2=scipy.interpolate.RectBivariateSpline(xi, yi, mz.T, kx=1, ky=1)
        D = f1.ev(x,y)
        mz1 = np.ceil(f2.ev(x,y)).astype(mz.dtype)
        if (kwargs['type'] != 'z'):
            f3=scipy.interpolate.RectBivariateSpline(xi, yi, mu.T, kx=1, ky=1)
            f4=scipy.interpolate.RectBivariateSpline(xi, yi, mv.T, kx=1, ky=1)
            mu1 = np.ceil(f3.ev(x,y)).astype(mu.dtype)
            mv1 = np.ceil(f4.ev(x,y)).astype(mv.dtype)
    else:
        # use scipy regular grid to interpolate values for a given method
        r1 = scipy.interpolate.RegularGridInterpolator((yi,xi), hz,
            method=kwargs['method'], bounds_error=False)
        r2 = scipy.interpolate.RegularGridInterpolator((yi,xi), mz,
            method=kwargs['method'], bounds_error=False, fill_value=0)
        D = r1.__call__(np.c_[y,x])
        mz1 = np.ceil(r2.__call__(np.c_[y,x])).astype(mz.dtype)
        if (kwargs['type'] != 'z'):
            r3 = scipy.interpolate.RegularGridInterpolator((yi,xi), mu,
                method=kwargs['method'], bounds_error=False, fill_value=0)
            r4 = scipy.interpolate.RegularGridInterpolator((yi,xi), mv,
                method=kwargs['method'], bounds_error=False, fill_value=0)
            mu1 = np.ceil(r3.__call__(np.c_[y,x])).astype(mu.dtype)
            mv1 = np.ceil(r4.__call__(np.c_[y,x])).astype(mv.dtype)

    # u and v: velocities in cm/s
    if kwargs['type'] in ('v','u'):
        unit_conv = (D/100.0)
    # U and V: transports in m^2/s
    elif kwargs['type'] in ('V','U'):
        unit_conv = 1.0

    # read and interpolate each constituent
    if isinstance(model_file,list):
        constituents = [read_constituents(m)[0].pop() for m in model_file]
        nc = len(constituents)
    else:
        constituents,nc = read_constituents(model_file, grid=kwargs['grid'])
    # number of output data points
    npts = len(D)
    amplitude = np.ma.zeros((npts,nc))
    amplitude.mask = np.zeros((npts,nc), dtype=bool)
    ph = np.ma.zeros((npts,nc))
    ph.mask = np.zeros((npts,nc), dtype=bool)
    for i,c in enumerate(constituents):
        if (kwargs['type'] == 'z'):
            # read constituent from elevation file
            if (kwargs['grid'] == 'ATLAS'):
                z0,zlocal = read_atlas_elevation(model_file, i, c)
                xi,yi,z = combine_atlas_model(x0, y0, z0, pmask, zlocal,
                    variable='z')
            elif (kwargs['grid'] == 'ESR'):
                z = read_netcdf_file(model_file, i, variable='z')
                # apply flexure scaling
                if kwargs['apply_flexure']:
                    z *= sf
            elif isinstance(model_file,list):
                z = read_elevation_file(model_file[i], 0)
            else:
                z = read_elevation_file(model_file, i)
            # replace original values with extend matrices
            if global_grid:
                z = extend_matrix(z)
            # copy mask to elevation
            z.mask |= mz.astype(bool)
            # interpolate amplitude and phase of the constituent
            z1 = np.ma.zeros((npts), dtype=z.dtype)
            if (kwargs['method'] == 'bilinear'):
                # replace zero values with nan
                z[(z==0) | z.mask] = np.nan
                # use quick bilinear to interpolate values
                z1.data[:] = bilinear_interp(xi, yi, z, x, y,
                    dtype=np.longcomplex)
                # replace nan values with fill_value
                z1.mask = (np.isnan(z1.data) | (mz1.astype(bool)))
                z1.data[z1.mask] = z1.fill_value
            elif (kwargs['method'] == 'spline'):
                # use scipy bivariate splines to interpolate values
                f1 = scipy.interpolate.RectBivariateSpline(xi,yi,
                    z.real.T, kx=1, ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(xi,yi,
                    z.imag.T, kx=1, ky=1)
                z1.data.real = f1.ev(x,y)
                z1.data.imag = f2.ev(x,y)
                # replace zero values with fill_value
                z1.mask = (mz1.astype(bool))
                z1.data[z1.mask] = z1.fill_value
            else:
                # use scipy regular grid to interpolate values
                r1 = scipy.interpolate.RegularGridInterpolator((yi,xi), z,
                    method=kwargs['method'],
                    bounds_error=False,
                    fill_value=z1.fill_value)
                z1 = np.ma.zeros((npts), dtype=z.dtype)
                z1.data[:] = r1.__call__(np.c_[y,x])
                # replace invalid values with fill_value
                z1.mask = (z1.data == z1.fill_value) | (mz1.astype(bool))
                z1.data[z1.mask] = z1.fill_value
            # extrapolate data using nearest-neighbors
            if kwargs['extrapolate'] and np.any(z1.mask):
                # find invalid data points
                inv, = np.nonzero(z1.mask)
                # replace zero values with nan
                z[(z==0) | z.mask] = np.nan
                # extrapolate points within cutoff of valid model points
                z1[inv] = nearest_extrap(xi, yi, z, x[inv], y[inv],
                    dtype=np.longcomplex,
                    cutoff=kwargs['cutoff'],
                    EPSG=EPSG)
            # amplitude and phase of the constituent
            amplitude.data[:,i] = np.abs(z1.data)
            amplitude.mask[:,i] = np.copy(z1.mask)
            ph.data[:,i] = np.arctan2(-np.imag(z1.data), np.real(z1.data))
            ph.mask[:,i] = np.copy(z1.mask)
        elif kwargs['type'] in ('U','u'):
            # read constituent from transport file
            if (kwargs['grid'] == 'ATLAS'):
                u0,v0,uvlocal = read_atlas_transport(model_file, i, c)
                xi,yi,u = combine_atlas_model(x0, y0, u0, pmask, uvlocal,
                    variable='u')
            elif (kwargs['grid'] == 'ESR'):
                u = read_netcdf_file(model_file, i, variable='u')
            elif isinstance(model_file,list):
                u,v = read_transport_file(model_file[i], 0)
            else:
                u,v = read_transport_file(model_file, i)
            # replace original values with extend matrices
            if global_grid:
                u = extend_matrix(u)
            # copy mask to u transports
            u.mask |= mu.astype(bool)
            # x coordinates for u transports
            xu = xi - dx/2.0
            # interpolate amplitude and phase of the constituent
            u1 = np.ma.zeros((npts), dtype=u.dtype)
            if (kwargs['method'] == 'bilinear'):
                # replace zero values with nan
                u[(u==0) | u.mask] = np.nan
                # use quick bilinear to interpolate values
                u1.data[:] = bilinear_interp(xu, yi, u, x, y,
                    dtype=np.longcomplex)
                # replace nan values with fill_value
                u1.mask = (np.isnan(u1.data) | (mu1.astype(bool)))
                u1.data[u1.mask] = u1.fill_value
            elif (kwargs['method'] == 'spline'):
                f1 = scipy.interpolate.RectBivariateSpline(xu, yi,
                    u.real.T, kx=1, ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(xu, yi,
                    u.imag.T, kx=1, ky=1)
                u1.data.real = f1.ev(x,y)
                u1.data.imag = f2.ev(x,y)
                # replace zero values with fill_value
                u1.mask = (mu1.astype(bool))
                u1.data[u1.mask] = u1.fill_value
            else:
                # use scipy regular grid to interpolate values
                r1 = scipy.interpolate.RegularGridInterpolator((yi,xu), u,
                    method=kwargs['method'], bounds_error=False,
                    fill_value=u1.fill_value)
                u1.data[:] = r1.__call__(np.c_[y,x])
                # replace invalid values with fill_value
                u1.mask = (u1.data == u1.fill_value) | (mu1.astype(bool))
                u1.data[u1.mask] = u1.fill_value
            # extrapolate data using nearest-neighbors
            if kwargs['extrapolate'] and np.any(u1.mask):
                # find invalid data points
                inv, = np.nonzero(u1.mask)
                # replace zero values with nan
                u[(u==0) | u.mask] = np.nan
                # extrapolate points within cutoff of valid model points
                u1[inv] = nearest_extrap(xu, yi, u, x[inv], y[inv],
                    dtype=np.longcomplex,
                    cutoff=kwargs['cutoff'],
                    EPSG=EPSG)
            # convert units
            # amplitude and phase of the constituent
            amplitude.data[:,i] = np.abs(u1.data)/unit_conv
            amplitude.mask[:,i] = np.copy(u1.mask)
            ph.data[:,i] = np.arctan2(-np.imag(u1), np.real(u1))
            ph.mask[:,i] = np.copy(u1.mask)
        elif kwargs['type'] in ('V','v'):
            # read constituent from transport file
            if (kwargs['grid'] == 'ATLAS'):
                u0,v0,uvlocal = read_atlas_transport(model_file, i, c)
                xi,yi,v = combine_atlas_model(x0, y0, v0, pmask, uvlocal,
                    variable='v')
            elif (kwargs['grid'] == 'ESR'):
                v = read_netcdf_file(model_file, i, type='v')
            elif isinstance(model_file,list):
                u,v = read_transport_file(model_file[i], 0)
            else:
                u,v = read_transport_file(model_file, i)
            # replace original values with extend matrices
            if global_grid:
                v = extend_matrix(v)
            # copy mask to v transports
            v.mask |= mv.astype(bool)
            # y coordinates for v transports
            yv = yi - dy/2.0
            # interpolate amplitude and phase of the constituent
            v1 = np.ma.zeros((npts), dtype=v.dtype)
            if (kwargs['method'] == 'bilinear'):
                # replace zero values with nan
                v[(v==0) | v.mask] = np.nan
                # use quick bilinear to interpolate values
                v1.data[:] = bilinear_interp(xi, yv, v, x, y,
                    dtype=np.longcomplex)
                # replace nan values with fill_value
                v1.mask = (np.isnan(v1.data) | (mv1.astype(bool)))
                v1.data[v1.mask] = v1.fill_value
            elif (kwargs['method'] == 'spline'):
                f1 = scipy.interpolate.RectBivariateSpline(xi, yv,
                    v.real.T, kx=1, ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(xi, yv,
                    v.imag.T, kx=1, ky=1)
                v1.data.real = f1.ev(x,y)
                v1.data.imag = f2.ev(x,y)
                # replace zero values with fill_value
                v1.mask = (mv1.astype(bool))
                v1.data[v1.mask] = v1.fill_value
            else:
                # use scipy regular grid to interpolate values
                r1 = scipy.interpolate.RegularGridInterpolator((yv,xi), v,
                    method=kwargs['method'],
                    bounds_error=False,
                    fill_value=v1.fill_value)
                v1.data[:] = r1.__call__(np.c_[y,x])
                # replace invalid values with fill_value
                v1.mask = (v1.data == v1.fill_value) | (mv1.astype(bool))
                v1.data[v1.mask] = v1.fill_value
            # extrapolate data using nearest-neighbors
            if kwargs['extrapolate'] and np.any(v1.mask):
                # find invalid data points
                inv, = np.nonzero(v1.mask)
                # replace zero values with nan
                v[(v==0) | v.mask] = np.nan
                # extrapolate points within cutoff of valid model points
                v1[inv] = nearest_extrap(xi, yv, v, x[inv], y[inv],
                    dtype=np.longcomplex,
                    cutoff=kwargs['cutoff'],
                    EPSG=EPSG)
            # convert units
            # amplitude and phase of the constituent
            amplitude.data[:,i] = np.abs(v1.data)/unit_conv
            amplitude.mask[:,i] = np.copy(v1.mask)
            ph.data[:,i] = np.arctan2(-np.imag(v1), np.real(v1))
            ph.mask[:,i] = np.copy(v1.mask)
        # update mask to invalidate points outside model domain
        ph.mask[:,i] |= invalid
        amplitude.mask[:,i] |= invalid

    # convert phase to degrees
    phase = ph*180.0/np.pi
    phase.data[phase.data < 0] += 360.0
    # replace data for invalid mask values
    amplitude.data[amplitude.mask] = amplitude.fill_value
    phase.data[phase.mask] = phase.fill_value
    # return the interpolated values
    return (amplitude,phase,D,constituents)

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

# PURPOSE: read tide grid file
def read_tide_grid(input_file):
    """
    Read grid file to extract model coordinates, bathymetry, masks and indices

    Parameters
    ----------
    input_file: str
        input grid file

    Returns
    -------
    x: float
        x-coordinates of input grid
    y: float
        y-coordinates of input grid
    hz: float
        model bathymetry
    mz: int
        land/water mask
    iob: int
        open boundary index
    dt: float
        time step
    """
    # open the file
    fid = open(os.path.expanduser(input_file),'rb')
    fid.seek(4,0)
    # read data as big endian
    # get model dimensions and limits
    nx, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    ny, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    # extract x and y limits (these could be latitude and longitude)
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    dt, = np.fromfile(fid, dtype=np.dtype('>f4'), count=1)
    # convert longitudinal limits (if x == longitude)
    if (xlim[0] < 0) & (xlim[1] < 0) & (dt > 0):
        xlim += 360.0
    # create x and y arrays arrays (these could be lon and lat values)
    dx = (xlim[1] - xlim[0])/nx
    dy = (ylim[1] - ylim[0])/ny
    x = np.linspace(xlim[0]+dx/2.0, xlim[1]-dx/2.0, nx)
    y = np.linspace(ylim[0]+dy/2.0, ylim[1]-dy/2.0, ny)
    # read nob and iob from file
    nob, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    if (nob == 0):
        fid.seek(20,1)
        iob = []
    else:
        fid.seek(8,1)
        iob=np.fromfile(fid, dtype=np.dtype('>i4'), count=2*nob).reshape(nob,2)
        fid.seek(8,1)
    # read hz matrix
    hz = np.fromfile(fid, dtype=np.dtype('>f4'), count=nx*ny).reshape(ny,nx)
    fid.seek(8,1)
    # read mz matrix
    mz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nx*ny).reshape(ny,nx)
    # close the file
    fid.close()
    # return values
    return (x,y,hz,mz,iob,dt)

# PURPOSE: read tide grid file with localized solutions
def read_atlas_grid(input_file):
    """
    Read ATLAS grid file to extract model coordinates, bathymetry, masks and
    indices for both global and local solutions

    Parameters
    ----------
    input_file: str
        input ATLAS grid file

    Returns
    -------
    x: float
        x-coordinates of input ATLAS grid
    y: float
        y-coordinates of input ATLAS grid
    hz: float
        model bathymetry
    mz: int
        land/water mask
    iob: int
        open boundary index
    dt: float
        time step
    pmask: int
        global mask
    local: dict
        dictionary of local tidal solutions for grid variables

            - ``'depth'``: model bathymetry
    """
    # read the input file to get file information
    fd = os.open(os.path.expanduser(input_file), os.O_RDONLY)
    file_info = os.fstat(fd)
    fid = os.fdopen(fd, 'rb')
    fid.seek(4,0)
    # read data as big endian
    # get model dimensions and limits
    nx, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    ny, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    # extract latitude and longitude limits
    lats = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    lons = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    dt, = np.fromfile(fid, dtype=np.dtype('>f4'), count=1)
    # create lon and lat arrays
    dlon = (lons[1] - lons[0])/nx
    dlat = (lats[1] - lats[0])/ny
    x = np.linspace(lons[0]+dlon/2.0,lons[1]-dlon/2.0,nx)
    y = np.linspace(lats[0]+dlat/2.0,lats[1]-dlat/2.0,ny)
    # read nob and iob from file
    nob, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    if (nob == 0):
        fid.seek(20,1)
        iob = []
    else:
        fid.seek(8,1)
        iob=np.fromfile(fid, dtype=np.dtype('>i4'), count=2*nob).reshape(nob,2)
        fid.seek(8,1)
    # read hz matrix
    hz = np.fromfile(fid, dtype=np.dtype('>f4'), count=nx*ny).reshape(ny,nx)
    fid.seek(8,1)
    # read mz matrix
    mz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nx*ny).reshape(ny,nx)
    fid.seek(8,1)
    # read pmask matrix
    pmask = np.fromfile(fid, dtype=np.dtype('>i4'), count=nx*ny).reshape(ny,nx)
    fid.seek(4,1)
    # read local models
    nmod = 0
    local = {}
    # while the file position is not at the end of file
    while (fid.tell() < file_info.st_size):
        # add 1 to number of models
        fid.seek(4,1)
        nmod += 1
        # get local model dimensions and limits
        nx1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        ny1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nd, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        # extract latitude and longitude limits of local model
        lt1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        ln1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        # extract name
        name = fid.read(20).strip()
        fid.seek(8,1)
        iz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nd)
        jz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nd)
        fid.seek(8,1)
        depth = np.ma.zeros((ny1,nx1))
        depth.mask = np.ones((ny1,nx1), dtype=bool)
        depth.data[jz-1,iz-1] = np.fromfile(fid, dtype=np.dtype('>f4'), count=nd)
        depth.mask[jz-1,iz-1] = False
        fid.seek(4,1)
        # save to dictionary
        local[name] = dict(lon=ln1, lat=lt1, depth=depth)
    # close the file
    fid.close()
    # return values
    return (x,y,hz,mz,iob,dt,pmask,local)

# PURPOSE: read grid file
def read_netcdf_grid(input_file):
    """
    Read netCDF4 grid file to extract model coordinates, bathymetry,
    masks and flexure scaling factors

    Parameters
    ----------
    input_file: str
        input grid file

    Returns
    -------
    x: float
        x-coordinates of input grid
    y: float
        y-coordinates of input grid
    hz: float
        model bathymetry
    mz: int
        land/water mask
    sf: float
        scaling factor for applying ice flexure
    """
    # read the netcdf format tide grid file
    fileID = netCDF4.Dataset(os.path.expanduser(input_file),'r')
    # read coordinates and flip y orientation
    x = fileID.variables['x'][:].copy()
    y = fileID.variables['y'][::-1].copy()
    # read water column thickness and flip y orientation
    hz = fileID.variables['wct'][::-1,:].copy()
    # read mask and flip y orientation
    mz = fileID.variables['mask'][::-1,:].copy()
    # read flexure and convert from percent to scale factor
    sf = fileID.variables['flexure'][::-1,:]/100.0
    # update bathymetry and scale factor masks
    hz.mask = (hz.data == 0.0)
    sf.mask = (sf.data == 0.0)
    # close the grid file
    fileID.close()
    # return values
    return (x,y,hz,mz,sf)

# PURPOSE: read list of constituents from an elevation or transport file
def read_constituents(input_file, grid='OTIS'):
    """
    Read the list of constituents from an elevation or transport file

    Parameters
    ----------
    input_file: str
        input tidal file
    grid: str, default 'OTIS'
        Tide model file type to read

            - ``'ATLAS'``: reading a global solution with localized solutions
            - ``'ESR'``: combined global or local netCDF4 solution
            - ``'OTIS'``: combined global or local solution

    Returns
    -------
    constituents: list
        list of tidal constituent IDs
    nc: int
        number of constituents
    """
    # check that model file is accessible
    if not os.access(os.path.expanduser(input_file), os.F_OK):
        raise FileNotFoundError(os.path.expanduser(input_file))
    if (grid == 'ESR'):
        # open the netCDF4 file
        fid = netCDF4.Dataset(os.path.expanduser(input_file),'r')
        constituents = fid.variables['constituents'].constituent_order.split()
        nc = len(constituents)
        fid.close()
    else:
        # open the file
        fid = open(os.path.expanduser(input_file),'rb')
        ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
        fid.seek(16,1)
        constituents = [c.decode("utf8").rstrip() for c in fid.read(nc*4).split()]
        fid.close()
    return (constituents,nc)

# PURPOSE: read elevation file to extract real and imaginary components for
# constituent
def read_elevation_file(input_file,ic):
    """
    Read elevation file to extract real and imaginary components for constituent

    Parameters
    ----------
    input_file: str
        input elevation file
    ic: int
        index of consituent

    Returns
    -------
    h: float
        tidal elevation
    """
    # open the file
    fid = open(os.path.expanduser(input_file),'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    # extract x and y limits
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    # skip records to constituent
    nskip = ic*(nx*ny*8+8) + 8 + ll - 28
    fid.seek(nskip,1)
    # real and imaginary components of elevation
    h = np.ma.zeros((ny,nx), dtype=np.complex64)
    h.mask = np.zeros((ny,nx), dtype=bool)
    for i in range(ny):
        temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nx)
        h.data.real[i,:] = temp[0:2*nx-1:2]
        h.data.imag[i,:] = temp[1:2*nx:2]
    # update mask for nan values
    h.mask[np.isnan(h.data)] = True
    # replace masked values with fill value
    h.data[h.mask] = h.fill_value
    # close the file
    fid.close()
    # return the elevation
    return h

# PURPOSE: read elevation file with localized solutions to extract real and
# imaginary components for constituent
def read_atlas_elevation(input_file, ic, constituent):
    """
    Read elevation file with localized solutions to extract real and imaginary
    components for constituent

    Parameters
    ----------
    input_file: str
        input ATLAS elevation file
    ic: int
        index of consituent
    constituent: str
        tidal constituent ID

    Returns
    -------
    h: float
        global tidal elevation
    local: dict
        dictionary of local tidal solutions for elevation variables

            - ``'z'``: tidal elevation
    """
    # read the input file to get file information
    fd = os.open(os.path.expanduser(input_file), os.O_RDONLY)
    file_info = os.fstat(fd)
    fid = os.fdopen(fd, 'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    # extract x and y limits
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    # skip records to constituent
    nskip = 8 + nc*4 + ic*(nx*ny*8 + 8)
    fid.seek(nskip,1)
    # real and imaginary components of elevation
    h = np.ma.zeros((ny,nx), dtype=np.complex64)
    h.mask = np.zeros((ny,nx), dtype=bool)
    for i in range(ny):
        temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nx)
        h.data.real[i,:] = temp[0:2*nx-1:2]
        h.data.imag[i,:] = temp[1:2*nx:2]
    # skip records after constituent
    nskip = (nc-ic-1)*(nx*ny*8 + 8) + 4
    fid.seek(nskip,1)
    # read local models to find constituent
    nmod = 0
    local = {}
    # while the file position is not at the end of file
    while (fid.tell() < file_info.st_size):
        # add 1 to number of models
        fid.seek(4,1)
        nmod += 1
        # get local model dimensions and limits
        nx1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        ny1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nc1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nz, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        # extract latitude and longitude limits of local model
        lt1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        ln1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        # extract constituents for localized solution
        cons = fid.read(nc1*4).strip().decode("utf8").split()
        # check if constituent is in list of localized solutions
        if (constituent in cons):
            ic1, = [i for i,c in enumerate(cons) if (c == constituent)]
            # extract name
            name = fid.read(20).strip()
            fid.seek(8,1)
            iz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nz)
            jz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nz)
            # skip records to constituent
            nskip = 8 + ic1*(8*nz + 8)
            fid.seek(nskip,1)
            # real and imaginary components of elevation
            h1 = np.ma.zeros((ny1,nx1), fill_value=np.nan, dtype=np.complex64)
            h1.mask = np.ones((ny1,nx1), dtype=bool)
            temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nz)
            h1.data.real[jz-1,iz-1] = temp[0:2*nz-1:2]
            h1.data.imag[jz-1,iz-1] = temp[1:2*nz:2]
            h1.mask[jz-1,iz-1] = False
            # save constituent to dictionary
            local[name] = dict(lon=ln1,lat=lt1,z=h1)
            # skip records after constituent
            nskip = (nc1-ic1-1)*(8*nz + 8) + 4
            fid.seek(nskip,1)
        else:
            # skip records for local model if constituent not in list
            nskip = 40 + 16*nz + (nc1-1)*(8*nz + 8)
            fid.seek(nskip,1)
    # close the file
    fid.close()
    # return the elevation
    return (h,local)

# PURPOSE: read transport file to extract real and imaginary components for
# constituent
def read_transport_file(input_file,ic):
    """
    Read transport file to extract real and imaginary components for constituent

    Parameters
    ----------
    input_file: str
        input transport file
    ic: int
        index of consituent

    Returns
    -------
    u: float
        zonal tidal transport
    v: float
        meridional zonal transport
    """
    # open the file
    fid = open(os.path.expanduser(input_file),'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    # extract x and y limits
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    # skip records to constituent
    nskip = ic*(nx*ny*16+8) + 8 + ll - 28
    fid.seek(nskip,1)
    # real and imaginary components of transport
    u = np.ma.zeros((ny,nx), dtype=np.complex64)
    u.mask = np.zeros((ny,nx), dtype=bool)
    v = np.ma.zeros((ny,nx), dtype=np.complex64)
    v.mask = np.zeros((ny,nx), dtype=bool)
    for i in range(ny):
        temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=4*nx)
        u.data.real[i,:] = temp[0:4*nx-3:4]
        u.data.imag[i,:] = temp[1:4*nx-2:4]
        v.data.real[i,:] = temp[2:4*nx-1:4]
        v.data.imag[i,:] = temp[3:4*nx:4]
    # update mask for nan values
    u.mask[np.isnan(u.data)] = True
    v.mask[np.isnan(v.data)] = True
    # replace masked values with fill value
    u.data[u.mask] = u.fill_value
    v.data[v.mask] = v.fill_value
    # close the file
    fid.close()
    # return the transport components
    return (u,v)

# PURPOSE: read transport file with localized solutions to extract real and
# imaginary components for constituent
def read_atlas_transport(input_file, ic, constituent):
    """
    Read transport file with localized solutions to extract real and imaginary
    components for constituent

    Parameters
    ----------
    input_file: str
        input ATLAS transport file
    ic: int
        index of consituent
    constituent: str
        tidal constituent ID

    Returns
    -------
    u: float
        global zonal tidal transport
    v: float
        global meridional zonal transport
    local: dict
        dictionary of local tidal solutions for transport variables

            - ``'u'``: zonal tidal transport
            - ``'v'``: meridional zonal transport
    """
    # read the input file to get file information
    fd = os.open(os.path.expanduser(input_file), os.O_RDONLY)
    file_info = os.fstat(fd)
    fid = os.fdopen(fd, 'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    # extract x and y limits
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    # skip records to constituent
    nskip = 8 + nc*4 + ic*(nx*ny*16 + 8)
    fid.seek(nskip,1)
    # real and imaginary components of transport
    u = np.ma.zeros((ny,nx), dtype=np.complex64)
    u.mask = np.zeros((ny,nx), dtype=bool)
    v = np.ma.zeros((ny,nx), dtype=np.complex64)
    v.mask = np.zeros((ny,nx), dtype=bool)
    for i in range(ny):
        temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=4*nx)
        u.data.real[i,:] = temp[0:4*nx-3:4]
        u.data.imag[i,:] = temp[1:4*nx-2:4]
        v.data.real[i,:] = temp[2:4*nx-1:4]
        v.data.imag[i,:] = temp[3:4*nx:4]
    # skip records after constituent
    nskip = (nc-ic-1)*(nx*ny*16 + 8) + 4
    fid.seek(nskip,1)
    # read local models to find constituent
    nmod = 0
    local = {}
    # while the file position is not at the end of file
    while (fid.tell() < file_info.st_size):
        # add 1 to number of models
        fid.seek(4,1)
        nmod += 1
        # get local model dimensions and limits
        nx1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        ny1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nc1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nu, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nv, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        # extract latitude and longitude limits of local model
        lt1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        ln1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        # extract constituents for localized solution
        cons = fid.read(nc1*4).strip().decode("utf8").split()
        # check if constituent is in list of localized solutions
        if (constituent in cons):
            ic1, = [i for i,c in enumerate(cons) if (c == constituent)]
            # extract name
            name = fid.read(20).strip()
            fid.seek(8,1)
            iu = np.fromfile(fid, dtype=np.dtype('>i4'), count=nu)
            ju = np.fromfile(fid, dtype=np.dtype('>i4'), count=nu)
            fid.seek(8,1)
            iv = np.fromfile(fid, dtype=np.dtype('>i4'), count=nv)
            jv = np.fromfile(fid, dtype=np.dtype('>i4'), count=nv)
            # skip records to constituent
            nskip = 8 + ic1*(8*nu + 8*nv + 16)
            fid.seek(nskip,1)
            # real and imaginary components of u transport
            u1 = np.ma.zeros((ny1,nx1), fill_value=np.nan, dtype=np.complex64)
            u1.mask = np.ones((ny1,nx1), dtype=bool)
            tmpu = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nu)
            u1.data.real[ju-1,iu-1] = tmpu[0:2*nu-1:2]
            u1.data.imag[ju-1,iu-1] = tmpu[1:2*nu:2]
            u1.mask[ju-1,iu-1] = False
            fid.seek(8,1)
            # real and imaginary components of v transport
            v1 = np.ma.zeros((ny1,nx1), fill_value=np.nan, dtype=np.complex64)
            v1.mask = np.ones((ny1,nx1), dtype=bool)
            tmpv = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nv)
            v1.data.real[jv-1,iv-1] = tmpv[0:2*nv-1:2]
            v1.data.imag[jv-1,iv-1] = tmpv[1:2*nv:2]
            v1.mask[jv-1,iv-1] = False
            # save constituent to dictionary
            local[name] = dict(lon=ln1,lat=lt1,u=u1,v=v1)
            # skip records after constituent
            nskip = (nc1-ic1-1)*(8*nu + 8*nv + 16) + 4
            fid.seek(nskip,1)
        else:
            # skip records for local model if constituent not in list
            nskip = 56 + 16*nu + 16*nv + (nc1-1)*(8*nu + 8*nv + 16)
            fid.seek(nskip,1)
    # close the file
    fid.close()
    # return the transport components
    return (u,v,local)

# PURPOSE: create a 2 arc-minute grid mask from mz and depth variables
def create_atlas_mask(xi, yi, mz, local, variable=None):
    """
    Creates a high-resolution grid mask from model variables

    Parameters
    ----------
    xi: float
        input x-coordinates of global tide model
    yi: float
        input y-coordinates of global tide model
    mz: int
        global land/water mask
    local: dict
        dictionary of local tidal solutions
    variable: str or NoneType, default None
        key for variable within each local solution

            - ``'depth'``: model bathymetry
            - ``'z'``: tidal elevation
            - ``'u'``: zonal tidal transport
            - ``'v'``: meridional zonal transport

    Returns
    -------
    x30: float
        x-coordinates of high-resolution tide model
    y30: float
        y-coordinates of high-resolution tide model
    m30: int
        high-resolution land/water mask
    """
    # create 2 arc-minute grid dimensions
    d30 = 1.0/30.0
    x30 = np.arange(d30/2.0, 360.0+d30/2.0, d30)
    y30 = np.arange(-90.0+d30/2.0, 90.0+d30/2.0, d30)
    # interpolate global mask to create initial 2 arc-minute mask
    xcoords=np.clip((len(xi)-1)*(x30-xi[0])/(xi[-1]-xi[0]),0,len(xi)-1)
    ycoords=np.clip((len(yi)-1)*(y30-yi[0])/(yi[-1]-yi[0]),0,len(yi)-1)
    IY,IX = np.meshgrid(np.around(ycoords), np.around(xcoords), indexing='ij')
    # interpolate with nearest-neighbors
    m30 = np.ma.zeros((len(y30),len(x30)), dtype=np.int8,fill_value=0)
    m30.data[:,:] = mz[IY.astype(np.int32), IX.astype(np.int32)]
    # iterate over localized solutions to fill in high-resolution coastlines
    for key,val in local.items():
        # shape of local variable
        ny,nx = np.shape(val[variable])
        # correct limits for local grid
        lon0 = np.floor(val['lon'][0]/d30)*d30
        lat0 = np.floor(val['lat'][0]/d30)*d30
        # create latitude and longitude for local model
        xi = lon0 + np.arange(nx)*d30
        yi = lat0 + np.arange(ny)*d30
        IX,IY = np.meshgrid(xi, yi)
        # local model output
        validy,validx = np.nonzero(np.logical_not(val[variable].mask))
        # check if any model longitudes are -180:180
        X = np.where(IX[validy,validx] <= 0.0,
            IX[validy,validx] + 360.0, IX[validy,validx])
        # grid indices of local model
        ii = ((X - x30[0])//d30).astype('i')
        jj = ((IY[validy,validx] - y30[0])//d30).astype('i')
        # fill global mask with regional solution
        m30[jj,ii] = 1
    # return the 2 arc-minute mask
    m30.mask = (m30.data == m30.fill_value)
    return m30

# PURPOSE: resample global solution to higher-resolution
def interpolate_atlas_model(xi, yi, zi, spacing=1.0/30.0):
    """
    Interpolates global ATLAS tidal solutions into a
    higher-resolution sampling

    Parameters
    ----------
    xi: float
        input x-coordinates of global tide model
    yi: float
        input y-coordinates of global tide model
    zi: float
        global tide model data
    spacing: float
        output grid spacing

    Returns
    -------
    xs: float
        x-coordinates of high-resolution tide model
    ys: float
        y-coordinates of high-resolution tide model
    zs: float
        high-resolution tidal solution for variable
    """
    # create resampled grid dimensions
    xs = np.arange(spacing/2.0, 360.0+spacing/2.0, spacing)
    ys = np.arange(-90.0+spacing/2.0, 90.0+spacing/2.0, spacing)
    # interpolate global solution
    zs = np.ma.zeros((len(ys),len(xs)), dtype=zi.dtype)
    zs.mask = np.zeros((len(ys),len(xs)), dtype=bool)
    # test if combining elevation/transport variables with complex components
    if np.iscomplexobj(zs):
        f1 = scipy.interpolate.RectBivariateSpline(xi, yi, zi.real.T, kx=1,ky=1)
        f2 = scipy.interpolate.RectBivariateSpline(xi, yi, zi.imag.T, kx=1,ky=1)
        zs.data.real[:,:] = f1(xs,ys).T
        zs.data.imag[:,:] = f2(xs,ys).T
    else:
        f = scipy.interpolate.RectBivariateSpline(xi, yi, zi.T, kx=1,ky=1)
        zs.data[:,:] = f(xs,ys).T
    # return resampled solution and coordinates
    return (xs,ys,zs)

# PURPOSE: combines global and local atlas solutions
def combine_atlas_model(xi, yi, zi, pmask, local, variable=None):
    """
    Combines global and local ATLAS tidal solutions into a single
    high-resolution solution

    Parameters
    ----------
    xi: float
        input x-coordinates of global tide model
    yi: float
        input y-coordinates of global tide model
    zi: float
        global tide model data
    pmask: int
        global mask
    local: dict
        dictionary of local tidal solutions
    variable: str or NoneType, default None
        key for variable within each local solution

            - ``'depth'``: model bathymetry
            - ``'z'``: tidal elevation
            - ``'u'``: zonal tidal transport
            - ``'v'``: meridional zonal transport

    Returns
    -------
    x30: float
        x-coordinates of high-resolution tide model
    y30: float
        y-coordinates of high-resolution tide model
    z30: float
        combined high-resolution tidal solution for variable
    """
    # create 2 arc-minute grid dimensions
    d30 = 1.0/30.0
    # interpolate global solution to 2 arc-minute solution
    x30,y30,z30 = interpolate_atlas_model(xi, yi, zi, spacing=d30)
    # iterate over localized solutions
    for key,val in local.items():
        # shape of local variable
        ny,nx = np.shape(val[variable])
        # correct limits for local grid
        lon0 = np.floor(val['lon'][0]/d30)*d30
        lat0 = np.floor(val['lat'][0]/d30)*d30
        # create latitude and longitude for local model
        xi = lon0 + np.arange(nx)*d30
        yi = lat0 + np.arange(ny)*d30
        IX,IY = np.meshgrid(xi,yi)
        # local model output
        validy,validx = np.nonzero(np.logical_not(val[variable].mask))
        # check if any model longitudes are -180:180
        X = np.where(IX[validy,validx] <= 0.0,
            IX[validy,validx] + 360.0, IX[validy,validx])
        # grid indices of local model
        ii = ((X - x30[0])//d30).astype('i')
        jj = ((IY[validy,validx] - y30[0])//d30).astype('i')
        # fill global mask with regional solution
        z30.data[jj,ii] = val[variable][validy,validx]
    # return 2 arc-minute solution and coordinates
    return (x30,y30,z30)

# PURPOSE: read netCDF4 file to extract real and imaginary components for
# constituent
def read_netcdf_file(input_file, ic, variable=None):
    """
    Read netCDF4 file to extract real and imaginary components for constituent

    Parameters
    ----------
    input_file: str
        input transport file
    ic: int
        index of consituent
    variable: str or NoneType, default None
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'U'``: horizontal depth-averaged transport
            - ``'v'``: vertical transport velocities
            - ``'V'``: vertical depth-averaged transport

    Returns
    -------
    hc: complex
        complex form of tidal constituent oscillation
    """
    # read the netcdf format tide grid file
    fileID = netCDF4.Dataset(os.path.expanduser(input_file), 'r')
    # variable dimensions
    nx = fileID.dimensions['x'].size
    ny = fileID.dimensions['y'].size
    # real and imaginary components of tidal constituent
    hc = np.ma.zeros((ny,nx), dtype=np.complex64)
    hc.mask = np.zeros((ny,nx), dtype=bool)
    # extract constituent and flip y orientation
    if (variable == 'z'):
        hc.data.real[:,:] = fileID.variables['hRe'][ic,::-1,:]
        hc.data.imag[:,:] = -fileID.variables['hIm'][ic,::-1,:]
    elif variable in ('U','u'):
        hc.data.real[:,:] = fileID.variables['uRe'][ic,::-1,:]
        hc.data.imag[:,:] = -fileID.variables['uIm'][ic,::-1,:]
    elif variable in ('V','v'):
        hc.data.real[:,:] = fileID.variables['vRe'][ic,::-1,:]
        hc.data.imag[:,:] = -fileID.variables['vIm'][ic,::-1,:]
    # close the file
    fileID.close()
    # return output variables
    return hc

# For a rectangular bathymetry grid:
# construct masks for zeta, u and v nodes on a C-grid
def Muv(hz):
    """
    Construct masks for zeta, u and v nodes on a C-grid
    """
    ny,nx = np.shape(hz)
    mz = (hz > 0).astype(int)
    # x-indices
    indx = np.zeros((nx), dtype=int)
    indx[:-1] = np.arange(1,nx)
    indx[-1] = 0
    # y-indices
    indy = np.zeros((ny), dtype=int)
    indy[:-1] = np.arange(1,ny)
    indy[-1] = 0
    # calculate mu and mv
    mu = np.zeros((ny,nx), dtype=int)
    mv = np.zeros((ny,nx), dtype=int)
    mu[indy,:] = mz*mz[indy,:]
    mv[:,indx] = mz*mz[:,indx]
    return (mu,mv,mz)

# PURPOSE: Interpolate bathymetry to zeta, u and v nodes on a C-grid
def Huv(hz):
    """
    Interpolate bathymetry to zeta, u and v nodes on a C-grid
    """
    ny,nx = np.shape(hz)
    mu,mv,mz = Muv(hz)
    # x-indices
    indx = np.zeros((nx), dtype=int)
    indx[0] = nx-1
    indx[1:] = np.arange(1,nx)
    # y-indices
    indy = np.zeros((ny), dtype=int)
    indy[0] = ny-1
    indy[1:] = np.arange(1,ny)
    # calculate hu and hv
    hu = mu*(hz + hz[indy,:])/2.0
    hv = mv*(hz + hz[:,indx])/2.0
    return (hu,hv)
