#!/usr/bin/env python
u"""
read_tide_model.py
Written by Tyler Sutterley (12/2022)

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
    Updated 12/2022: refactored tide read programs under io
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
import pyTMD.io

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)

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

    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.extract_constants(ilon, ilat,
        grid_file=grid_file, model_file=model_file, EPSG=EPSG,
        **kwargs)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.extend_array(input_array, step_size)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.extend_matrix(input_matrix)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_otis_grid(input_file)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_atlas_grid(input_file)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_netcdf_grid(input_file)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_constituents(input_file, grid=grid)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_otis_elevation(input_file, ic)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_atlas_elevation(input_file, ic, constituent)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_otis_transport(input_file, ic)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_atlas_transport(input_file, ic, constituent)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.create_atlas_mask(xi, yi, mz, local,
        variable=variable)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.interpolate_atlas_model(xi, yi, zi, spacing=spacing)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.combine_atlas_model(xi, yi, zi, pmask, local,
        variable=variable)

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
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.read_netcdf_file(input_file, ic, variable=variable)

# For a rectangular bathymetry grid:
# construct masks for zeta, u and v nodes on a C-grid
def Muv(hz):
    """
    Construct masks for zeta, u and v nodes on a C-grid
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.Muv(hz)

# PURPOSE: Interpolate bathymetry to zeta, u and v nodes on a C-grid
def Huv(hz):
    """
    Interpolate bathymetry to zeta, u and v nodes on a C-grid
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.OTIS.Huv(hz)
