#!/usr/bin/env python
u"""
read_FES_model.py (07/2020)
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
    directory: data directory for tide data files
    model_files: list of model files for each constituent (can be gzipped)
    TYPE: tidal variable to run
        z: heights
        u: horizontal transport velocities
        v: vertical transport velocities
    VERSION: model version to run
        FES1999
        FES2004
        FES2012
        FES2014

OPTIONS:
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, cubic, nearest: scipy griddata interpolations
    GZIP: input ascii or netCDF4 files are compressed
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
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    bilinear_interp.py: bilinear interpolation of data to specified coordinates

UPDATE HISTORY:
    Written 07/2020
"""
import os
import gzip
import netCDF4
import numpy as np
import scipy.interpolate
from pyTMD.bilinear_interp import bilinear_interp

#-- PURPOSE: extract tidal harmonic constants from tide models at coordinates
def extract_FES_constants(ilon, ilat, directory, model_files,
    TYPE='z', VERSION=None, METHOD='spline', GZIP=True, SCALE=1):
    """
    Reads files for an ascii or netCDF4 tidal model
    Makes initial calculations to run the tide program
    Spatially interpolates tidal constituents to input coordinates

    Arguments
    ---------
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    directory: data directory for tide data files
    grid_file: grid file for model (can be gzipped)
    model_files: list of model files for each constituent (can be gzipped)

    Keyword arguments
    -----------------
    TYPE: tidal variable to run
        z: heights
        u: horizontal transport velocities
        v: vertical transport velocities
    VERSION: model version to run
        FES1999
        FES2004
        FES2012
        FES2014
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, cubic, nearest: scipy griddata interpolations
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
    #-- number of constituents
    nc = len(model_files)
    #-- amplitude and phase
    amplitude = np.ma.zeros((npts,nc))
    amplitude.mask = np.zeros((npts,nc),dtype=np.bool)
    phase = np.ma.zeros((npts,nc))
    phase.mask = np.zeros((npts,nc),dtype=np.bool)
    #-- read and interpolate each constituent
    for i,fi in enumerate(model_files):
        #-- read constituent from elevation file
        if VERSION in ('FES1999','FES2004'):
            hc,lon,lat = read_ascii_file(os.path.join(directory,fi),
                GZIP=GZIP,TYPE=TYPE,VERSION=VERSION)
        elif VERSION in ('FES2012','FES2014'):
            hc,lon,lat = read_netcdf_file(os.path.join(directory,fi),
                GZIP=GZIP,TYPE=TYPE,VERSION=VERSION)
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
            hci.mask[:] = np.isnan(hci.data)
            hci.data[hci.mask] = hci.fill_value
        elif (METHOD == 'spline'):
            #-- interpolate complex form of the constituent with scipy
            f1=scipy.interpolate.RectBivariateSpline(lon,lat,hc.data.real.T,kx=1,ky=1)
            f2=scipy.interpolate.RectBivariateSpline(lon,lat,hc.data.imag.T,kx=1,ky=1)
            f3=scipy.interpolate.RectBivariateSpline(lon,lat,hc.mask.T,kx=1,ky=1)
            hci.data.real[:] = f1.ev(ilon,ilat)
            hci.data.imag[:] = f2.ev(ilon,ilat)
            hci.mask[:] = f3.ev(ilon,ilat).astype(np.bool)
        else:
            #-- create mesh grids of latitude and longitude
            X,Y = np.meshgrid(lon,lat)
            interp_points = np.c_[X.flatten(),Y.flatten()]
            #-- replace invalid values with nan
            hc[hc.mask] = np.nan
            #-- interpolate complex form of the constituent with scipy
            hci.data[:] = scipy.interpolate.griddata(interp_points,
                hc.flatten(), np.c_[ilon,ilat], method=METHOD)
            #-- mask invalid values
            hci.mask[:] = np.isnan(hci.data[:])
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

#-- PURPOSE: read FES ascii tide model grid files
def read_ascii_file(input_file,GZIP=False,TYPE=None,VERSION=None):
    """
    Read FES (Finite Element Solution) tide model file

    Arguments
    ---------
    input_file: model file

    Keyword arguments
    -----------------
    GZIP: input files are compressed
    VERSION: model version

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
    #-- longitude range (lonmin, lonmax)
    lonmin,lonmax = np.array(file_contents[0].split(), dtype=np.float)
    #-- latitude range (latmin, latmax)
    latmin,latmax = np.array(file_contents[1].split(), dtype=np.float)
    #-- grid step size (dlon, dlat)
    dlon,dlat = np.array(file_contents[2].split(), dtype=np.float)
    #-- grid dimensions (nlon, nlat)
    nlon,nlat = np.array(file_contents[3].split(), dtype=np.int)
    #-- mask fill value
    masked_values = file_contents[4].split()
    fill_value = np.float(masked_values[0])
    #-- create output variables
    lat = np.linspace(latmin, latmax, nlat)
    lon = np.linspace(lonmin,lonmax,nlon)
    amp = np.ma.zeros((nlat,nlon),fill_value=fill_value,dtype=np.float32)
    ph = np.ma.zeros((nlat,nlon),fill_value=fill_value,dtype=np.float32)
    #-- create masks for output variables (0=valid)
    amp.mask = np.zeros((nlat,nlon),dtype=np.bool)
    ph.mask = np.zeros((nlat,nlon),dtype=np.bool)
    #-- starting line to fill amplitude and phase variables
    i1 = 5
    #-- for each latitude
    for i in range(nlat):
        for j in range(nlon//30):
            j1 = j*30
            amp.data[i,j1:j1+30]=np.array(file_contents[i1].split(),dtype='f')
            ph.data[i,j1:j1+30]=np.array(file_contents[i1+1].split(),dtype='f')
            i1 += 2
        #-- add last tidal variables
        j1 = (j+1)*30
        j2 = nlon % 30
        amp.data[i,j1:j1+j2] = np.array(file_contents[i1].split(),dtype='f')
        ph.data[i,j1:j1+j2] = np.array(file_contents[i1+1].split(),dtype='f')
        i1 += 2
    #-- calculate complex form of constituent oscillation
    hc = amp*np.exp(-1j*ph*np.pi/180.0)
    #-- set masks
    hc.mask = (amp.data == amp.fill_value) | (ph.data == ph.fill_value)
    #-- return output variables
    return (hc,lon,lat)

#-- PURPOSE: read FES netCDF4 tide model files
def read_netcdf_file(input_file,GZIP=False,TYPE=None,VERSION=None):
    """
    Read FES (Finite Element Solution) tide model netCDF4 file

    Arguments
    ---------
    input_file: model file

    Keyword arguments
    -----------------
    GZIP: input files are compressed
    VERSION: model version
    TYPE: tidal variable to run
        z: heights
        u: horizontal transport velocities
        v: vertical transport velocities

    Returns
    -------
    hc: complex form of tidal constituent oscillation
    lon: longitude of tidal model
    lat: latitude of tidal model
    """
    #-- read the netcdf format tide elevation file
    #-- reading a combined global solution with localized solutions
    if GZIP:
        f = gzip.open(input_file,'rb')
        fileID = netCDF4.Dataset(input_file,'r',memory=f.read())
    else:
        fileID = netCDF4.Dataset(input_file,'r')
    #-- variable dimensions for each model
    if (VERSION == 'FES2012'):
        lon = fileID.variables['longitude'][:]
        lat = fileID.variables['latitude'][:]
    elif (VERSION == 'FES2014'):
        lon = fileID.variables['lon'][:]
        lat = fileID.variables['lat'][:]
    #-- amplitude and phase components for each type
    if (TYPE == 'z'):
        amp = fileID.variables['amplitude'][:]
        ph = fileID.variables['phase'][:]
    elif (TYPE == 'u'):
        amp = fileID.variables['Ua'][:]
        ph = fileID.variables['Ug'][:]
    elif (TYPE == 'v'):
        amp = fileID.variables['Va'][:]
        ph = fileID.variables['Vg'][:]
    #-- close the file
    fileID.close()
    f.close() if GZIP else None
    #-- calculate complex form of constituent oscillation
    hc = amp*np.exp(-1j*ph*np.pi/180.0)
    #-- set masks
    hc.mask = (amp.data == amp.fill_value) | (ph.data == ph.fill_value)
    #-- return output variables
    return (hc,lon,lat)
