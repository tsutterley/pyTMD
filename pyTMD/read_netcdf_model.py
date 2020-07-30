#!/usr/bin/env python
u"""
read_netcdf_model.py (07/2020)
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
    directory: data directory for tide data files
    grid_file: grid file for model (can be gzipped)
    model_files: list of model files for each constituent (can be gzipped)

OPTIONS:
    TYPE: tidal variable to run
        z: heights
        u: horizontal transport velocities
        U: horizontal depth-averaged transport
        v: vertical transport velocities
        V: vertical depth-averaged transport
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, cubic, nearest: scipy griddata interpolations
    GZIP: input netCDF4 files are compressed
    SCALE: scaling factor for converting to output units

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
    bilinear_interp.py: bilinear interpolation of data to specified coordinates

UPDATE HISTORY:
    Updated 07/2020: added function docstrings. separate bilinear interpolation
        changed TYPE variable to keyword argument. update griddata interpolation
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Written 09/2019
"""
import os
import gzip
import netCDF4
import numpy as np
import scipy.interpolate
from pyTMD.bilinear_interp import bilinear_interp

#-- PURPOSE: extract tidal harmonic constants from tide models at coordinates
def extract_netcdf_constants(ilon, ilat, directory, grid_file, model_files,
    TYPE='z', METHOD='spline', GZIP=True, SCALE=1):
    """
    Reads files for a netCDF4 tidal model
    Makes initial calculations to run the tide program
    Spatially interpolates tidal constituents to input coordinates

    Arguments
    ---------
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    directory: data directory for tide data files
    grid_file: grid file for model (can be gzipped)
    model_files: list of model files for each constituent (can be gzipped)
    TYPE: tidal variable to run
        z: heights
        u: horizontal transport velocities
        U: horizontal depth-averaged transport
        v: vertical transport velocities
        V: vertical depth-averaged transport

    Keyword arguments
    -----------------
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, cubic, nearest: scipy griddata interpolations
    GZIP: input netCDF4 files are compressed
    SCALE: scaling factor for converting to output units

    Returns
    -------
    amplitude: amplitudes of tidal constituents
    phase: phases of tidal constituents
    D: bathymetry of tide model
    constituents: list of model constituents
    """

    #-- read the netcdf format tide grid file
    #-- reading a combined global solution with localized solutions
    if GZIP:
        #-- open remote file with netCDF4
        #-- read GZIP file
        f = gzip.open(os.path.join(directory,grid_file),'rb')
        fileID = netCDF4.Dataset(grid_file,'r',memory=f.read())
    else:
        fileID = netCDF4.Dataset(os.path.join(directory,grid_file),'r')
    #-- variable dimensions
    nx = fileID.dimensions['nx'].size
    ny = fileID.dimensions['ny'].size
    #-- allocate numpy masked array for bathymetry
    bathymetry = np.ma.zeros((ny,nx))
    #-- read bathmetry and coordinates for variable type
    if (TYPE == 'z'):
        #-- get bathymetry at nodes
        bathymetry.data[:,:] = fileID.variables['hz'][:,:].T
        #-- read latitude and longitude at z-nodes
        lon = fileID.variables['lon_z'][:].copy()
        lat = fileID.variables['lat_z'][:].copy()
    elif TYPE in ('U','u'):
        #-- get bathymetry at u-nodes
        bathymetry.data[:,:] = fileID.variables['hu'][:,:].T
        #-- read latitude and longitude at u-nodes
        lon = fileID.variables['lon_u'][:].copy()
        lat = fileID.variables['lat_u'][:].copy()
    elif TYPE in ('V','v'):
        #-- get bathymetry at v-nodes
        bathymetry.data[:,:] = fileID.variables['hv'][:,:].T
        #-- read latitude and longitude at v-nodes
        lon = fileID.variables['lon_v'][:].copy()
        lat = fileID.variables['lat_v'][:].copy()
    #-- close the grid file
    fileID.close()
    f.close() if GZIP else None

    #-- grid step size of tide model
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]
    #-- replace original values with extend arrays/matrices
    lon = extend_array(lon, dlon)
    bathymetry = extend_matrix(bathymetry)
    #-- create masks
    bathymetry.mask = (bathymetry.data == 0)
    #-- create meshes from latitude and longitude
    gridlon,gridlat = np.meshgrid(lon,lat)

    #-- adjust longitudinal convention of input latitude and longitude
    #-- to fit tide model convention
    lt0, = np.nonzero(ilon < 0)
    ilon[lt0] += 360.0
    #-- number of points
    npts = len(ilon)

    #-- interpolate bathymetry and mask to output points
    D = np.ma.zeros((npts))
    D.mask = np.zeros((npts),dtype=np.bool)
    if (METHOD == 'bilinear'):
        #-- use quick bilinear to interpolate values
        D.data[:] = bilinear_interp(lon,lat,bathymetry.data,ilon,ilat)
        D.mask[:] = bilinear_interp(lon,lat,bathymetry.mask,ilon,ilat)
    elif (METHOD == 'spline'):
        #-- use scipy bivariate splines to interpolate values
        f1=scipy.interpolate.RectBivariateSpline(lon,lat,
            bathymetry.data.T,kx=1,ky=1)
        f2=scipy.interpolate.RectBivariateSpline(lon,lat,
            bathymetry.mask.T,kx=1,ky=1)
        D.data[:] = f1.ev(ilon,ilat)
        D.mask[:] = f2.ev(ilon,ilat).astype(np.bool)
    else:
        #-- use scipy griddata to interpolate values
        interp_points = np.c_[lon.flatten(),lat.flatten()]
        D.data[:] = scipy.interpolate.griddata(interp_points,
            bathymetry.flatten(), np.c_[ilon,ilat], method=METHOD)
        D.mask[:] = scipy.interpolate.griddata(interp_points,
            mask.flatten(), np.c_[ilon,ilat], method=METHOD).astype(np.bool)

    #-- u and v are velocities in cm/s
    if TYPE in ('v','u'):
        unit_conv = (D.data*100.0)
    #-- U and V are transports in m^2/s
    elif TYPE in ('V','U'):
        unit_conv = 1.0

    #-- number of constituents
    nc = len(model_files)
    #-- list of constituents
    constituents = []
    #-- amplitude and phase
    ampl = np.ma.zeros((npts,nc))
    ampl.mask = np.zeros((npts,nc),dtype=np.bool)
    phase = np.ma.zeros((npts,nc))
    phase.mask = np.zeros((npts,nc),dtype=np.bool)
    #-- read and interpolate each constituent
    for i,fi in enumerate(model_files):
        if (TYPE == 'z'):
            #-- read constituent from elevation file
            z,con = read_elevation_file(os.path.join(directory,fi),GZIP)
            #-- append constituent to list
            constituents.append(con)
            #-- replace original values with extend matrices
            z = extend_matrix(z)
            #-- interpolate amplitude and phase of the constituent
            z1 = np.ma.zeros((npts),dtype=z.dtype)
            z1.mask = np.zeros((npts),dtype=np.bool)
            if (METHOD == 'bilinear'):
                z1.data[:] = bilinear_interp(lon,lat,z,ilon,ilat,dtype=z.dtype)
                #-- mask zero values
                z1.mask[z1.data == 0] = True
                z1.data[z1.mask] = z1.fill_value
            elif (METHOD == 'spline'):
                f1 = scipy.interpolate.RectBivariateSpline(lon,lat,
                    z.data.real.T,kx=1,ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(lon,lat,
                    z.data.imag.T,kx=1,ky=1)
                f3 = scipy.interpolate.RectBivariateSpline(lon,lat,
                    z.mask.T,kx=1,ky=1)
                z1.data.real = f1.ev(ilon,ilat)
                z1.data.imag = f2.ev(ilon,ilat)
                z1.mask = f3.ev(ilon,ilat).astype(np.bool)
                #-- mask invalid values
                z1.data[z1.mask] = z1.fill_value
            else:
                #-- replace zero values with nan
                z[z==0] = np.nan
                #-- use scipy griddata to interpolate values
                z1.data = scipy.interpolate.griddata(interp_points,
                    z.flatten(), np.c_[ilon,ilat], method=METHOD)
                z1.mask[np.isnan(z1.data)] = True
                z1.data[z1.mask] = z1.fill_value
            #-- amplitude and phase of the constituent
            ampl[:,i] = np.abs(z1)
            phase[:,i] = np.arctan2(-np.imag(z1),np.real(z1))
        elif TYPE in ('U','u','V','v'):
            #-- read constituent from transport file
            tr,con = read_transport_file(os.path.join(directory,fi),TYPE,GZIP)
            #-- append constituent to list
            constituents.append(con)
            #-- replace original values with extend matrices
            tr = extend_matrix(tr)
            #-- interpolate values
            tr1 = np.ma.zeros((npts),dtype=tr.dtype)
            tr1.mask = np.zeros((npts),dtype=np.bool)
            if (METHOD == 'bilinear'):
                tr1.data[:]=bilinear_interp(lon,lat,tr,ilon,ilat,dtype=tr.dtype)
                #-- mask zero values
                tr1.mask[tr1.data == 0] = True
                tr1.data[tr1.mask] = tr1.fill_value
            elif (METHOD == 'spline'):
                f1 = scipy.interpolate.RectBivariateSpline(lon,lat,
                    tr.data.real.T,kx=1,ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(lon,lat,
                    tr.data.imag.T,kx=1,ky=1)
                f3 = scipy.interpolate.RectBivariateSpline(lon,lat,
                    tr.mask.T,kx=1,ky=1)
                tr1.data.real = f1.ev(ilon,ilat)
                tr1.data.imag = f2.ev(ilon,ilat)
                tr1.mask = f3.ev(ilon,ilat).astype(np.bool)
                #-- mask invalid values
                tr1.data[tr1.mask] = z1.fill_value
            else:
                #-- replace zero values with nan
                tr[tr==0] = np.nan
                #-- use scipy interpolate to interpolate values
                tr1 = scipy.interpolate.griddata(interp_points,
                    tr.flatten(), np.c_[ilon,ilat], method=METHOD)
                tr1.mask[np.isnan(tr1.data)] = True
                tr1.data[tr1.mask] = tr1.fill_value
            #-- convert units
            tr1 = tr1/unit_conv
            #-- amplitude and phase of the constituent
            ampl[:,i] = np.abs(tr1)
            phase[:,i] = np.arctan2(-np.imag(tr1),np.real(tr1))

    #-- convert amplitude from input units to meters
    amplitude = ampl*SCALE
    #-- convert phase to degrees
    phase = phase*180.0/np.pi
    phase[phase < 0] += 360.0
    #-- return the interpolated values
    return (amplitude,phase,D,constituents)

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
    temp = np.zeros((n+2),dtype=input_array.dtype)
    #-- extended array [x-1,x0,...,xN,xN+1]
    temp[0] = input_array[0] - step_size
    temp[1:-1] = input_array[:]
    temp[-1] = input_array[-1] + step_size
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
    temp = np.ma.zeros((ny,nx+2),dtype=input_matrix.dtype)
    temp[:,0] = input_matrix[:,-1]
    temp[:,1:-1] = input_matrix[:,:]
    temp[:,-1] = input_matrix[:,0]
    return temp

#-- PURPOSE: read elevation file to extract real and imaginary components for
#-- constituent
def read_elevation_file(input_file,GZIP):
    """
    Read elevation file to extract real and imaginary components for constituent

    Arguments
    ---------
    input_file: input elevation file

    Keyword arguments
    -----------------
    GZIP: input netCDF4 files are compressed

    Returns
    -------
    h: tidal elevation
    con: tidal constituent ID
    """
    #-- read the netcdf format tide elevation file
    #-- reading a combined global solution with localized solutions
    if GZIP:
        f = gzip.open(input_file,'rb')
        fileID = netCDF4.Dataset(input_file,'r',memory=f.read())
    else:
        fileID = netCDF4.Dataset(input_file,'r')
    #-- constituent name
    con = fileID.variables['con'][:].tostring().decode('utf-8')
    #-- variable dimensions
    nx = fileID.dimensions['nx'].size
    ny = fileID.dimensions['ny'].size
    #-- real and imaginary components of elevation
    h = np.ma.zeros((ny,nx),dtype=np.complex64)
    h.mask = np.zeros((ny,nx),dtype=np.bool)
    h.data.real[:,:] = fileID.variables['hRe'][:,:].T
    h.data.imag[:,:] = fileID.variables['hIm'][:,:].T
    #-- close the file
    fileID.close()
    f.close() if GZIP else None
    #-- return the elevation and constituent
    return (h,con.strip())

#-- PURPOSE: read transport file to extract real and imaginary components for
#-- constituent
def read_transport_file(input_file,TYPE,GZIP):
    """
    Read transport file to extract real and imaginary components for constituent

    Arguments
    ---------
    input_file: input transport file

    Keyword arguments
    -----------------
    TYPE: tidal variable to run
        u: horizontal transport velocities
        U: horizontal depth-averaged transport
        v: vertical transport velocities
        V: vertical depth-averaged transport
    GZIP: input netCDF4 files are compressed

    Returns
    -------
    tr: tidal transport
    con: tidal constituent ID
    """
    #-- read the netcdf format tide grid file
    #-- reading a combined global solution with localized solutions
    if GZIP:
        f = gzip.open(input_file,'rb')
        fileID = netCDF4.Dataset(input_file,'r',memory=f.read())
    else:
        fileID = netCDF4.Dataset(input_file,'r')
    #-- constituent name
    con = fileID.variables['con'][:].tostring().decode('utf-8')
    #-- variable dimensions
    nx = fileID.dimensions['nx'].size
    ny = fileID.dimensions['ny'].size
    #-- real and imaginary components of transport
    tr = np.ma.zeros((ny,nx),dtype=np.complex64)
    tr.mask = np.zeros((ny,nx),dtype=np.bool)
    if TYPE in ('U','u'):
        tr.data.real[:,:] = fileID.variables['uRe'][:,:].T
        tr.data.imag[:,:] = fileID.variables['uIm'][:,:].T
    elif TYPE in ('V','v'):
        tr.data.real[:,:] = fileID.variables['vRe'][:,:].T
        tr.data.imag[:,:] = fileID.variables['vIm'][:,:].T
    #-- close the file
    fileID.close()
    f.close() if GZIP else None
    #-- return the transport components and constituent
    return (tr,con.strip())
