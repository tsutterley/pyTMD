#!/usr/bin/env python
u"""
read_tide_model.py (07/2020)
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
    GRID: binary file type to read
        ATLAS: reading a global solution with localized solutions
        OTIS: combined global solution

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

PROGRAM DEPENDENCIES:
    convert_ll_xy.py: converts lat/lon points to and from projected coordinates
    bilinear_interp.py: bilinear interpolation of data to specified coordinates

UPDATE HISTORY:
    Updated 07/2020: added function docstrings. separate bilinear interpolation
        update griddata interpolation. changed TYPE variable to keyword argument
    Updated 06/2020: output currents as numpy masked arrays
        use argmin and argmax in bilinear interpolation
    Updated 11/2019: interpolate heights and fluxes to numpy masked arrays
    Updated 09/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 01/2019: decode constituents for Python3 compatibility
    Updated 08/2018: added option GRID for using ATLAS outputs that
        combine both global and localized tidal solutions
        added multivariate spline interpolation option
    Updated 07/2018: added different interpolation methods
    Updated 09/2017: Adapted for Python
"""
import os
import numpy as np
import scipy.interpolate
from pyTMD.convert_ll_xy import convert_ll_xy
from pyTMD.bilinear_interp import bilinear_interp

#-- PURPOSE: extract tidal harmonic constants from tide models at coordinates
def extract_tidal_constants(ilon, ilat, grid_file, model_file, EPSG, TYPE='z',
    METHOD='spline', GRID='OTIS'):
    """
    Reads files for an OTIS-formatted tidal model
    Makes initial calculations to run the tide program
    Spatially interpolates tidal constituents to input coordinates

    Arguments
    ---------
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    grid_file: grid file for model
    model_file: model file containing each constituent
    EPSG: projection of tide model data
    TYPE: tidal variable to run
        z: heights
        u: horizontal transport velocities
        v: vertical transport velocities

    Keyword arguments
    -----------------
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, cubic, nearest: scipy griddata interpolations
    GRID: binary file type to read
        ATLAS: reading a global solution with localized solutions
        OTIS: combined global solution

    Returns
    -------
    amplitude: amplitudes of tidal constituents
    phase: phases of tidal constituents
    D: bathymetry of tide model
    constituents: list of model constituents
    """

    #-- read the OTIS-format tide grid file
    if (GRID == 'ATLAS'):
        #-- if reading a global solution with localized solutions
        x0,y0,hz0,mz0,iob,dt,pmask,local = read_atlas_grid(grid_file)
        xi,yi,hz = combine_atlas_model(x0,y0,hz0,pmask,local,VARIABLE='depth')
        mz = create_atlas_mask(x0,y0,mz0,local,VARIABLE='depth')
    else:
        #-- if reading a pure global solution
        xi,yi,hz,mz,iob,dt = read_tide_grid(grid_file)
    #-- run wrapper function to convert coordinate systems of input lat/lon
    x,y = convert_ll_xy(ilon,ilat,EPSG,'F')
    #-- grid step size of tide model
    dx = xi[1] - xi[0]
    dy = yi[1] - yi[0]

    if (TYPE != 'z'):
        mz,mu,mv = Muv(hz)
        hu,hv = Huv(hz)

    #-- if global: extend limits
    GLOBAL = False
    #-- replace original values with extend arrays/matrices
    if ((xi[-1] - xi[0]) == (360.0 - dx)) & (EPSG == '4326'):
        xi = extend_array(xi, dx)
        hz = extend_matrix(hz)
        mz = extend_matrix(mz)
        #-- set global flag
        GLOBAL = True

    #-- adjust longitudinal convention of input latitude and longitude
    #-- to fit tide model convention
    xmin = np.min(x)
    xmax = np.max(y)
    if (xmin < xi[0]) & (EPSG == '4326'):
        lt0, = np.nonzero(x < 0)
        x[lt0] += 360.0
    if (xmax > xi[-1]) & (EPSG == '4326'):
        gt180, = np.nonzero(x > 180)
        x[gt180] -= 360.0

    #-- create meshes from latitude and longitude
    ux = xi - dx/2.0
    vy = yi - dy/2.0
    X,Y = np.meshgrid(xi,yi)
    Xu,Yu = np.meshgrid(ux,yi)
    Xv,Yv = np.meshgrid(xi,vy)

    #-- masks zero values
    hz = np.ma.array(hz,mask=(hz==0))
    if (TYPE != 'z'):
        #-- replace original values with extend matrices
        if GLOBAL:
            hu = extend_matrix(hu)
            hv = extend_matrix(hv)
            mu = extend_matrix(mu)
            mv = extend_matrix(mv)
        #-- masks zero values
        hu = np.ma.array(hu,mask=(hu==0))
        hv = np.ma.array(hv,mask=(hv==0))

    #-- interpolate depth and mask to output points
    if (METHOD == 'bilinear'):
        #-- use quick bilinear to interpolate values
        D = bilinear_interp(xi,yi,hz,x,y)
        mz1 = bilinear_interp(xi,yi,mz,x,y)
        mz1 = np.ceil(mz1).astype(mz.dtype)
        if (TYPE != 'z'):
            mu1 = bilinear_interp(xi,yi,mu,x,y)
            mu1 = np.ceil(mu1).astype(mu.dtype)
            mv1 = bilinear_interp(xi,yi,mv,x,y)
            mv1 = np.ceil(mv1).astype(mz.dtype)
    elif (METHOD == 'spline'):
        #-- use scipy bivariate splines to interpolate values
        f1=scipy.interpolate.RectBivariateSpline(xi,yi,hz.T,kx=1,ky=1)
        f2=scipy.interpolate.RectBivariateSpline(xi,yi,mz.T,kx=1,ky=1)
        D = f1.ev(x,y)
        mz1 = np.ceil(f2.ev(x,y)).astype(mz.dtype)
        if (TYPE != 'z'):
            f3=scipy.interpolate.RectBivariateSpline(xi,yi,mu.T,kx=1,ky=1)
            f4=scipy.interpolate.RectBivariateSpline(xi,yi,mv.T,kx=1,ky=1)
            mu1 = np.ceil(f3.ev(x,y)).astype(mu.dtype)
            mv1 = np.ceil(f4.ev(x,y)).astype(mv.dtype)
    else:
        #-- use scipy griddata to interpolate values
        interp_points = np.c_[X.flatten(),Y.flatten()]
        D = scipy.interpolate.griddata(interp_points,
            hz.flatten(), np.c_[x,y], method=METHOD)
        mz1 = scipy.interpolate.griddata(interp_points,
            mz.real.flatten(), np.c_[x,y], method=METHOD)
        mz1 = np.ceil(mz1).astype(mz.dtype)
        if (TYPE != 'z'):
            mu1 = scipy.interpolate.griddata(interp_points,
                mu.flatten(), np.c_[x,y], method=METHOD)
            mv1 = scipy.interpolate.griddata(interp_points,
                mv.flatten(), np.c_[x,y], method=METHOD)
            mu1 = np.ceil(mu1).astype(mu.dtype)
            mv1 = np.ceil(mv1).astype(mv.dtype)

    #-- u and v are velocities in cm/s
    if TYPE in ('v','u'):
        unit_conv = (D*100.0)
    #-- U and V are transports in m^2/s
    elif TYPE in ('V','U'):
        unit_conv = 1.0

    #-- read and interpolate each constituent
    constituents,nc = read_constituents(model_file)
    npts = len(D)
    amplitude = np.ma.zeros((npts,nc))
    amplitude.mask = np.zeros((npts,nc),dtype=np.bool)
    ph = np.ma.zeros((npts,nc))
    ph.mask = np.zeros((npts,nc),dtype=np.bool)
    for i,c in enumerate(constituents):
        if (TYPE == 'z'):
            #-- read constituent from elevation file
            if (GRID == 'ATLAS'):
                z0,zlocal = read_atlas_elevation(model_file,i,c)
                xi,yi,z=combine_atlas_model(x0,y0,z0,pmask,zlocal,VARIABLE='z')
            else:
                z = read_elevation_file(model_file,i)
            #-- replace original values with extend matrices
            if GLOBAL:
                z = extend_matrix(z)
            #-- interpolate amplitude and phase of the constituent
            if (METHOD == 'bilinear'):
                #-- replace zero values with nan
                z[z==0] = np.nan
                #-- use quick bilinear to interpolate values
                z1 = np.ma.zeros((npts),dtype=z.dtype)
                z1.data = bilinear_interp(xi,yi,z,x,y,dtype=np.complex128)
                #-- replace nan values with fill_value
                z1.mask = (np.isnan(z1.data) | (~mz1.astype(np.bool)))
                z1.data[z1.mask] = z1.fill_value
            elif (METHOD == 'spline'):
                #-- use scipy bivariate splines to interpolate values
                f1 = scipy.interpolate.RectBivariateSpline(xi,yi,
                    z.real.T,kx=1,ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(xi,yi,
                    z.imag.T,kx=1,ky=1)
                z1 = np.ma.zeros((npts),dtype=z.dtype)
                z1.data.real = f1.ev(x,y)
                z1.data.imag = f2.ev(x,y)
                #-- replace zero values with fill_value
                z1.mask = (~mz1.astype(np.bool))
                z1.data[z1.mask] = z1.fill_value
            else:
                #-- replace zero values with nan
                z[z==0] = np.nan
                #-- use scipy griddata to interpolate values
                z1 = np.ma.zeros((npts),dtype=z.dtype)
                z1.data=scipy.interpolate.griddata(interp_points,
                    z.flatten(), np.c_[x,y], method=METHOD)
                #-- replace nan values with fill_value
                z1.mask = (np.isnan(z1.data) | (~mz1.astype(np.bool)))
                z1.data[z1.mask] = z1.fill_value
            #-- amplitude and phase of the constituent
            amplitude.data[:,i] = np.abs(z1.data)
            amplitude.mask[:,i] = np.copy(z1.mask)
            ph.data[:,i] = np.arctan2(-np.imag(z1.data),np.real(z1.data))
            ph.mask[:,i] = np.copy(z1.mask)
        elif TYPE in ('U','u'):
            #-- read constituent from transport file
            if (GRID == 'ATLAS'):
                u0,v0,uvlocal = read_atlas_transport(model_file,i,c)
                xi,yi,u=combine_atlas_model(x0,y0,u0,pmask,uvlocal,VARIABLE='u')
            else:
                u,v = read_transport_file(model_file,i)
            #-- replace original values with extend matrices
            if GLOBAL:
                u = extend_matrix(u)
            #-- interpolate values
            if (METHOD == 'bilinear'):
                #-- replace zero values with nan
                u[u==0] = np.nan
                #-- use quick bilinear to interpolate values
                u1 = np.ma.zeros((npts),dtype=u.dtype)
                u1.data = bilinear_interp(xi,yi,u,x,y,dtype=np.complex128)
                #-- replace nan values with fill_value
                u1.mask = (np.isnan(u1.data) | (~mu1.astype(np.bool)))
                u1.data[u1.mask] = u1.fill_value
            elif (METHOD == 'spline'):
                f1 = scipy.interpolate.RectBivariateSpline(xi,yi,
                    u.real.T,kx=1,ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(xi,yi,
                    u.imag.T,kx=1,ky=1)
                u1 = np.ma.zeros((npts),dtype=u.dtype)
                u1.data.real = f1.ev(x,y)
                u1.data.imag = f2.ev(x,y)
                #-- replace zero values with fill_value
                u1.mask = (~mu1.astype(np.bool))
                u1.data[u1.mask] = u1.fill_value
            else:
                #-- replace zero values with nan
                u[u==0] = np.nan
                #-- use scipy griddata to interpolate values
                u1 = np.ma.zeros((npts),dtype=u.dtype)
                u1.data=scipy.interpolate.griddata(interp_points,
                    u.flatten(), np.c_[x,y], method=METHOD)
                #-- replace nan values with fill_value
                u1.mask = (np.isnan(u1.data) | (~mu1.astype(np.bool)))
                u1.data[u1.mask] = u1.fill_value
            #-- convert units
            u1 = u1/unit_conv
            #-- amplitude and phase of the constituent
            amplitude.data[:,i] = np.abs(u1)
            amplitude.mask[:,i] = np.copy(u1.mask)
            ph.data[:,i] = np.arctan2(-np.imag(u1),np.real(u1))
            ph.mask[:,i] = np.copy(u1.mask)
        elif TYPE in ('V','v'):
            #-- read constituent from transport file
            if (GRID == 'ATLAS'):
                u0,v0,uvlocal = read_atlas_transport(model_file,i,c)
                xi,yi,v = combine_atlas_model(x0,y0,v0,pmask,local,VARIABLE='v')
            else:
                u,v = read_transport_file(model_file,i)
            #-- replace original values with extend matrices
            if GLOBAL:
                v = extend_matrix(v)
            #-- interpolate values
            if (METHOD == 'bilinear'):
                #-- replace zero values with nan
                v[v==0] = np.nan
                #-- use quick bilinear to interpolate values
                v1 = np.ma.zeros((npts),dtype=v.dtype)
                v1.data = bilinear_interp(xi,yi,v,x,y,dtype=np.complex128)
                #-- replace nan values with fill_value
                v1.mask = (np.isnan(v1.data) | (~mv1.astype(np.bool)))
                v1.data[v1.mask] = v1.fill_value
            elif (METHOD == 'spline'):
                f1 = scipy.interpolate.RectBivariateSpline(xi,yi,
                    v.real.T,kx=1,ky=1)
                f2 = scipy.interpolate.RectBivariateSpline(xi,yi,
                    v.imag.T,kx=1,ky=1)
                v1 = np.ma.zeros((npts),dtype=v.dtype)
                v1.data.real = f1.ev(x,y)
                v1.data.imag = f2.ev(x,y)
                #-- replace zero values with fill_value
                v1.mask = (~mv1.astype(np.bool))
                v1.data[v1.mask] = v1.fill_value
            else:
                #-- replace zero values with nan
                v[u==0] = np.nan
                #-- use scipy griddata to interpolate values
                v1 = np.ma.zeros((npts),dtype=v.dtype)
                v1.data=scipy.interpolate.griddata(interp_points,
                    v.flatten(), np.c_[x,y], method=METHOD)
                #-- replace nan values with fill_value
                v1.mask = (np.isnan(v1.data) | (~mv1.astype(np.bool)))
                v1.data[v1.mask] = v1.fill_value
            #-- convert units
            v1 = v1/unit_conv
            #-- amplitude and phase of the constituent
            amplitude.data[:,i] = np.abs(v1)
            amplitude.mask[:,i] = np.copy(v1.mask)
            ph.data[:,i] = np.arctan2(-np.imag(v1),np.real(v1))
            ph.mask[:,i] = np.copy(v1.mask)

    #-- convert phase to degrees
    phase = ph*180.0/np.pi
    phase.data[phase.data < 0] += 360.0
    #-- replace data for invalid mask values
    amplitude.data[amplitude.mask] = amplitude.fill_value
    phase.data[phase.mask] = phase.fill_value
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

#-- PURPOSE: read tide grid file
def read_tide_grid(input_file):
    """
    Read grid file to extract model coordinates, bathymety, masks and indices

    Arguments
    ---------
    input_file: input grid file

    Returns
    -------
    x: x-coordinates of input grid
    y: y-coordinates of input grid
    hz: model bathymety
    mz: land/water mask
    iob: open boundary index
    dt: time step
    """
    #-- open the file
    fid = open(os.path.expanduser(input_file),'rb')
    fid.seek(4,0)
    #-- read data as big endian
    #-- get model dimensions and limits
    nx, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    ny, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    #-- extract x and y limits (these could be latitude and longitude)
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    dt, = np.fromfile(fid, dtype=np.dtype('>f4'), count=1)
    #-- convert longitudinal limits (if x == longitude)
    if (xlim[0] < 0) & (xlim[1] < 0) & (dt > 0):
        xlim += 360.0
    #-- create x and y arrays arrays (these could be lon and lat values)
    dx = (xlim[1] - xlim[0])/nx
    dy = (ylim[1] - ylim[0])/ny
    x = np.linspace(xlim[0]+dx/2.0,xlim[1]-dx/2.0,nx)
    y = np.linspace(ylim[0]+dy/2.0,ylim[1]-dy/2.0,ny)
    #-- read nob and iob from file
    nob, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    if (nob == 0):
        fid.seek(20,1)
        iob = []
    else:
        fid.seek(8,1)
        iob=np.fromfile(fid, dtype=np.dtype('>i4'), count=2*nob).reshape(nob,2)
        fid.seek(8,1)
    #-- read hz matrix
    hz = np.fromfile(fid, dtype=np.dtype('>f4'), count=nx*ny).reshape(ny,nx)
    fid.seek(8,1)
    #-- read mz matrix
    mz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nx*ny).reshape(ny,nx)
    #-- close the file
    fid.close()
    #-- return values
    return (x,y,hz,mz,iob,dt)

#-- PURPOSE: read tide grid file with localized solutions
def read_atlas_grid(input_file):
    """
    Read ATLAS grid file to extract model coordinates, bathymety, masks and
    indices for both global and local solutions

    Arguments
    ---------
    input_file: input ATLAS grid file

    Returns
    -------
    x: x-coordinates of input ATLAS grid
    y: y-coordinates of input ATLAS grid
    hz: model bathymety
    mz: land/water mask
    iob: open boundary index
    dt: time step
    pmask: global mask
    local: dictionary of local tidal solutions for grid variables
        depth: model bathymety
    """
    #-- read the input file to get file information
    fd = os.open(os.path.expanduser(input_file),os.O_RDONLY)
    file_info = os.fstat(fd)
    fid = os.fdopen(fd, 'rb')
    fid.seek(4,0)
    #-- read data as big endian
    #-- get model dimensions and limits
    nx, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    ny, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    #-- extract latitude and longitude limits
    lats = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    lons = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    dt, = np.fromfile(fid, dtype=np.dtype('>f4'), count=1)
    #-- create lon and lat arrays
    dlon = (lons[1] - lons[0])/nx
    dlat = (lats[1] - lats[0])/ny
    x = np.linspace(lons[0]+dlon/2.0,lons[1]-dlon/2.0,nx)
    y = np.linspace(lats[0]+dlat/2.0,lats[1]-dlat/2.0,ny)
    #-- read nob and iob from file
    nob, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    if (nob == 0):
        fid.seek(20,1)
        iob = []
    else:
        fid.seek(8,1)
        iob=np.fromfile(fid, dtype=np.dtype('>i4'), count=2*nob).reshape(nob,2)
        fid.seek(8,1)
    #-- read hz matrix
    hz = np.fromfile(fid, dtype=np.dtype('>f4'), count=nx*ny).reshape(ny,nx)
    fid.seek(8,1)
    #-- read mz matrix
    mz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nx*ny).reshape(ny,nx)
    fid.seek(8,1)
    #-- read pmask matrix
    pmask = np.fromfile(fid, dtype=np.dtype('>i4'), count=nx*ny).reshape(ny,nx)
    fid.seek(4,1)
    #-- read local models
    nmod = 0
    local = {}
    #-- while the file position is not at the end of file
    while (fid.tell() < file_info.st_size):
        #-- add 1 to number of models
        fid.seek(4,1)
        nmod += 1
        #-- get local model dimensions and limits
        nx1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        ny1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nd, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        #-- extract latitude and longitude limits of local model
        lt1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        ln1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        #-- extract name
        name = fid.read(20).strip()
        fid.seek(8,1)
        iz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nd)
        jz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nd)
        fid.seek(8,1)
        depth = np.full((ny1,nx1),np.nan)
        depth[jz-1,iz-1] = np.fromfile(fid, dtype=np.dtype('>f4'), count=nd)
        fid.seek(4,1)
        #-- save to dictionary
        local[name] = dict(lon=ln1,lat=lt1,depth=depth)
    #-- close the file
    fid.close()
    #-- return values
    return (x,y,hz,mz,iob,dt,pmask,local)

#-- PURPOSE: read list of constituents from an elevation or transport file
def read_constituents(input_file):
    """
    Read the list of constituents from an elevation or transport file

    Arguments
    ---------
    input_file: input tidal file

    Returns
    -------
    constituents: list of tidal constituent IDs
    nc: number of constituents
    """
    #-- open the file
    fid = open(os.path.expanduser(input_file),'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    fid.seek(16,1)
    constituents = [c.decode("utf-8").rstrip() for c in fid.read(nc*4).split()]
    fid.close()
    return (constituents,nc)

#-- PURPOSE: read elevation file to extract real and imaginary components for
#-- constituent
def read_elevation_file(input_file,ic):
    """
    Read elevation file to extract real and imaginary components for constituent

    Arguments
    ---------
    input_file: input elevation file
    ic: index of consituent

    Returns
    -------
    h: tidal elevation
    """
    #-- open the file
    fid = open(os.path.expanduser(input_file),'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    #-- extract x and y limits
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    #-- skip records to constituent
    nskip = ic*(nx*ny*8+8) + 8 + ll - 28
    fid.seek(nskip,1)
    #-- real and imaginary components of elevation
    h = np.ma.zeros((ny,nx),dtype=np.complex64)
    h.mask = np.zeros((ny,nx),dtype=np.bool)
    for i in range(ny):
        temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nx)
        h.data.real[i,:] = temp[0:2*nx-1:2]
        h.data.imag[i,:] = temp[1:2*nx:2]
    #-- close the file
    fid.close()
    #-- return the elevation
    return h

#-- PURPOSE: read elevation file with localized solutions to extract real and
#-- imaginary components for constituent
def read_atlas_elevation(input_file,ic,constituent):
    """
    Read elevation file with localized solutions to extract real and imaginary
    components for constituent

    Arguments
    ---------
    input_file: input ATLAS elevation file
    ic: index of consituent
    constituent: tidal constituent ID

    Returns
    -------
    h: global tidal elevation
    local: dictionary of local tidal solutions for elevation variables
        z: tidal elevation
    """
    #-- read the input file to get file information
    fd = os.open(os.path.expanduser(input_file),os.O_RDONLY)
    file_info = os.fstat(fd)
    fid = os.fdopen(fd, 'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    #-- extract x and y limits
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    #-- skip records to constituent
    nskip = 8 + nc*4 + ic*(nx*ny*8 + 8)
    fid.seek(nskip,1)
    #-- real and imaginary components of elevation
    h = np.ma.zeros((ny,nx),dtype=np.complex64)
    h.mask = np.zeros((ny,nx),dtype=np.bool)
    for i in range(ny):
        temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nx)
        h.data.real[i,:] = temp[0:2*nx-1:2]
        h.data.imag[i,:] = temp[1:2*nx:2]
    #-- skip records after constituent
    nskip = (nc-ic-1)*(nx*ny*8 + 8) + 4
    fid.seek(nskip,1)
    #-- read local models to find constituent
    nmod = 0
    local = {}
    #-- while the file position is not at the end of file
    while (fid.tell() < file_info.st_size):
        #-- add 1 to number of models
        fid.seek(4,1)
        nmod += 1
        #-- get local model dimensions and limits
        nx1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        ny1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nc1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nz, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        #-- extract latitude and longitude limits of local model
        lt1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        ln1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        #-- extract constituents for localized solution
        cons = fid.read(nc1*4).strip().split()
        #-- check if constituent is in list of localized solutions
        if (constituent in cons):
            ic1, = [i for i,c in enumerate(cons) if (c == constituent)]
            #-- extract name
            name = fid.read(20).strip()
            fid.seek(8,1)
            iz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nz)
            jz = np.fromfile(fid, dtype=np.dtype('>i4'), count=nz)
            #-- skip records to constituent
            nskip = 8 + ic1*(8*nz + 8)
            fid.seek(nskip,1)
            #-- real and imaginary components of elevation
            h1 = np.ma.zeros((ny1,nx1),fill_value=np.nan,dtype=np.complex64)
            h1.mask = np.zeros((ny1,nx1),dtype=np.bool)
            temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nz)
            h1.data.real[jz-1,iz-1] = temp[0:2*nz-1:2]
            h1.data.imag[jz-1,iz-1] = temp[1:2*nz:2]
            #-- save constituent to dictionary
            local[name] = dict(lon=ln1,lat=lt1,z=h1)
            #-- skip records after constituent
            nskip = (nc1-ic1-1)*(8*nz + 8) + 4
            fid.seek(nskip,1)
        else:
            #-- skip records for local model if constituent not in list
            nskip = 40 + 16*nz + (nc1-1)*(8*nz + 8)
            fid.seek(nskip,1)
    #-- close the file
    fid.close()
    #-- return the elevation
    return (h,local)

#-- PURPOSE: read transport file to extract real and imaginary components for
#-- constituent
def read_transport_file(input_file,ic):
    """
    Read transport file to extract real and imaginary components for constituent

    Arguments
    ---------
    input_file: input transport file
    ic: index of consituent

    Returns
    -------
    u: zonal tidal transport
    v: meridional zonal transport
    """
    #-- open the file
    fid = open(os.path.expanduser(input_file),'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    #-- extract x and y limits
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    #-- skip records to constituent
    nskip = ic*(nx*ny*16+8) + 8 + ll - 28
    fid.seek(nskip,1)
    #-- real and imaginary components of transport
    u = np.ma.zeros((ny,nx),dtype=np.complex64)
    u.mask = np.zeros((ny,nx),dtype=np.bool)
    v = np.ma.zeros((ny,nx),dtype=np.complex64)
    v.mask = np.zeros((ny,nx),dtype=np.bool)
    for i in range(ny):
        temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=4*nx)
        u.data.real[i,:] = temp[0:4*nx-3:4]
        u.data.imag[i,:] = temp[1:4*nx-2:4]
        v.data.real[i,:] = temp[2:4*nx-1:4]
        v.data.imag[i,:] = temp[3:4*nx:4]
    #-- close the file
    fid.close()
    #-- return the transport components
    return (u,v)

#-- PURPOSE: read transport file with localized solutions to extract real and
#-- imaginary components for constituent
def read_atlas_transport(input_file,ic,constituent):
    """
    Read transport file with localized solutions to extract real and imaginary
    components for constituent

    Arguments
    ---------
    input_file: input ATLAS transport file
    ic: index of consituent
    constituent: tidal constituent ID

    Returns
    -------
    u: global zonal tidal transport
    v: global meridional zonal transport
    local: dictionary of local tidal solutions for transport variables
        u: zonal tidal transport
        v: meridional zonal transport
    """
    #-- read the input file to get file information
    fd = os.open(os.path.expanduser(input_file),os.O_RDONLY)
    file_info = os.fstat(fd)
    fid = os.fdopen(fd, 'rb')
    ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
    nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
    #-- extract x and y limits
    ylim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    xlim = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
    #-- skip records to constituent
    nskip = 8 + nc*4 + ic*(nx*ny*16 + 8)
    fid.seek(nskip,1)
    #-- real and imaginary components of transport
    u = np.ma.zeros((ny,nx),dtype=np.complex64)
    u.mask = np.zeros((ny,nx),dtype=np.bool)
    v = np.ma.zeros((ny,nx),dtype=np.complex64)
    v.mask = np.zeros((ny,nx),dtype=np.bool)
    for i in range(ny):
        temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=4*nx)
        u.data.real[i,:] = temp[0:4*nx-3:4]
        u.data.imag[i,:] = temp[1:4*nx-2:4]
        v.data.real[i,:] = temp[2:4*nx-1:4]
        v.data.imag[i,:] = temp[3:4*nx:4]
    #-- skip records after constituent
    nskip = (nc-ic-1)*(nx*ny*16 + 8) + 4
    fid.seek(nskip,1)
    #-- read local models to find constituent
    nmod = 0
    local = {}
    #-- while the file position is not at the end of file
    while (fid.tell() < file_info.st_size):
        #-- add 1 to number of models
        fid.seek(4,1)
        nmod += 1
        #-- get local model dimensions and limits
        nx1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        ny1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nc1, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nu, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        nv, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
        #-- extract latitude and longitude limits of local model
        lt1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        ln1 = np.fromfile(fid, dtype=np.dtype('>f4'), count=2)
        #-- extract constituents for localized solution
        cons = fid.read(nc1*4).strip().split()
        #-- check if constituent is in list of localized solutions
        if (constituent in cons):
            ic1, = [i for i,c in enumerate(cons) if (c == constituent)]
            #-- extract name
            name = fid.read(20).strip()
            fid.seek(8,1)
            iu = np.fromfile(fid, dtype=np.dtype('>i4'), count=nu)
            ju = np.fromfile(fid, dtype=np.dtype('>i4'), count=nu)
            fid.seek(8,1)
            iv = np.fromfile(fid, dtype=np.dtype('>i4'), count=nv)
            jv = np.fromfile(fid, dtype=np.dtype('>i4'), count=nv)
            #-- skip records to constituent
            nskip = 8 + ic1*(8*nu + 8*nv + 16)
            fid.seek(nskip,1)
            #-- real and imaginary components of u transport
            u1 = np.ma.zeros((ny1,nx1),fill_value=np.nan,dtype=np.complex64)
            u1.mask = np.zeros((ny1,nx1),dtype=np.bool)
            tmpu = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nu)
            u1.data.real[ju-1,iu-1] = tmpu[0:2*nu-1:2]
            u1.data.imag[ju-1,iu-1] = tmpu[1:2*nu:2]
            fid.seek(8,1)
            #-- real and imaginary components of v transport
            v1 = np.ma.zeros((ny1,nx1),fill_value=np.nan,dtype=np.complex64)
            v1.mask = np.zeros((ny1,nx1),dtype=np.bool)
            tmpv = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nv)
            v1.data.real[jv-1,iv-1] = tmpv[0:2*nv-1:2]
            v1.data.imag[jv-1,iv-1] = tmpv[1:2*nv:2]
            #-- save constituent to dictionary
            local[name] = dict(lon=ln1,lat=lt1,u=u1,v=v1)
            #-- skip records after constituent
            nskip = (nc1-ic1-1)*(8*nu + 8*nv + 16) + 4
            fid.seek(nskip,1)
        else:
            #-- skip records for local model if constituent not in list
            nskip = 56 + 16*nu + 16*nv + (nc1-1)*(8*nu + 8*nv + 16)
            fid.seek(nskip,1)
    #-- close the file
    fid.close()
    #-- return the transport components
    return (u,v,local)

#-- PURPOSE: create a 2 arc-minute grid mask from mz and depth variables
def create_atlas_mask(xi,yi,mz,local,VARIABLE=None):
    """
    Creates a high-resolution grid mask from model variables

    Arguments
    ---------
    xi: input x-coordinates of global tide model
    yi: input y-coordinates of global tide model
    mz: global land/water mask
    local: dictionary of local tidal solutions

    Keyword arguments
    -----------------
    VARIABLE: key for variable within each local solution
        depth: model bathymety

    Returns
    -------
    x30: x-coordinates of high-resolution tide model
    y30: y-coordinates of high-resolution tide model
    m30: high-resolution land/water mask
    """
    #-- create 2 arc-minute grid dimensions
    d30 = 1.0/30.0
    x30 = np.arange(d30/2.0, 360.0+d30/2.0, d30)
    y30 = np.arange(-90.0+d30/2.0, 90.0+d30/2.0, d30)
    #-- interpolate global mask to create initial 2 arc-minute mask
    m30 = np.zeros((len(y30),len(x30)),dtype=mz.dtype)
    f = scipy.interpolate.RectBivariateSpline(xi,yi,mz.T,kx=1,ky=1)
    m30[:,:] = f(x30,y30).T
    #-- iterate over localized solutions to fill in high-resolution coastlines
    for key,val in local.items():
        #-- local model output
        zlocal = val[VARIABLE][:]
        validy,validx = np.nonzero(np.isfinite(zlocal.real))
        #-- create latitude and longitude for local model
        ilon = np.arange(val['lon'][0]+d30/2.0,val['lon'][1]+d30/2.0,d30)
        ilat = np.arange(val['lat'][0]+d30/2.0,val['lat'][1]+d30/2.0,d30)
        X,Y = np.meshgrid(ilon,ilat)
        for indy,indx in zip(validy,validx):
            #-- check if model is -180:180
            lon30 = (X[indy,indx]+360.) if (X[indy,indx]<=0.0) else X[indy,indx]
            ii = np.int((lon30 - x30[0])/d30)
            jj = np.int((Y[indy,indx] - y30[0])/d30)
            #-- fill global mask with regional solution
            m30[jj,ii] = 1
    #-- return the 2 arc-minute mask
    return m30

#-- PURPOSE: combines global and local atlas solutions
def combine_atlas_model(xi,yi,zi,pmask,local,VARIABLE=None):
    """
    Combines global and local ATLAS tidal solutions into a single
    high-resolution solution

    Arguments
    ---------
    xi: input x-coordinates of global tide model
    yi: input y-coordinates of global tide model
    zi: global tide model data
    pmask: global mask
    local: dictionary of local tidal solutions

    Keyword arguments
    -----------------
    VARIABLE: key for variable within each local solution
        depth: model bathymety
        z: tidal elevation
        u: zonal tidal transport
        v: meridional zonal transport

    Returns
    -------
    x30: x-coordinates of high-resolution tide model
    y30: y-coordinates of high-resolution tide model
    z30: combined high-resolution tidal solution for variable
    """
    #-- create 2 arc-minute grid dimensions
    d30 = 1.0/30.0
    x30 = np.arange(d30/2.0, 360.0+d30/2.0, d30)
    y30 = np.arange(-90.0+d30/2.0, 90.0+d30/2.0, d30)
    #-- interpolate global solution to 2 arc-minute solution
    z30 = np.ma.zeros((len(y30),len(x30)),dtype=zi.dtype)
    z30.mask = np.zeros((len(y30),len(x30)),dtype=np.bool)
    #-- test if combining elevation/transport variables with complex components
    if np.iscomplexobj(z30):
        f1 = scipy.interpolate.RectBivariateSpline(xi, yi, zi.real.T, kx=1,ky=1)
        f2 = scipy.interpolate.RectBivariateSpline(xi, yi, zi.imag.T, kx=1,ky=1)
        z30.data.real[:,:] = f1(x30,y30).T
        z30.data.imag[:,:] = f2(x30,y30).T
    else:
        f = scipy.interpolate.RectBivariateSpline(xi, yi, zi.T, kx=1,ky=1)
        z30.data[:,:] = f(x30,y30).T
    #-- iterate over localized solutions
    for key,val in local.items():
        #-- local model output
        zlocal = val[VARIABLE][:]
        validy,validx = np.nonzero(np.isfinite(zlocal.data.real))
        #-- create latitude and longitude for local model
        ilon = np.arange(val['lon'][0]+d30/2.0,val['lon'][1]+d30/2.0,d30)
        ilat = np.arange(val['lat'][0]+d30/2.0,val['lat'][1]+d30/2.0,d30)
        X,Y = np.meshgrid(ilon,ilat)
        for indy,indx in zip(validy,validx):
            #-- check if model is -180:180
            lon30 = (X[indy,indx]+360.) if (X[indy,indx]<=0.0) else X[indy,indx]
            ii = np.int((lon30 - x30[0])/d30)
            jj = np.int((Y[indy,indx] - y30[0])/d30)
            #-- fill global model with regional solution
            z30.data[jj,ii] = zlocal[indy,indx]
    #-- return 2 arc-minute solution and coordinates
    return (x30,y30,z30)

#-- For a rectangular bathymetry grid:
#-- construct masks for zeta, u and v nodes on a C-grid
def Muv(hz):
    """
    Construct masks for zeta, u and v nodes on a C-grid
    """
    ny,nx = np.shape(hz)
    mz = (hz > 0).astype(np.int)
    #-- x-indices
    indx = np.zeros((nx),dtype=np.int)
    indx[:-1] = np.arange(1,nx)
    indx[-1] = 0
    #-- y-indices
    indy = np.zeros((ny),dtype=np.int)
    indy[:-1] = np.arange(1,ny)
    indy[-1] = 0
    #-- calculate mu and mv
    mu = np.zeros((ny,nx),dtype=np.int)
    mv = np.zeros((ny,nx),dtype=np.int)
    mu[indy,:] = mz*mz[indy,:]
    mv[:,indx] = mz*mz[:,indx]
    return (mu,mv,mz)

#-- PURPOSE: Interpolate bathymety to zeta, u and v nodes on a C-grid
def Huv(hz):
    """
    Interpolate bathymety to zeta, u and v nodes on a C-grid
    """
    ny,nx = np.shape(hz)
    mu,mv,mz = Muv(hz)
    #-- x-indices
    indx = np.zeros((nx),dtype=np.int)
    indx[0] = nx-1
    indx[1:] = np.arange(1,nx)
    #-- y-indices
    indy = np.zeros((ny),dtype=np.int)
    indy[0] = ny-1
    indy[1:] = np.arange(1,ny)
    #-- calculate hu and hv
    hu = mu*(hz + hz[indy,:])/2.0
    hv = mv*(hz + hz[:,indx])/2.0
    return (hu,hv)
