#!/usr/bin/env python
u"""
read_netcdf_model.py (09/2019)
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
    type: tidal variable to run
        z: heights
        u: horizontal transport velocities
        v: vertical transport velocities

OPTIONS:
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

UPDATE HISTORY:
    Written 09/2019
"""
import os
import gzip
import netCDF4
import numpy as np
import scipy.interpolate

#-- extract tidal harmonic constants from OTIS tide models at coordinates
def extract_netcdf_constants(ilon, ilat, directory, grid_file, model_files,
    type, METHOD='spline', GZIP=True, SCALE=1):

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
    if (type == 'z'):
        #-- get bathymetry at nodes
        bathymetry.data[:,:] = fileID.variables['hz'][:,:].T
        #-- read latitude and longitude at z-nodes
        lon = fileID.variables['lon_z'][:].copy()
        lat = fileID.variables['lat_z'][:].copy()
    elif type in ('U','u'):
        #-- get bathymetry at u-nodes
        bathymetry.data[:,:] = fileID.variables['hu'][:,:].T
        #-- read latitude and longitude at u-nodes
        lon = fileID.variables['lon_u'][:].copy()
        lat = fileID.variables['lat_u'][:].copy()
    elif type in ('V','v'):
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
        D.data[:] = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
            bathymetry.flatten(), zip(x,y), method=METHOD)
        D.mask[:] = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
            mask.flatten(), zip(x,y), method=METHOD).astype(np.bool)

    #-- u and v are velocities in cm/s
    if type in ('v','u'):
        unit_conv = (D.data*100.0)
    #-- U and V are transports in m^2/s
    elif type in ('V','U'):
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
        if (type == 'z'):
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
                z1.data[:] = bilinear_interp(lon,lat,z,ilon,ilat)
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
                z1.data = scipy.interpolate.griddata(zip(lon.flatten(),
                    lon.flatten()), z.flatten(), zip(ilon,ilat),
                    method=METHOD)
                z1.mask[np.isnan(z1.data)] = True
                z1.data[z1.mask] = z1.fill_value
            #-- amplitude and phase of the constituent
            ampl[:,i] = np.abs(z1)
            phase[:,i] = np.arctan2(-np.imag(z1),np.real(z1))
        elif type in ('U','u','V','v'):
            #-- read constituent from transport file
            tr,con = read_transport_file(os.path.join(directory,fi),type,GZIP)
            #-- append constituent to list
            constituents.append(con)
            #-- replace original values with extend matrices
            tr = extend_matrix(tr)
            #-- interpolate values
            tr1 = np.ma.zeros((npts),dtype=tr.dtype)
            tr1.mask = np.zeros((npts),dtype=np.bool)
            if (METHOD == 'bilinear'):
                tr1.data[:] = bilinear_interp(lon,lat,tr,ilon,ilat)
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
                tr1 = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
                    tr.flatten(), zip(x,y), method=METHOD)
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

#-- wrapper function to extend an array
def extend_array(input_array,step_size):
    n = len(input_array)
    temp = np.zeros((n+2),dtype=input_array.dtype)
    temp[0] = input_array[0] - step_size
    temp[1:-1] = input_array[:]
    temp[-1] = input_array[-1] + step_size
    return temp

#-- wrapper function to extend a matrix
def extend_matrix(input_matrix):
    ny,nx = np.shape(input_matrix)
    temp = np.ma.zeros((ny,nx+2),dtype=input_matrix.dtype)
    temp[:,0] = input_matrix[:,-1]
    temp[:,1:-1] = input_matrix[:,:]
    temp[:,-1] = input_matrix[:,0]
    return temp

#-- read elevation file to extract real and imaginary components for constituent
def read_elevation_file(input_file,GZIP):
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

#-- read transport file to extract real and imaginary components for constituent
def read_transport_file(input_file,type,GZIP):
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
    if type in ('U','u'):
        tr.data.real[:,:] = fileID.variables['uRe'][:,:].T
        tr.data.imag[:,:] = fileID.variables['uIm'][:,:].T
    elif type in ('V','v'):
        tr.data.real[:,:] = fileID.variables['vRe'][:,:].T
        tr.data.imag[:,:] = fileID.variables['vIm'][:,:].T
    #-- close the file
    fileID.close()
    f.close() if GZIP else None
    #-- return the transport components and constituent
    return (tr,con.strip())

#-- PURPOSE: bilinear interpolation of input data to output data
def bilinear_interp(ilon,ilat,idata,lon,lat):
    #-- degrees to radians
    dtr = np.pi/180.0
    #-- grid step size of tide model
    dlon = np.abs(ilon[1] - ilon[0])
    dlat = np.abs(ilat[1] - ilat[0])
    #-- Convert input coordinates to radians
    phi = ilon*dtr
    th = (90.0 - ilat)*dtr
    #-- Convert output data coordinates to radians
    xphi = lon*dtr
    xth = (90.0 - lat)*dtr
    #-- interpolate gridded data values to data
    data = np.zeros_like(lon,dtype=np.complex128)
    for i,l in enumerate(lon):
        #-- calculating the indices for the original grid
        dx = (ilon - np.floor(lon[i]/dlon)*dlon)**2
        dy = (ilat - np.floor(lat[i]/dlat)*dlat)**2
        iph = np.min(np.nonzero(dx == np.min(dx)))
        ith = np.min(np.nonzero(dy == np.min(dy)))
        #-- if on corner value: use exact
        if ((lat[i] == ilat[ith]) & (lon[i] == ilon[iph])):
            data[i] = idata[ith,iph]
        elif ((lat[i] == ilat[ith+1]) & (lon[i] == ilon[iph])):
            data[i] = idata[ith+1,iph]
        elif ((lat[i] == ilat[ith]) & (lon[i] == ilon[iph+1])):
            data[i] = idata[ith,iph+1]
        elif ((lat[i] == ilat[ith+1]) & (lon[i] == ilon[iph+1])):
            data[i] = idata[ith+1,iph+1]
        else:
            #-- corner weight values for i,j
            Wa = (xphi[i]-phi[iph])*(xth[i]-th[ith])
            Wb = (phi[iph+1]-xphi[i])*(xth[i]-th[ith])
            Wc = (xphi[i]-phi[iph])*(th[ith+1]-xth[i])
            Wd = (phi[iph+1]-xphi[i])*(th[ith+1]-xth[i])
            #-- divisor weight value
            W = (phi[iph+1]-phi[iph])*(th[ith+1]-th[ith])
            #-- corner data values for i,j
            Ia = idata[ith,iph]#-- (0,0)
            Ib = idata[ith,iph+1]#-- (1,0)
            Ic = idata[ith+1,iph]#-- (0,1)
            Id = idata[ith+1,iph+1]#-- (1,1)
            #-- calculate interpolated value for i
            data[i] = (Ia*Wa + Ib*Wb + Ic*Wc + Id*Wd)/W
    #-- return interpolated values
    return data
