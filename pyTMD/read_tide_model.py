#!/usr/bin/env python
u"""
read_tide_model.py (07/2018)
Reads files for a tidal model and makes initial calculations to run tide program
Includes functions to extract tidal harmonic constants out of a tidal model for
	given locations

Reads OTIS format tidal solutions provided by Ohio State University and ESR
	http://volkov.oce.orst.edu/tides/region.html
	https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
	ftp://ftp.esr.org/pub/datasets/tmd/

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users
	scipy: Scientific Tools for Python
		http://www.scipy.org/

PROGRAM DEPENDENCIES:
	convert_xy_ll.py: converts lat/lon points to and from projected coordinates

UPDATE HISTORY:
	Updated 07/2018: added different interpolation methods
	Updated 09/2017: Adapted for Python
"""
import os
import numpy as np
import scipy.interpolate
from convert_xy_ll import convert_xy_ll

#-- extract tidal harmonic constants out of a tidal model at coordinates
def extract_tidal_constants(ilon, ilat, grid_file, model_file, EPSG, type,
	METHOD='linear'):
	#-- read the OTIS-format tide grid file
	xi,yi,hz,mz,iob,dt = read_tide_grid(grid_file)
	#-- run wrapper function to convert coordinate systems of input lat/lon
	x,y = convert_xy_ll(ilon,ilat,EPSG,'F')
	#-- grid step size of tide model
	dx = xi[1] - xi[0]
	dy = yi[1] - yi[0]

	if (type != 'z'):
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

	#-- replace zero values with nan
	hz[hz==0] = np.nan
	if (type != 'z'):
		#-- replace original values with extend matrices
		if GLOBAL:
			hu = extend_matrix(hu)
			hv = extend_matrix(hv)
			mu = extend_matrix(mu)
			mv = extend_matrix(mv)
		#-- replace zero values with nan
		hu[hu==0] = np.nan
		hv[hv==0] = np.nan

	if (METHOD == 'bilinear'):
		D = bilinear_interp(xi,yi,hz,x,y)
		mz1 = bilinear_interp(xi,yi,mz,x,y)
	else:
		#-- use scipy interpolate to interpolate values
		D = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
			hz.flatten(), zip(x,y), method=METHOD)
		mz1 = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
			mz.real.flatten(), zip(x,y), method=METHOD)

	#-- u and v are velocities in cm/s
	if type in ('v','u'):
		unit_conv = (D*100.0)
	#-- U and V are transports in m^2/s
	elif type in ('V','U'):
		unit_conv = 1.0

	#-- read and interpolate each constituent
	constituents,nc = read_constituents(model_file)
	npts = len(D)
	amplitude = np.zeros((npts,nc))
	phase = np.zeros((npts,nc))
	for i,c in enumerate(constituents):
		if (type == 'z'):
			#-- read constituent from elevation file
			z = read_elevation_file(model_file,i)
			#-- replace original values with extend matrices
			if GLOBAL:
				z = extend_matrix(z)
			#-- replace zero values with nan
			z[z==0] = np.nan
			#-- interpolate amplitude and phase of the constituent
			if (METHOD == 'bilinear'):
				z1 = bilinear_interp(xi,yi,z,x,y)
			else:
				#-- use scipy interpolate to interpolate values
				z1 = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
					z.flatten(), zip(x,y), method=METHOD)
			#-- amplitude and phase of the constituent
			amplitude[:,i] = np.abs(z1)
			phase[:,i] = np.arctan2(-np.imag(z1),np.real(z1))
		elif type in ('U','u'):
			#-- read constituent from transport file
			u,v = read_transport_file(model_file,i)
			#-- replace original values with extend matrices
			if GLOBAL:
				u = extend_matrix(u)
			#-- replace zero values with nan
			u[u==0] = np.nan
			#-- interpolate values
			if (METHOD == 'bilinear'):
				u1 = bilinear_interp(xi,yi,u,x,y)
				mu1 = bilinear_interp(xi,yi,mu,x,y)
			else:
				#-- use scipy interpolate to interpolate values
				u1 = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
					u.flatten(), zip(x,y), method=METHOD)
				mu1 = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
					mu.real.flatten(), zip(x,y), method=METHOD)
			#-- convert units
			u1 = u1/unit_conv
			#-- amplitude and phase of the constituent
			amplitude[:,i] = np.abs(z1)
			phase[:,i] = np.arctan2(-np.imag(z1),np.real(z1))
		elif type in ('V','v'):
			#-- read constituent from transport file
			u,v = read_transport_file(input_file,i)
			#-- replace original values with extend matrices
			if GLOBAL:
				v = extend_matrix(v)
			#-- replace zero values with nan
			v[v==0] = np.nan
			#-- interpolate values
			if (METHOD == 'bilinear'):
				v1 = bilinear_interp(xi,yi,v,x,y)
				mv1 = bilinear_interp(xi,yi,mv,x,y)
			else:
				#-- use scipy interpolate to interpolate values
				v1 = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
					v.flatten(), zip(x,y), method=METHOD)
				mv1 = scipy.interpolate.griddata(zip(X.flatten(),Y.flatten()),
					mv.real.flatten(), zip(x,y), method=METHOD)
			#-- convert units
			v1 = v1/unit_conv
			#-- amplitude and phase of the constituent
			amplitude[:,i] = np.abs(v1)
			phase[:,i] = np.arctan2(-np.imag(v1),np.real(v1))
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
	temp = np.zeros((ny,nx+2),dtype=input_matrix.dtype)
	temp[:,0] = input_matrix[:,-1]
	temp[:,1:-1] = input_matrix[:,:]
	temp[:,-1] = input_matrix[:,0]
	return temp

#-- read tide grid file
def read_tide_grid(input_file):
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
	x = np.arange(xlim[0]+dx/2.0,xlim[1]+dx/2.0,dx)
	y = np.arange(ylim[0]+dy/2.0,ylim[1]+dy/2.0,dy)
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

#-- read list of constituents from an elevation or transport file
def read_constituents(input_file):
	#-- open the file
	fid = open(os.path.expanduser(input_file),'rb')
	ll, = np.fromfile(fid, dtype=np.dtype('>i4'), count=1)
	nx,ny,nc = np.fromfile(fid, dtype=np.dtype('>i4'), count=3)
	fid.seek(16,1)
	constituents = [c.rstrip() for c in fid.read(nc*4).split()]
	fid.close()
	return constituents,nc

#-- read elevation file to extract real and imaginary components for constituent
def read_elevation_file(input_file,ic):
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
	h = np.zeros((ny,nx),dtype=np.complex64)
	for i in range(ny):
		temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=2*nx)
		h.real[i,:] = temp[0:2*nx-1:2]
		h.imag[i,:] = temp[1:2*nx:2]
	#-- close the file
	fid.close()
	#-- return the elevation
	return h

#-- read transport file to extract real and imaginary components for constituent
def read_transport_file(input_file,ic):
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
	u = np.zeros((ny,nx),dtype=np.complex64)
	v = np.zeros((ny,nx),dtype=np.complex64)
	for i in range(ny):
		temp = np.fromfile(fid, dtype=np.dtype('>f4'), count=4*nx)
		u.real[i,:] = temp[0:4*nx-3:4]
		u.imag[i,:] = temp[1:4*nx-2:4]
		v.real[i,:] = temp[2:4*nx-1:4]
		v.imag[i,:] = temp[3:4*nx:4]
	#-- close the file
	fid.close()
	#-- return the transport components
	return (u,v)

#-- For a rectangular bathymetry grid:
#-- construct masks for zeta, u and v nodes on a C-grid
def Muv(hz):
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
	mu[indy,:] = mz*mz[indy,:]
	mv[:,indx] = mz*mz[:,indx]
	return (mu,mv,mz)

def Huv(hz):
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
	#-- interpolate gridded ur values to data
	data = np.zeros_like(lon)
	for i,l in enumerate(lon):
		#-- calculating the indices for the original grid
		dx = (ilon - np.floor(lon[i]/dlon)*dlon)**2
		dy = (ilat - np.floor(lat[i]/dlat)*dlat)**2
		iph, = np.nonzero(dx == np.min(dx))[0]
		ith, = np.nonzero(dy == np.min(dy))[0]
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
