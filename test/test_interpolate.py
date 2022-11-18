#!/usr/bin/env python
u"""
test_interpolate.py (03/2021)
Test the interpolation and extrapolation routines

UPDATE HISTORY:
    Updated 11/2022: use f-strings for formatting verbose or ascii output
	Written 03/2021
"""
import os
import inspect
import pytest
import numpy as np
import scipy.io
import scipy.interpolate
import pyTMD.spatial
import pyTMD.utilities
import pyTMD.bilinear_interp
import pyTMD.nearest_extrap

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

# PURPOSE: Download max determinant nodes from spherepts
# https://github.com/gradywright/spherepts
@pytest.fixture(scope="module", autouse=True)
def download_nodes(N=324):
	matfile = f'md{N:05d}.mat'
	HOST = ['https://github.com','gradywright','spherepts','raw',
		'master','nodes','max_determinant',matfile]
	pyTMD.utilities.from_http(HOST,
		local=os.path.join(filepath,matfile),
		verbose=True)
	yield
	# remove the node file
	os.remove(os.path.join(filepath,matfile))

# Franke's 3D evaluation function
def franke_3d(x,y,z):
	F1 = 0.75*np.exp(-((9.*x-2.)**2 + (9.*y-2.)**2 + (9.0*z-2.)**2)/4.)
	F2 = 0.75*np.exp(-((9.*x+1.)**2/49. + (9.*y+1.)/10. + (9.0*z+1.)/10.))
	F3 = 0.5*np.exp(-((9.*x-7.)**2 + (9.*y-3.)**2 + (9.*z-5)**2)/4.)
	F4 = 0.2*np.exp(-((9.*x-4.)**2 + (9.*y-7.)**2 + (9.*z-5.)**2))
	F = F1 + F2 + F3 - F4
	return F

# use max determinant nodes from spherepts
def test_cartesian(N=324):
	# read the node file
	matfile = f'md{N:05d}.mat'
	xd = scipy.io.loadmat(os.path.join(filepath,matfile))
	x,y,z = (xd['x'][:,0],xd['x'][:,1],xd['x'][:,2])
	# convert from cartesian to sphere
	lon,lat,r = pyTMD.spatial.to_sphere(x,y,z)
	X,Y,Z = pyTMD.spatial.to_cartesian(lon,lat,a_axis=r,flat=0.0)
	# verify that coordinates are within tolerance
	assert np.all(np.isclose(x,X))
	assert np.all(np.isclose(y,Y))
	assert np.all(np.isclose(z,Z))

# use max determinant nodes from spherepts
def test_geodetic(N=324):
	# read the node file
	matfile = f'md{N:05d}.mat'
	xd = scipy.io.loadmat(os.path.join(filepath,matfile))
	# convert from cartesian to sphere
	ln,lt,_ = pyTMD.spatial.to_sphere(xd['x'][:,0],
		xd['x'][:,1],xd['x'][:,2])
	# convert from sphere to cartesian
	X,Y,Z = pyTMD.spatial.to_cartesian(ln,lt)
	# convert from cartesian to geodetic
	lon = np.zeros((N))
	lat = np.zeros((N))
	h = np.zeros((N))
	for i in range(N):
		lon[i],lat[i],h[i] = pyTMD.spatial.to_geodetic(X[i],Y[i],Z[i])
	# fix coordinates to be 0:360
	lon[lon < 0] += 360.0
	# verify that coordinates are within tolerance
	assert np.all(np.isclose(ln,lon))
	assert np.all(np.isclose(lt,lat))

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','linear','bilinear'])
# PURPOSE: test interpolation routines over a sphere
def test_interpolate(METHOD, N=324):
	# read the node file
	matfile = f'md{N:05d}.mat'
	xd = scipy.io.loadmat(os.path.join(filepath,matfile))
	x,y,z = (xd['x'][:,0],xd['x'][:,1],xd['x'][:,2])
	# convert from cartesian to sphere
	lon,lat,_ = pyTMD.spatial.to_sphere(x,y,z)
	# compute functional values at nodes
	val = franke_3d(x,y,z)
	# calculate output points (standard lat/lon grid)
	dlon,dlat = (1.0,1.0)
	LON = np.arange(0,360+dlon,dlon)
	LAT = np.arange(-90,90+dlat,dlat)
	ny,nx = (len(LAT),len(LON))
	gridlon,gridlat = np.meshgrid(LON,LAT)
	X,Y,Z = pyTMD.spatial.to_cartesian(gridlon,gridlat,
		a_axis=1.0,flat=0.0)
	# calculate functional values at output points
	FI = np.ma.zeros((ny,nx))
	FI.data[:] = franke_3d(X,Y,Z)
	FI.mask = np.zeros((ny,nx),dtype=bool)
	# use interpolation routines to get values
	if (METHOD == 'bilinear'):
		# use quick bilinear to interpolate values
		test = pyTMD.bilinear_interp(LON,LAT,FI,lon,lat)
	elif (METHOD == 'spline'):
		# use scipy bivariate splines to interpolate values
		f1 = scipy.interpolate.RectBivariateSpline(LON,LAT,
			FI.T,kx=1,ky=1)
		test = f1.ev(lon,lat)
	else:
		# use scipy regular grid to interpolate values
		r1 = scipy.interpolate.RegularGridInterpolator((LAT,LON),FI,
			method=METHOD,bounds_error=False)
		test = r1.__call__(np.c_[lat,lon])
	# verify that coordinates are within tolerance
	eps = np.finfo(np.float16).eps
	assert np.all(np.isclose(val,test,atol=eps))

# PURPOSE: test extrapolation over a sphere
def test_extrapolate(N=324):
	# read the node file
	matfile = f'md{N:05d}.mat'
	xd = scipy.io.loadmat(os.path.join(filepath,matfile))
	x,y,z = (xd['x'][:,0],xd['x'][:,1],xd['x'][:,2])
	# convert from cartesian to sphere
	lon,lat,_ = pyTMD.spatial.to_sphere(x,y,z)
	# compute functional values at nodes
	val = franke_3d(x,y,z)
	# calculate output points (standard lat/lon grid)
	dlon,dlat = (1.0,1.0)
	LON = np.arange(0,360+dlon,dlon)
	LAT = np.arange(90,-90-dlat,-dlat)
	ny,nx = (len(LAT),len(LON))
	gridlon,gridlat = np.meshgrid(LON,LAT)
	X,Y,Z = pyTMD.spatial.to_cartesian(gridlon,gridlat,
		a_axis=1.0,flat=0.0)
	# calculate functional values at output points
	FI = np.ma.zeros((ny,nx))
	FI.data[:] = franke_3d(X,Y,Z)
	FI.mask = np.zeros((ny,nx),dtype=bool)
	# use nearest neighbors extrapolation to points
	test = pyTMD.nearest_extrap(LON,LAT,FI,lon,lat,EPSG='4326')
	# verify that coordinates are within tolerance
	assert np.all(np.isclose(val,test,atol=0.1))

# PURPOSE: test that extrapolation will not occur if invalid
def test_extrapolation_checks(N=324):
	# read the node file
	matfile = f'md{N:05d}.mat'
	xd = scipy.io.loadmat(os.path.join(filepath,matfile))
	x,y,z = (xd['x'][:,0],xd['x'][:,1],xd['x'][:,2])
	# convert from cartesian to sphere
	lon,lat,_ = pyTMD.spatial.to_sphere(x,y,z)
	# calculate output points (standard lat/lon grid)
	dlon,dlat = (1.0,1.0)
	LON = np.arange(0,360+dlon,dlon)
	LAT = np.arange(90,-90-dlat,-dlat)
	ny,nx = (len(LAT),len(LON))
	# calculate functional values at output points
	FI = np.ma.zeros((ny,nx))
	FI.mask = np.ones((ny,nx),dtype=bool)
	# use nearest neighbors extrapolation to points
	# in case where there are no valid grid points
	test = pyTMD.nearest_extrap(LON,LAT,FI,lon,lat,EPSG='4326')
	assert(np.all(test.mask))
	# use nearest neighbors extrapolation
	# in case where there are no points to be extrapolated
	test = pyTMD.nearest_extrap(LON,LAT,FI,[],[])
	assert np.logical_not(test)
