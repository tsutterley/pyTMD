#!/usr/bin/env python
u"""
test_spatial.py (11/2020)
Verify file read and write with spatial utilities
"""
import os
import ssl
import pytest
import warnings
import inspect
import numpy as np
import pyTMD.spatial
import pyTMD.utilities

#-- current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

#-- PURPOSE: test the read and write of ascii files
def test_ascii():
    #-- number of data points
    n_time = 30
    #-- create a test dataset
    output = {}
    output['y'] = np.random.randint(-90,90,size=n_time)
    output['x'] = np.random.randint(-180,180,size=n_time)
    output['data'] = np.random.randn(n_time)
    output['time'] = np.random.randint(0,31557600,size=n_time)
    #-- output netCDF4 and HDF5 file attributes
    #-- will be added to YAML header in csv files
    attrib = {}
    #-- latitude
    attrib['y'] = {}
    attrib['y']['long_name'] = 'Latitude'
    attrib['y']['units'] = 'Degrees_North'
    #-- longitude
    attrib['x'] = {}
    attrib['x']['long_name'] = 'Longitude'
    attrib['x']['units'] = 'Degrees_East'
    #-- long-period equilibrium tides
    attrib['data'] = {}
    attrib['data']['long_name'] = 'Height_above_WGS84_ellipsoid'
    attrib['data']['units'] = 'meters'
    #-- time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'seconds since 2018-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    #-- create test ascii file
    output_file = os.path.join(filepath,'test.csv')
    pyTMD.spatial.to_ascii(output, attrib, output_file, delimiter=',',
        columns=['time','y','x','data'], header=True, verbose=True)
    #-- read test ascii file (change case to test find function)
    input_file = os.path.join(filepath,'TEST.csv')
    test = pyTMD.spatial.from_ascii(input_file, header='YAML',
        columns=['time','y','x','data'], verbose=True)
    #-- check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    #-- read test ascii file as bytes
    fid = open(output_file,'r')
    test = pyTMD.spatial.from_ascii(fid, compression='bytes', header='YAML',
        columns=['time','y','x','data'])
    #-- check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    #-- remove the test file
    os.remove(output_file)

#-- PURPOSE: test the read and write of netCDF4 files
def test_netCDF4():
    #-- number of data points
    n_time = 3000
    #-- create a test dataset
    output = {}
    output['y'] = np.random.randint(-90,90,size=n_time)
    output['x'] = np.random.randint(-180,180,size=n_time)
    output['data'] = np.random.randn(n_time)
    output['time'] = np.random.randint(0,31557600,size=n_time)
    #-- output netCDF4 and HDF5 file attributes
    #-- will be added to YAML header in csv files
    attrib = {}
    #-- latitude
    attrib['y'] = {}
    attrib['y']['long_name'] = 'Latitude'
    attrib['y']['units'] = 'Degrees_North'
    #-- longitude
    attrib['x'] = {}
    attrib['x']['long_name'] = 'Longitude'
    attrib['x']['units'] = 'Degrees_East'
    #-- long-period equilibrium tides
    attrib['data'] = {}
    attrib['data']['long_name'] = 'Height_above_WGS84_ellipsoid'
    attrib['data']['units'] = 'meters'
    #-- time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'seconds since 2018-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    #-- create test netCDF4 file
    output_file = os.path.join(filepath,'test.nc')
    pyTMD.spatial.to_netCDF4(output, attrib, output_file, verbose=True)
    #-- read test netCDF4 file (change case to find test function)
    input_file = os.path.join(filepath,'TEST.nc')
    test = pyTMD.spatial.from_netCDF4(input_file, timename='time',
        xname='x', yname='y', varname='data', verbose=True)
    #-- check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    #-- read test netCDF4 file as bytes
    fid = open(output_file, 'rb')
    test = pyTMD.spatial.from_netCDF4(fid, compression='bytes',
        timename='time', xname='x', yname='y', varname='data')
    #-- check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    #-- remove the test file
    os.remove(output_file)

#-- PURPOSE: test the read and write of HDF5 files
def test_HDF5():
    #-- number of data points
    n_time = 3000
    #-- create a test dataset
    output = {}
    output['y'] = np.random.randint(-90,90,size=n_time)
    output['x'] = np.random.randint(-180,180,size=n_time)
    output['data'] = np.random.randn(n_time)
    output['time'] = np.random.randint(0,31557600,size=n_time)
    #-- output netCDF4 and HDF5 file attributes
    #-- will be added to YAML header in csv files
    attrib = {}
    #-- latitude
    attrib['y'] = {}
    attrib['y']['long_name'] = 'Latitude'
    attrib['y']['units'] = 'Degrees_North'
    #-- longitude
    attrib['x'] = {}
    attrib['x']['long_name'] = 'Longitude'
    attrib['x']['units'] = 'Degrees_East'
    #-- long-period equilibrium tides
    attrib['data'] = {}
    attrib['data']['long_name'] = 'Height_above_WGS84_ellipsoid'
    attrib['data']['units'] = 'meters'
    #-- time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'seconds since 2018-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    #-- create test HDF5 file
    output_file = os.path.join(filepath,'test.H5')
    pyTMD.spatial.to_HDF5(output, attrib, output_file, verbose=True)
    #-- read test HDF5 file (change case to test find function)
    input_file = os.path.join(filepath,'TEST.H5')
    test = pyTMD.spatial.from_HDF5(input_file, timename='time',
        xname='x', yname='y', varname='data', verbose=True)
    #-- check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    #-- read test HDF5 file as bytes
    fid = open(output_file, 'rb')
    test = pyTMD.spatial.from_HDF5(fid, compression='bytes',
        timename='time', xname='x', yname='y', varname='data')
    #-- check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    #-- remove the test file
    os.remove(output_file)

#-- PURPOSE: test the read and write of geotiff files
def test_geotiff(username, password):
    #-- build urllib2 opener for NSIDC with NASA Earthdata credentials
    pyTMD.utilities.build_opener(username, password, context=ssl.SSLContext(),
        password_manager=True, get_ca_certs=False, redirect=False,
        authorization_header=True, urs='https://urs.earthdata.nasa.gov')
    #-- download NASA Operation IceBridge DMS L3 Photogrammetric DEM
    HOST = ['https://n5eil01u.ecs.nsidc.org','ICEBRIDGE','IODEM3.001',
        '2009.10.25','IODEM3_20091025_212618_02720_DEM.tif']
    input_file = os.path.join(filepath,HOST[-1])
    remote_buffer = pyTMD.utilities.from_http(HOST, local=input_file,
        context=None, verbose=True, mode=0o775)
    dinput = pyTMD.spatial.from_geotiff(input_file, verbose=True)
    #-- copy global geotiff attributes for projection and grid parameters
    attrib = {a:dinput['attributes'][a] for a in ['wkt','spacing','extent']}
    #-- copy variable attributes for data
    attrib['data'] = {}
    for key,val in dinput['attributes']['data'].items():
        if isinstance(val,np.float32):
            attrib['data'][key] = np.float(val)
        else:
            attrib['data'][key] = np.copy(val)
    #-- create test geotiff file
    output_file = os.path.join(filepath,'test.tif')
    output = {'data':dinput['data'].astype(np.float)}
    pyTMD.spatial.to_geotiff(output, attrib, output_file, verbose=True)
    #-- check that data is valid
    test = pyTMD.spatial.from_geotiff(output_file, verbose=True)
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in dinput.items())
    #-- check that data is valid from in-memory object
    test = pyTMD.spatial.from_geotiff(remote_buffer, compression='bytes')
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in dinput.items())
    #-- remove the test files
    os.remove(input_file)
    os.remove(output_file)

#-- PURPOSE: test ellipsoid conversion
def test_convert_ellipsoid():
    #-- semimajor axis (a) and flattening (f) for TP and WGS84 ellipsoids
    atop,ftop = (6378136.3,1.0/298.257)
    awgs,fwgs = (6378137.0,1.0/298.257223563)
    #-- create latitude and height array in WGS84
    lat_WGS84 = 90.0 - np.arange(181,dtype=np.float)
    elev_WGS84 = 3000.0 + np.zeros((181),dtype=np.float)
    #-- convert from WGS84 to Topex/Poseidon Ellipsoids
    lat_TPX,elev_TPX = pyTMD.spatial.convert_ellipsoid(lat_WGS84, elev_WGS84,
        awgs, fwgs, atop, ftop, eps=1e-12, itmax=10)
    #-- check minimum and maximum are expected from NSIDC IDL program
    #-- ftp://ftp.nsidc.org/DATASETS/icesat/tools/idl/ellipsoid/test_ce.pro
    minlat = np.min(lat_TPX-lat_WGS84)
    maxlat = np.max(lat_TPX-lat_WGS84)
    minelev = 100.0*np.min(elev_TPX-elev_WGS84)
    maxelev = 100.0*np.max(elev_TPX-elev_WGS84)
    assert np.isclose([minlat,maxlat],[-1.2305653e-7,1.2305653e-7]).all()
    assert np.isclose([minelev,maxelev],[70.000000,71.368200]).all()
    #-- convert back from Topex/Poseidon to WGS84 Ellipsoids
    phi_WGS84,h_WGS84 = pyTMD.spatial.convert_ellipsoid(lat_TPX, elev_TPX,
        atop, ftop, awgs, fwgs, eps=1e-12, itmax=10)
    #-- check minimum and maximum are expected from NSIDC IDL program
    #-- ftp://ftp.nsidc.org/DATASETS/icesat/tools/idl/ellipsoid/test_ce.pro
    minlatdel = np.min(phi_WGS84-lat_WGS84)
    maxlatdel = np.max(phi_WGS84-lat_WGS84)
    minelevdel = 100.0*np.min(h_WGS84-elev_WGS84)
    maxelevdel = 100.0*np.max(h_WGS84-elev_WGS84)
    assert np.isclose([minlatdel,maxlatdel],[-2.1316282e-14,2.1316282e-14]).all()
    assert np.isclose([minelevdel,maxelevdel],[-1.3287718e-7,1.6830199e-7]).all()
