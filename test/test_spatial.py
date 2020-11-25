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
    #-- read test ascii file
    input_file = os.path.join(filepath,'TEST.csv')
    test = pyTMD.spatial.from_ascii(input_file, header='YAML',
        columns=['time','y','x','data'], verbose=True)
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
    #-- read test netCDF4 file
    input_file = os.path.join(filepath,'TEST.nc')
    test = pyTMD.spatial.from_netCDF4(input_file, timename='time',
        xname='x', yname='y', varname='data', verbose=True)
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
    #-- read test HDF5 file
    input_file = os.path.join(filepath,'TEST.H5')
    test = pyTMD.spatial.from_HDF5(input_file, timename='time',
        xname='x', yname='y', varname='data', verbose=True)
    #-- check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    #-- remove the test file
    os.remove(output_file)

#-- PURPOSE: test the read and write of geotiff files
def test_geotiff(username, password):
    #-- build urllib2 opener with credentials
    pyTMD.utilities.build_opener(username, password, context=ssl.SSLContext(),
        password_manager=True, get_ca_certs=False, redirect=False,
        authorization_header=True, urs='https://urs.earthdata.nasa.gov')
    #-- download NASA Operation IceBridge DMS L3 Photogrammetric DEM
    HOST = ['https://n5eil01u.ecs.nsidc.org','ICEBRIDGE','IODEM3.001',
        '2009.10.25','IODEM3_20091025_212618_02720_DEM.tif']
    input_file = os.path.join(filepath,HOST[-1])
    pyTMD.utilities.from_http(HOST, local=input_file, context=None,
        verbose=True, mode=0o775)
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
    #-- remove the test files
    os.remove(input_file)
    os.remove(output_file)
