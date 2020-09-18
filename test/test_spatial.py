#!/usr/bin/env python
u"""
test_spatial.py (09/2020)
Verify file read and write with spatial utilities
"""
import os
import pytest
import warnings
import inspect
import numpy as np
import pyTMD.spatial

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
