#!/usr/bin/env python
u"""
test_fes_predict.py (09/2021)
Tests that FES2014 data can be downloaded from AWS S3 bucket
Tests the read program to verify that constituents are being extracted
Tests that interpolated results are comparable to FES2014 program

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    boto3: Amazon Web Services (AWS) SDK for Python
        https://boto3.amazonaws.com/v1/documentation/api/latest/index.html

UPDATE HISTORY:
    Updated 09/2021: update check tide points to add compression flags
    Updated 05/2021: added test for check point program
    Updated 03/2021: use pytest fixture to setup and teardown model data
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Written 08/2020
"""
import os
import io
import boto3
import shutil
import pytest
import inspect
import warnings
import posixpath
import numpy as np
import pyTMD.time
import pyTMD.model
import pyTMD.utilities
import pyTMD.read_FES_model
import pyTMD.predict_tide_drift
import pyTMD.infer_minor_corrections
import pyTMD.check_tide_points
import pyTMD.calc_delta_time

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

# PURPOSE: Download FES2014 constituents from AWS S3 bucket
@pytest.fixture(scope="module", autouse=True)
def download_model(aws_access_key_id,aws_secret_access_key,aws_region_name):
    # get aws session object
    session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name)
    # get s3 object and bucket object for pytmd data
    s3 = session.resource('s3')
    bucket = s3.Bucket('pytmd')

    # model parameters for FES2014
    model = pyTMD.model(filepath,compressed=True,
        verify=False).elevation('FES2014')
    # recursively create model directory
    os.makedirs(model.model_directory)
    # retrieve each model file from s3
    for model_file in model.model_file:
        # retrieve constituent file
        f = os.path.basename(model_file)
        obj = bucket.Object(key=posixpath.join('fes2014','ocean_tide',f))
        response = obj.get()
        # save constituent data
        with open(model_file, 'wb') as destination:
            shutil.copyfileobj(response['Body'], destination)
        assert os.access(model_file, os.F_OK)
    # run tests
    yield
    # clean up model
    shutil.rmtree(model.model_directory)

# PURPOSE: Tests check point program
def test_check_FES2014():
    lons = np.zeros((10)) + 178.0
    lats = -45.0 - np.arange(10)*5.0
    obs = pyTMD.check_tide_points(lons, lats, DIRECTORY=filepath,
        MODEL='FES2014', GZIP=True, EPSG=4326)
    exp = np.array([True, True, True, True, True,
        True, True, True, False, False])
    assert np.all(obs == exp)

# PURPOSE: Tests that interpolated results are comparable to FES program
def test_verify_FES2014():
    # model parameters for FES2014
    model_directory = os.path.join(filepath,'fes2014','ocean_tide')
    # constituent files included in test
    model_files = ['2n2.nc.gz','k1.nc.gz','k2.nc.gz','m2.nc.gz','m4.nc.gz',
        'mf.nc.gz','mm.nc.gz','msqm.nc.gz','mtm.nc.gz','n2.nc.gz','o1.nc.gz',
        'p1.nc.gz','q1.nc.gz','s1.nc.gz','s2.nc.gz']
    model_file = [os.path.join(model_directory,m) for m in model_files]
    c = ['2n2','k1','k2','m2','m4','mf','mm','msqm','mtm','n2','o1',
        'p1','q1','s1','s2']
    model_format = 'FES'
    VERSION = 'FES2014'
    TYPE = 'z'
    SCALE = 1.0/100.0

    # read validation dataset
    # extract time (Modified Julian Days), latitude, longitude, and tide data
    names = ('CNES','Hour','Latitude','Longitude','Short_tide','LP_tide',
        'Pure_tide','Geo_tide','Rad_tide')
    formats = ('f','i','f','f','f','f','f','f','f')
    file_contents = np.loadtxt(os.path.join(filepath,'fes_slev.txt.gz'),
        skiprows=1,dtype=dict(names=names,formats=formats))
    longitude = file_contents['Longitude']
    latitude = file_contents['Latitude']
    # convert short tide estimates to meters
    validation = file_contents['Short_tide']/100.0
    npts = len(file_contents)

    # convert time from CNES Julian Days to days since 1992-01-01T00:00:00
    # CNES Julian Days = Days relative to 1950-01-01 (MJD:33282)
    tide_time = file_contents['CNES'] - 15340.0

    # extract amplitude and phase from tide model
    amp,ph = pyTMD.read_FES_model.extract_FES_constants(longitude,
        latitude, model_file, type=TYPE, version=VERSION,
        method='spline', compressed=True, scale=SCALE,)
    # interpolate delta times from calendar dates to tide time
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
    deltat = pyTMD.calc_delta_time(delta_file, tide_time)
    # calculate complex phase in radians for Euler's
    # calculate constituent oscillations
    hc = amp*np.exp(-1j*ph*np.pi/180.0)

    # allocate for out tides at point
    tide = np.ma.zeros((npts))
    tide.mask = np.zeros((npts),dtype=bool)
    # predict tidal elevations at time and infer minor corrections
    tide.mask[:] = np.any(hc.mask, axis=1)
    tide.data[:] = pyTMD.predict_tide_drift(tide_time, hc, c,
        deltat=deltat, corrections=model_format)
    minor = pyTMD.infer_minor_corrections(tide_time, hc, c,
        deltat=deltat, corrections=model_format)
    tide.data[:] += minor.data[:]

    # will verify differences between model outputs are within tolerance
    eps = 0.05
    # calculate differences between fes2014 and python version
    difference = np.ma.zeros((npts))
    difference.data[:] = tide.data - validation
    difference.mask = np.copy(tide.mask)
    if not np.all(difference.mask):
        assert np.all(np.abs(difference) <= eps)

# PURPOSE: test definition file functionality
@pytest.mark.parametrize("MODEL", ['FES2014'])
def test_definition_file(MODEL):
    # get model parameters
    model = pyTMD.model(filepath,compressed=True).elevation(MODEL)
    # create model definition file
    fid = io.StringIO()
    attrs = ['name','format','model_file','compressed','type','scale','version']
    for attr in attrs:
        val = getattr(model,attr)
        if isinstance(val,list):
            fid.write('{0}\t{1}\n'.format(attr,','.join(val)))
        else:
            fid.write('{0}\t{1}\n'.format(attr,val))
    fid.seek(0)
    # use model definition file as input
    m = pyTMD.model().from_file(fid)
    for attr in attrs:
        assert getattr(model,attr) == getattr(m,attr)
