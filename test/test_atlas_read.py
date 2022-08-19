#!/usr/bin/env python
u"""
test_atlas_read.py (09/2021)
Tests that ATLAS compact and netCDF4 data can be downloaded from AWS S3 bucket
Tests the read program to verify that constituents are being extracted

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
    Updated 09/2021: added test for model definition files
    Updated 03/2021: use pytest fixture to setup and teardown model data
        simplified netcdf inputs to be similar to binary OTIS read program
        replaced numpy bool/int to prevent deprecation warnings
    Written 09/2020
"""
import os
import re
import io
import gzip
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
import pyTMD.read_tide_model
import pyTMD.read_netcdf_model

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

# PURPOSE: Download TPXO8 ATLAS compact constituents from AWS S3 bucket
# @pytest.fixture(scope="module", autouse=True)
def download_TPXO8(aws_access_key_id,aws_secret_access_key,aws_region_name):
    # get aws session object
    session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name)
    # get s3 object and bucket object for pytmd data
    s3 = session.resource('s3')
    bucket = s3.Bucket('pytmd')

    # model parameters for TPXO8-ATLAS
    model_directory = os.path.join(filepath,'tpxo8_atlas')
    # recursively create model directory
    os.makedirs(model_directory)
    # retrieve each model file from s3
    for f in ['grid_tpxo8atlas_30_v1','hf.tpxo8_atlas_30_v1','uv.tpxo8_atlas_30_v1']:
        # retrieve constituent file
        obj = bucket.Object(key=posixpath.join('tpxo8_atlas',f))
        response = obj.get()
        # save constituent data
        with open(os.path.join(model_directory,f), 'wb') as destination:
            shutil.copyfileobj(response['Body'], destination)
        assert os.access(os.path.join(model_directory,f), os.F_OK)
    # run tests
    yield
    # clean up model
    shutil.rmtree(model_directory)

# PURPOSE: Download TPXO9 ATLAS V2 netCDF constituents from AWS S3 bucket
@pytest.fixture(scope="module", autouse=True)
def download_TPXO9_v2(aws_access_key_id,aws_secret_access_key,aws_region_name):
    # get aws session object
    session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name)
    # get s3 object and bucket object for pytmd data
    s3 = session.resource('s3')
    bucket = s3.Bucket('pytmd')

    # model parameters for TPXO9-atlas-v2
    model = pyTMD.model(filepath,format='netcdf',compressed=True,
        verify=False).elevation('TPXO9-atlas-v2')
    # recursively create model directory
    os.makedirs(model.model_directory)
    # retrieve grid file from s3
    f = os.path.basename(model.grid_file)
    obj = bucket.Object(key=posixpath.join('TPXO9_atlas_v2',f))
    response = obj.get()
    # save grid data
    with open(model.grid_file, 'wb') as destination:
        shutil.copyfileobj(response['Body'], destination)
    assert os.access(model.grid_file, os.F_OK)
    # retrieve each model file from s3
    for model_file in model.model_file:
        # retrieve constituent file
        f = os.path.basename(model_file)
        obj = bucket.Object(key=posixpath.join('TPXO9_atlas_v2',f))
        response = obj.get()
        # save constituent data
        with open(model_file, 'wb') as destination:
            shutil.copyfileobj(response['Body'], destination)
        assert os.access(model_file, os.F_OK)
    # run tests
    yield
    # clean up model
    shutil.rmtree(model.model_directory)

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','nearest'])
@pytest.mark.parametrize("EXTRAPOLATE", [False])
# PURPOSE: Tests that interpolated results are comparable to OTPSnc program
def test_read_TPXO9_v2(METHOD, EXTRAPOLATE):
    # model parameters for TPXO9-atlas-v2
    model_directory = os.path.join(filepath,'TPXO9_atlas_v2')
    # model grid file
    grid_file = os.path.join(model_directory,'grid_tpxo9_atlas_30_v2.nc.gz')
    # constituent files included in test
    model_files = ['h_m2_tpxo9_atlas_30_v2.nc.gz','h_s2_tpxo9_atlas_30_v2.nc.gz',
        'h_k1_tpxo9_atlas_30_v2.nc.gz','h_o1_tpxo9_atlas_30_v2.nc.gz']
    model_file = [os.path.join(model_directory,m) for m in model_files]
    constituents = ['m2','s2','k1','o1']
    TYPE = 'z'
    SCALE = 1.0/1000.0
    GZIP = True

    # read validation dataset (m2, s2, k1, o1)
    names = ('Lat', 'Lon', 'm2_amp', 'm2_ph', 's2_amp', 's2_ph',
        'k1_amp', 'k1_ph', 'o1_amp', 'o1_ph')
    formats = ('f','f','f','f','f','f','f','f','f','f')
    val = np.loadtxt(os.path.join(filepath,'extract_HC_sample_out.gz'),
        skiprows=3,dtype=dict(names=names,formats=formats))

    # extract amplitude and phase from tide model
    amp,ph,D,c = pyTMD.read_netcdf_model.extract_netcdf_constants(
        val['Lon'], val['Lat'], grid_file, model_file, type=TYPE,
        method=METHOD, extrapolate=EXTRAPOLATE, scale=SCALE,
        compressed=GZIP)
    # convert phase from 0:360 to -180:180
    ph[ph > 180] -= 360.0

    # will verify differences between model outputs are within tolerance
    amp_eps = 0.05
    ph_eps = 10.0
    # calculate differences between OTPSnc and python version
    for i,cons in enumerate(c):
        # verify constituents
        assert (cons == constituents[i])
        # calculate difference in amplitude and phase
        amp_diff = amp[:,i] - val['{0}_amp'.format(cons)]
        ph_diff = ph[:,i] - val['{0}_ph'.format(cons)]
        assert np.all(np.abs(amp_diff) <= amp_eps)
        assert np.all(np.abs(ph_diff) <= ph_eps)

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['bilinear'])
@pytest.mark.parametrize("EXTRAPOLATE", [False])
@pytest.mark.skip(reason='Need to validate over grounded point')
# PURPOSE: Tests that interpolated results are comparable to OTPS2 program
def test_verify_TPXO8(METHOD, EXTRAPOLATE):
    # model parameters for TPXO8-atlas
    model = pyTMD.model(filepath,compressed=False).elevation('TPXO8-atlas')
    # constituents for test
    constituents = ['m2','s2']

    # compile numerical expression operator
    rx = re.compile(r'[-+]?(?:(?:\d+\.\d+\.\d+)|(?:\d+\:\d+\:\d+)'
        r'|(?:\d*\.\d+)|(?:\d+\.?))')
    # read validation dataset (m2, s2)
    # Lat  Lon  mm.dd.yyyy hh:mm:ss  z(m)  Depth(m)
    with gzip.open(os.path.join(filepath,'predict_tide.out.gz'),'r') as f:
        file_contents = f.read().decode('ISO-8859-1').splitlines()
    # number of validation data points
    nval = len(file_contents) - 13
    # allocate for validation dataset
    val = dict(latitude=np.zeros((nval)),longitude=np.zeros((nval)),
        time=np.zeros((nval)),height=np.zeros((nval)))
    # counter for filling variables
    j = 0
    # for each line in the validation file
    for i,line in enumerate(file_contents):
        # extract numerical values
        line_contents = rx.findall(line)
        # skip line if not a data line
        if (len(line_contents) != 6):
            continue
        # save longitude, latitude and tide height
        val['latitude'][j] = np.float64(line_contents[0])
        val['longitude'][j] = np.float64(line_contents[1])
        val['height'][j] = np.float64(line_contents[4])
        # extract dates
        MM,DD,YY = np.array(line_contents[2].split('.'),dtype='f')
        hh,mm,ss = np.array(line_contents[3].split(':'),dtype='f')
        # convert from calendar dates into days since 1992-01-01T00:00:00
        val['time'][j] = pyTMD.time.convert_calendar_dates(YY, MM, DD,
            hour=hh, minute=mm, second=ss, epoch=(1992,1,1,0,0,0))
        # add to counter
        j += 1

    # extract amplitude and phase from tide model
    amp,ph,D,c = pyTMD.read_tide_model.extract_tidal_constants(
        val['longitude'], val['latitude'], model.grid_file,
        model.model_file, model.projection, type=model.type,
        method=METHOD, extrapolate=EXTRAPOLATE, grid=model.format)
    # delta time
    deltat = np.zeros_like(val['time'])
    # calculate complex phase in radians for Euler's
    # calculate constituent oscillations
    hc = amp*np.exp(-1j*ph*np.pi/180.0)
    # find index to reduce to list of wanted constituents
    i = [c.index(cons) for cons in constituents]

    # allocate for out tides at point
    tide = np.ma.zeros((nval))
    tide.mask = np.zeros((nval),dtype=bool)
    # predict tidal elevations at time
    tide.mask[:] = np.any(hc.mask, axis=1)
    tide.data[:] = pyTMD.predict_tide_drift(val['time'], hc[:,i],
        constituents, deltat=deltat, corrections=model.format)

    # will verify differences between model outputs are within tolerance
    eps = 0.03
    # calculate differences between OTPS2 and python version
    difference = np.ma.zeros((nval))
    difference.data[:] = tide.data - val['height']
    difference.mask = np.copy(tide.mask)
    if not np.all(difference.mask):
        assert np.all(np.abs(difference) <= eps)

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline'])
@pytest.mark.parametrize("EXTRAPOLATE", [False])
# PURPOSE: Tests that interpolated results are comparable to OTPSnc program
def test_verify_TPXO9_v2(METHOD, EXTRAPOLATE):
    # model parameters for TPXO9-atlas-v2
    model_directory = os.path.join(filepath,'TPXO9_atlas_v2')
    # model grid file
    grid_file = os.path.join(model_directory,'grid_tpxo9_atlas_30_v2.nc.gz')
    # constituent files included in test
    model_files = ['h_m2_tpxo9_atlas_30_v2.nc.gz','h_s2_tpxo9_atlas_30_v2.nc.gz',
        'h_k1_tpxo9_atlas_30_v2.nc.gz','h_o1_tpxo9_atlas_30_v2.nc.gz']
    model_file = [os.path.join(model_directory,m) for m in model_files]
    constituents = ['m2','s2','k1','o1']
    model_format = 'netcdf'
    TYPE = 'z'
    SCALE = 1.0/1000.0
    GZIP = True

    # compile numerical expression operator
    rx = re.compile(r'[-+]?(?:(?:\d+\.\d+\.\d+)|(?:\d+\:\d+\:\d+)'
        r'|(?:\d*\.\d+)|(?:\d+\.?))')
    # read validation dataset (m2, s2, k1, o1)
    # Lat  Lon  mm.dd.yyyy hh:mm:ss  z(m)  Depth(m)
    with gzip.open(os.path.join(filepath,'predict_tide_sample_out.gz'),'r') as f:
        file_contents = f.read().decode('ISO-8859-1').splitlines()
    # number of validation data points
    nval = len(file_contents) - 6
    # allocate for validation dataset
    val = dict(latitude=np.zeros((nval)),longitude=np.zeros((nval)),
        time=np.zeros((nval)),height=np.zeros((nval)))
    for i,line in enumerate(file_contents[6:]):
        # extract numerical values
        line_contents = rx.findall(line)
        val['latitude'][i] = np.float64(line_contents[0])
        val['longitude'][i] = np.float64(line_contents[1])
        val['height'][i] = np.float64(line_contents[4])
        # extract dates
        MM,DD,YY = np.array(line_contents[2].split('.'),dtype='f')
        hh,mm,ss = np.array(line_contents[3].split(':'),dtype='f')
        # convert from calendar dates into days since 1992-01-01T00:00:00
        val['time'][i] = pyTMD.time.convert_calendar_dates(YY, MM, DD,
            hour=hh, minute=mm, second=ss, epoch=(1992,1,1,0,0,0))

    # extract amplitude and phase from tide model
    amp,ph,D,c = pyTMD.read_netcdf_model.extract_netcdf_constants(
        val['longitude'], val['latitude'], grid_file, model_file,
        type=TYPE, method=METHOD, extrapolate=EXTRAPOLATE,
        scale=SCALE, compressed=GZIP)
    # delta time
    deltat = np.zeros_like(val['time'])
    # verify constituents
    assert (c == constituents)
    # calculate complex phase in radians for Euler's
    # calculate constituent oscillations
    hc = amp*np.exp(-1j*ph*np.pi/180.0)

    # allocate for out tides at point
    tide = np.ma.zeros((nval))
    tide.mask = np.zeros((nval),dtype=bool)
    # predict tidal elevations at time
    tide.mask[:] = np.any(hc.mask, axis=1)
    tide.data[:] = pyTMD.predict_tide_drift(val['time'], hc, c,
        deltat=deltat, corrections=model_format)

    # will verify differences between model outputs are within tolerance
    eps = 0.05
    # calculate differences between OTPSnc and python version
    difference = np.ma.zeros((nval))
    difference.data[:] = tide.data - val['height']
    difference.mask = np.copy(tide.mask)
    if not np.all(difference.mask):
        assert np.all(np.abs(difference) <= eps)

# parameterize ATLAS tide model
@pytest.mark.parametrize("MODEL", ['TPXO8-atlas','TPXO9-atlas-v2'])
# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','nearest'])
@pytest.mark.parametrize("EXTRAPOLATE", [False])
# PURPOSE: test the tide correction wrapper function
@pytest.mark.skip(reason='does not presently validate the ATLAS outputs')
def test_Ross_Ice_Shelf(MODEL, METHOD, EXTRAPOLATE):
    # create an image around the Ross Ice Shelf
    xlimits = np.array([-750000,550000])
    ylimits = np.array([-1450000,-300000])
    spacing = np.array([50e3,-50e3])
    # x and y coordinates
    x = np.arange(xlimits[0],xlimits[1]+spacing[0],spacing[0])
    y = np.arange(ylimits[1],ylimits[0]+spacing[1],spacing[1])
    xgrid,ygrid = np.meshgrid(x,y)
    # time dimension
    delta_time = 0.0
    # calculate tide map
    tide = pyTMD.compute_tide_corrections(xgrid, ygrid, delta_time,
        DIRECTORY=filepath, MODEL=MODEL, GZIP=True,
        EPOCH=(2000,1,1,0,0,0), TYPE='grid', TIME='TAI',
        EPSG=3031, METHOD=METHOD, EXTRAPOLATE=EXTRAPOLATE)
    assert np.any(tide)

# PURPOSE: test definition file functionality
@pytest.mark.parametrize("MODEL", ['TPXO9-atlas-v2'])
def test_definition_file(MODEL):
    # get model parameters
    model = pyTMD.model(filepath,compressed=True).elevation(MODEL)
    # create model definition file
    fid = io.StringIO()
    attrs = ['name','format','grid_file','model_file','compressed','type','scale']
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
