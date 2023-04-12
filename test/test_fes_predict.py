#!/usr/bin/env python
u"""
test_fes_predict.py (04/2023)
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
    Updated 04/2023: using pathlib to define and expand paths
    Updated 12/2022: add check for read and interpolate constants
    Updated 09/2021: update check tide points to add compression flags
    Updated 05/2021: added test for check point program
    Updated 03/2021: use pytest fixture to setup and teardown model data
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Written 08/2020
"""
import io
import boto3
import shutil
import pytest
import inspect
import pathlib
import posixpath
import numpy as np
import pyTMD.io
import pyTMD.time
import pyTMD.io.model
import pyTMD.utilities
import pyTMD.predict
import pyTMD.check_tide_points

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

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
    model = pyTMD.io.model(filepath,compressed=True,
        verify=False).elevation('FES2014')
    # recursively create model directory
    model.model_directory.mkdir(parents=True, exist_ok=True)
    # retrieve each model file from s3
    for model_file in model.model_file:
        # retrieve constituent file
        f = model_file.name
        obj = bucket.Object(key=posixpath.join('fes2014','ocean_tide',f))
        response = obj.get()
        # save constituent data
        with model_file.open(mode='wb') as destination:
            shutil.copyfileobj(response['Body'], destination)
        assert model_file.exists()
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
    model_directory = filepath.joinpath('fes2014','ocean_tide')
    # constituent files included in test
    model_files = ['2n2.nc.gz','k1.nc.gz','k2.nc.gz','m2.nc.gz','m4.nc.gz',
        'mf.nc.gz','mm.nc.gz','msqm.nc.gz','mtm.nc.gz','n2.nc.gz','o1.nc.gz',
        'p1.nc.gz','q1.nc.gz','s1.nc.gz','s2.nc.gz']
    model_file = [model_directory.joinpath(m) for m in model_files]
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
    file_contents = np.loadtxt(filepath.joinpath('fes_slev.txt.gz'),
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
    amp,ph = pyTMD.io.FES.extract_constants(longitude,
        latitude, model_file, type=TYPE, version=VERSION,
        method='spline', compressed=True, scale=SCALE)
    # interpolate delta times from calendar dates to tide time
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
    deltat = pyTMD.time.interpolate_delta_time(delta_file, tide_time)
    # calculate complex phase in radians for Euler's
    # calculate constituent oscillations
    hc = amp*np.exp(-1j*ph*np.pi/180.0)

    # allocate for out tides at point
    tide = np.ma.zeros((npts))
    tide.mask = np.zeros((npts),dtype=bool)
    # predict tidal elevations at time and infer minor corrections
    tide.mask[:] = np.any(hc.mask, axis=1)
    tide.data[:] = pyTMD.predict.drift(tide_time, hc, c,
        deltat=deltat, corrections=model_format)
    minor = pyTMD.predict.infer_minor(tide_time, hc, c,
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

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline'])
# PURPOSE: Tests that interpolated results are comparable
def test_compare_FES2014(METHOD):
    # model parameters for FES2014
    model_directory = filepath.joinpath('fes2014','ocean_tide')
    # constituent files included in test
    model_files = ['2n2.nc.gz','k1.nc.gz','k2.nc.gz','m2.nc.gz','m4.nc.gz',
        'mf.nc.gz','mm.nc.gz','msqm.nc.gz','mtm.nc.gz','n2.nc.gz','o1.nc.gz',
        'p1.nc.gz','q1.nc.gz','s1.nc.gz','s2.nc.gz']
    model_file = [model_directory.joinpath(m) for m in model_files]
    c = ['2n2','k1','k2','m2','m4','mf','mm','msqm','mtm','n2','o1',
        'p1','q1','s1','s2']
    VERSION = 'FES2014'
    TYPE = 'z'
    SCALE = 1.0

    # read validation dataset
    # extract time (Modified Julian Days), latitude, longitude, and tide data
    names = ('CNES','Hour','Latitude','Longitude','Short_tide','LP_tide',
        'Pure_tide','Geo_tide','Rad_tide')
    formats = ('f','i','f','f','f','f','f','f','f')
    file_contents = np.loadtxt(filepath.joinpath('fes_slev.txt.gz'),
        skiprows=1,dtype=dict(names=names,formats=formats))
    longitude = file_contents['Longitude']
    latitude = file_contents['Latitude']

    # extract amplitude and phase from tide model
    amp1, ph1 = pyTMD.io.FES.extract_constants(longitude,
        latitude, model_file, type=TYPE, version=VERSION,
        method=METHOD, compressed=True, scale=SCALE)
    # calculate complex form of constituent oscillation
    hc1 = amp1*np.exp(-1j*ph1*np.pi/180.0)

    # read and interpolate constituents from tide model
    constituents = pyTMD.io.FES.read_constants(model_file,
        type=TYPE, version=VERSION, compressed=True)
    amp2, ph2 = pyTMD.io.FES.interpolate_constants(longitude, latitude,
        constituents, method=METHOD, scale=SCALE)
    # calculate complex form of constituent oscillation
    hc2 = amp2*np.exp(-1j*ph2*np.pi/180.0)

    # will verify differences between model outputs are within tolerance
    eps = np.finfo(np.float16).eps
    # calculate differences between methods
    for i, cons in enumerate(c):
        # calculate difference in amplitude and phase
        difference = hc1[:,i] - hc2[:,i]
        assert np.all(np.abs(difference) <= eps)

    # validate iteration within constituents class
    for field, hc in constituents:
        # verify constituents
        assert np.ma.isMaskedArray(hc)
        # validate amplitude and phase functions
        amp = constituents.amplitude(field)
        phase = constituents.phase(field)
        assert np.ma.isMaskedArray(amp)
        assert np.ma.isMaskedArray(phase)
        # calculate complex form of constituent oscillation
        hcomplex = amp*np.exp(-1j*phase*np.pi/180.0)
        # calculate difference in amplitude and phase
        difference = hc - hcomplex
        assert np.all(np.abs(difference) <= eps)

# PURPOSE: test definition file functionality
@pytest.mark.parametrize("MODEL", ['FES2014'])
def test_definition_file(MODEL):
    # get model parameters
    model = pyTMD.io.model(filepath,compressed=True).elevation(MODEL)
    # create model definition file
    fid = io.StringIO()
    attrs = ['name','format','model_file','compressed','type','scale','version']
    for attr in attrs:
        val = getattr(model,attr)
        if isinstance(val,list):
            var = ','.join(str(v) for v in val)
            fid.write(f'{attr}\t{var}\n')
        else:
            fid.write(f'{attr}\t{val}\n')
    fid.seek(0)
    # use model definition file as input
    m = pyTMD.io.model().from_file(fid)
    for attr in attrs:
        assert getattr(model,attr) == getattr(m,attr)

# PURPOSE: test extend function
def test_extend_array():
    dlon = 1
    lon = np.arange(0, 360, dlon)
    valid = np.arange(-dlon, 360 + dlon, dlon)
    test = pyTMD.io.FES.extend_array(lon, dlon)
    assert np.all(test == valid)
