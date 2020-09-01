#!/usr/bin/env python
u"""
test_perth3_read.py (08/2020)
Tests that GOT4.7 data can be downloaded from AWS S3 bucket
Tests the read program to verify that constituents are being extracted
Tests that interpolated results are comparable to NASA PERTH3 program

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    boto3: Amazon Web Services (AWS) SDK for Python
        https://boto3.amazonaws.com/v1/documentation/api/latest/index.html

UPDATE HISTORY:
    Written 08/2020
"""
import os
import gzip
import boto3
import shutil
import pytest
import inspect
import warnings
import posixpath
import numpy as np
import pyTMD.time
import pyTMD.utilities
import pyTMD.read_GOT_model
import pyTMD.predict_tide_drift
import pyTMD.infer_minor_corrections
import pyTMD.compute_tide_corrections
import pyTMD.calc_delta_time

#-- current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

#-- PURPOSE: Download GOT4.7 constituents from AWS S3 bucket
def test_download_GOT47(aws_access_key_id,aws_secret_access_key,aws_region_name):
    #-- get aws session object
    session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name)
    #-- get s3 object and bucket object for pytmd data
    s3 = session.resource('s3')
    bucket = s3.Bucket('pytmd')

    #-- model parameters for GOT4.7
    model_directory = os.path.join(filepath,'GOT4.7','grids_oceantide')
    model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
        'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
    #-- recursively create model directory
    os.makedirs(model_directory)
    #-- retrieve each model file from s3
    for f in model_files:
        #-- retrieve constituent file
        obj = bucket.Object(key=posixpath.join('GOT4.7','grids_oceantide',f))
        response = obj.get()
        #-- save constituent data
        with open(os.path.join(model_directory,f), 'wb') as destination:
            shutil.copyfileobj(response['Body'], destination)
        assert os.access(os.path.join(model_directory,f), os.F_OK)

#-- parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','linear','bilinear'])
#-- PURPOSE: Tests that interpolated results are comparable to PERTH3 program
def test_verify_GOT47(METHOD):
    #-- model parameters for GOT4.7
    model_directory = os.path.join(filepath,'GOT4.7','grids_oceantide')
    #-- perth3 test program infers m4 tidal constituent
    model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
        'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz']
    c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1']
    model_format = 'GOT'
    SCALE = 1.0

    #-- read validation dataset
    with gzip.open(os.path.join(filepath,'perth_output_got4.7.gz'),'r') as fid:
        file_contents = fid.read().decode('ISO-8859-1').splitlines()
    #-- extract latitude, longitude, time (Modified Julian Days) and tide data
    npts = len(file_contents) - 2
    latitude = np.zeros((npts))
    longitude = np.zeros((npts))
    MJD = np.zeros((npts))
    validation = np.ma.zeros((npts))
    validation.mask = np.ones((npts),dtype=np.bool)
    for i in range(npts):
        line_contents = file_contents[i+2].split()
        latitude[i] = np.float(line_contents[0])
        longitude[i] = np.float(line_contents[1])
        MJD[i] = np.float(line_contents[2])
        if (len(line_contents) == 5):
            validation.data[i] = np.float(line_contents[3])
            validation.mask[i] = False

    #-- convert time from MJD to days since 1992-01-01T00:00:00
    tide_time = MJD - 48622.0

    #-- extract amplitude and phase from tide model
    amp,ph = pyTMD.read_GOT_model.extract_GOT_constants(longitude, latitude,
        model_directory, model_files, METHOD=METHOD, SCALE=SCALE)
    #-- interpolate delta times from calendar dates to tide time
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
    deltat = pyTMD.calc_delta_time(delta_file, tide_time)
    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- calculate constituent oscillations
    hc = amp*np.exp(cph)

    #-- allocate for out tides at point
    tide = np.ma.zeros((npts))
    tide.mask = np.zeros((npts),dtype=np.bool)
    #-- predict tidal elevations at time and infer minor corrections
    tide.mask[:] = np.any(hc.mask, axis=1)
    tide.data[:] = pyTMD.predict_tide_drift(tide_time, hc, c,
        DELTAT=deltat, CORRECTIONS=model_format)
    minor = pyTMD.infer_minor_corrections(tide_time, hc, c,
        DELTAT=deltat, CORRECTIONS=model_format)
    tide.data[:] += minor.data[:]

    #-- will verify differences between model outputs are within tolerance
    eps = 0.01
    #-- calculate differences between perth3 and python version
    difference = np.ma.zeros((npts))
    difference.data[:] = tide.data - validation
    difference.mask = (tide.mask | validation.mask)
    if not np.all(difference.mask):
        assert np.all(np.abs(difference) <= eps)

#-- parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','nearest','bilinear'])
#-- PURPOSE: test the tide correction wrapper function
def test_Ross_Ice_Shelf(METHOD):
    #-- create an image around the Ross Ice Shelf
    xlimits = np.array([-740000,520000])
    ylimits = np.array([-1430000,-300000])
    spacing = np.array([5e3,-5e3])
    #-- x and y coordinates
    x = np.arange(xlimits[0],xlimits[1]+spacing[0],spacing[0])
    y = np.arange(ylimits[1],ylimits[0]+spacing[1],spacing[1])
    xgrid,ygrid = np.meshgrid(x,y)
    #-- x and y dimensions
    nx = np.int((xlimits[1]-xlimits[0])/spacing[0])+1
    ny = np.int((ylimits[0]-ylimits[1])/spacing[1])+1
    #-- time dimension
    delta_time = 0.0
    #-- calculate tide map
    tide = pyTMD.compute_tide_corrections(xgrid, ygrid, delta_time,
        DIRECTORY=filepath, MODEL='GOT4.7', EPOCH=(2018,1,1,0,0,0),
        TYPE='grid', TIME='GPS', EPSG=3031, METHOD=METHOD)
    assert np.any(tide)
