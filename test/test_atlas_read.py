#!/usr/bin/env python
u"""
test_atlas_read.py (08/2020)
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
    Written 09/2020
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
import pyTMD.read_tide_model
import pyTMD.read_netcdf_model

#-- current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

#-- PURPOSE: Download TPXO8 ATLAS compact constituents from AWS S3 bucket
def test_download_TPXO8(aws_access_key_id,aws_secret_access_key,aws_region_name):
    #-- get aws session object
    session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name)
    #-- get s3 object and bucket object for pytmd data
    s3 = session.resource('s3')
    bucket = s3.Bucket('pytmd')

    #-- model parameters for TPXO8-ATLAS
    model_directory = os.path.join(filepath,'tpxo8_atlas')
    #-- recursively create model directory
    os.makedirs(model_directory)
    #-- retrieve each model file from s3
    for f in ['grid_tpxo8atlas_30_v1','hf.tpxo8_atlas_30_v1','uv.tpxo8_atlas_30_v1']:
        #-- retrieve constituent file
        obj = bucket.Object(key=posixpath.join('tpxo8_atlas',f))
        response = obj.get()
        #-- save constituent data
        with open(os.path.join(model_directory,f), 'wb') as destination:
            shutil.copyfileobj(response['Body'], destination)
        assert os.access(os.path.join(model_directory,f), os.F_OK)

#-- PURPOSE: Download TPXO9 ATLAS V2 netCDF constituents from AWS S3 bucket
def test_download_TPXO9_v2(aws_access_key_id,aws_secret_access_key,aws_region_name):
    #-- get aws session object
    session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name)
    #-- get s3 object and bucket object for pytmd data
    s3 = session.resource('s3')
    bucket = s3.Bucket('pytmd')

    #-- model parameters for TPXO9-atlas-v2
    model_directory = os.path.join(filepath,'TPXO9_atlas_v2')
    model_files = ['grid_tpxo9_atlas_30_v2.nc.gz','h_2n2_tpxo9_atlas_30_v2.nc.gz',
        'h_k1_tpxo9_atlas_30_v2.nc.gz','h_k2_tpxo9_atlas_30_v2.nc.gz',
        'h_m2_tpxo9_atlas_30_v2.nc.gz','h_m4_tpxo9_atlas_30_v2.nc.gz',
        'h_mn4_tpxo9_atlas_30_v2.nc.gz','h_ms4_tpxo9_atlas_30_v2.nc.gz',
        'h_n2_tpxo9_atlas_30_v2.nc.gz','h_o1_tpxo9_atlas_30_v2.nc.gz',
        'h_p1_tpxo9_atlas_30_v2.nc.gz','h_q1_tpxo9_atlas_30_v2.nc.gz',
        'h_s2_tpxo9_atlas_30_v2.nc.gz']
    #-- recursively create model directory
    os.makedirs(model_directory)
    #-- retrieve each model file from s3
    for f in model_files:
        #-- retrieve constituent file
        obj = bucket.Object(key=posixpath.join('TPXO9_atlas_v2',f))
        response = obj.get()
        #-- save constituent data
        with open(os.path.join(model_directory,f), 'wb') as destination:
            shutil.copyfileobj(response['Body'], destination)
        assert os.access(os.path.join(model_directory,f), os.F_OK)

#-- parameterize ATLAS tide model
@pytest.mark.parametrize("MODEL", ['TPXO8-atlas','TPXO9-atlas-v2'])
#-- parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','nearest','bilinear'])
@pytest.mark.parametrize("EXTRAPOLATE", [True])
#-- PURPOSE: test the tide correction wrapper function
def test_Ross_Ice_Shelf(MODEL, METHOD, EXTRAPOLATE):
    #-- create an image around the Ross Ice Shelf
    xlimits = np.array([-740000,520000])
    ylimits = np.array([-1430000,-300000])
    spacing = np.array([10e3,-10e3])
    #-- x and y coordinates
    x = np.arange(xlimits[0],xlimits[1]+spacing[0],spacing[0])
    y = np.arange(ylimits[1],ylimits[0]+spacing[1],spacing[1])
    xgrid,ygrid = np.meshgrid(x,y)
    #-- time dimension
    delta_time = 0.0
    #-- calculate tide map
    tide = pyTMD.compute_tide_corrections(xgrid, ygrid, delta_time,
        DIRECTORY=filepath, MODEL=MODEL, EPOCH=(2000,1,1,0,0,0),
        TYPE='grid', TIME='TAI', EPSG=3031, METHOD=METHOD,
        EXTRAPOLATE=EXTRAPOLATE)
    assert np.any(tide)
