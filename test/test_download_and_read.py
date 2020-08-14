#!/usr/bin/env python
u"""
test_download_and_read.py (08/2020)
Tests that CATS2008 data can be downloaded from the US Antarctic Program (USAP)
Tests the read program to verify that constituents are being extracted
Tests that interpolated results are comparable to Matlab TMD program
    https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5
"""
import os
import pytest
import inspect
import zipfile
import warnings
import posixpath
import numpy as np
import pyTMD.time
import pyTMD.utilities
import pyTMD.read_tide_model
import pyTMD.predict_tidal_ts
import pyTMD.infer_minor_corrections

#-- current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

#-- PURPOSE: Download CATS2008 from US Antarctic Program
def test_download_CATS2008():
    #-- download CATS2008 zip file and read as virtual file object
    HOST = ['https://www.usap-dc.org','dataset','usap-dc','601235',
        '2019-12-19T23:26:43.6Z','CATS2008.zip?dataset_id=601235']
    FILE = pyTMD.utilities.from_http(HOST)
    zfile = zipfile.ZipFile(FILE)
    print('{0} -->\n'.format(posixpath.join(*HOST)))
    #-- model files
    files = ['grid_CATS2008','hf.CATS2008.out','uv.CATS2008.out']
    assert all([f2 in [f1.filename for f1 in zfile.filelist] for f2 in files])
    #-- extract each member
    for member in zfile.filelist:
        local_file = os.path.join(filepath,member.filename)
        print('\t{0}\n'.format(local_file))
        zfile.extract(member, path=filepath)
    #-- close the zipfile object
    zfile.close()

#-- PURPOSE: Test read program that grids and constituents are as expected
def test_read_CATS2008(ny=2026,nx=1663):
    #-- model parameters for CATS2008
    grid_file = os.path.join(filepath,'grid_CATS2008')
    elevation_file = os.path.join(filepath,'hf.CATS2008.out')
    transport_file = os.path.join(filepath,'uv.CATS2008.out')
    #-- read CATS2008 grid file
    xi,yi,hz,mz,iob,dt = pyTMD.read_tide_model.read_tide_grid(grid_file)
    #-- check dimensions of input grids
    assert (hz.shape == (ny,nx))
    assert (mz.shape == (ny,nx))
    #-- check constituent list
    constituents,nc = pyTMD.read_tide_model.read_constituents(elevation_file)
    cons = ['m2','s2','n2','k2','k1','o1','p1','q1','mf','mm']
    assert all(c in constituents for c in cons)
    #-- check dimensions of input grids from elevation and transport files
    for i,c in enumerate(constituents):
        z = pyTMD.read_tide_model.read_elevation_file(elevation_file,i)
        u,v = pyTMD.read_tide_model.read_transport_file(transport_file,i)
        assert (z.shape == (ny,nx))
        assert (u.shape == (ny,nx))
        assert (v.shape == (ny,nx))

#-- PURPOSE: Tests that interpolated results are comparable to Matlab program
def test_verify_CATS2008():
    #-- model parameters for CATS2008
    grid_file = os.path.join(filepath,'grid_CATS2008')
    model_file = os.path.join(filepath,'hf.CATS2008.out')
    GRID = 'OTIS'
    EPSG = 'CATS2008'
    TYPE = 'z'

    #-- compare daily outputs at a single point
    ilon = [5.0]
    ilat = [-62.5]
    #-- calculate daily results for an entire year
    delta_time = np.arange(365)*86400.0
    #-- convert time to days since 1992-01-01T00:00:00
    tide_time = pyTMD.time.convert_delta_time(delta_time,
        epoch1=(1998,1,1,0,0,0), epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
    #-- presently not converting times to dynamic times
    deltat = np.zeros_like(tide_time)
    #-- extract amplitude and phase from tide model
    amp,ph,D,c = pyTMD.read_tide_model.extract_tidal_constants(ilon, ilat,
        grid_file, model_file, EPSG, TYPE=TYPE, METHOD='spline', GRID=GRID)

    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- calculate constituent oscillation
    hc = amp*np.exp(cph)

    #-- allocate for out tides at point
    tide = np.ma.zeros((365))
    tide.mask = np.zeros((365),dtype=np.bool)
    #-- predict tidal elevations at time and infer minor corrections
    tide.mask[:] = np.any(hc.mask)
    tide.data[:] = pyTMD.predict_tidal_ts(tide_time, hc, c,
        DELTAT=deltat, CORRECTIONS=GRID)
    minor = pyTMD.infer_minor_corrections(tide_time, hc, c,
        DELTAT=deltat, CORRECTIONS=GRID)
    tide.data[:] += minor.data[:]

    #-- read validation data from Matlab TMD program
    #-- https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5
    validation = np.loadtxt(os.path.join(filepath,'tmd_cats2008_output.txt.gz'))
    difference = tide.data - validation[:,3]
    #-- verify differences are within float32 tolerance
    eps = np.finfo(np.float32).eps
    assert all(np.abs(difference) < eps)
