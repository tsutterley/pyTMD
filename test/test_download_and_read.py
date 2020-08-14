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

#-- PURPOSE: Download Antarctic Tide Gauge Database from US Antarctic Program
def test_download_AntTG():
    #-- download Tide Gauge Database text file
    HOST = ['https://www.usap-dc.org','dataset','usap-dc','601358',
        '2020-07-10T19:50:08.8Z','AntTG_ocean_height_v1.txt?dataset_id=601358']
    local = os.path.join(filepath,'AntTG_ocean_height_v1.txt')
    pyTMD.utilities.from_http(HOST,local=local)

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

#-- PURPOSE: Tests that interpolated results are comparable to AntTG database
def test_compare_CATS2008():
    #-- model parameters for CATS2008
    grid_file = os.path.join(filepath,'grid_CATS2008')
    model_file = os.path.join(filepath,'hf.CATS2008.out')
    GRID = 'OTIS'
    EPSG = 'CATS2008'
    TYPE = 'z'

    #-- open Antarctic Tide Gauge (AntTG) database
    with open(os.path.join(filepath,'AntTG_ocean_height_v1.txt'),'r') as f:
        file_contents = f.read().splitlines()
    #-- counts the number of lines in the header
    count = 0
    HEADER = True
    #-- Reading over header text
    while HEADER:
        #-- check if file line at count starts with matlab comment string
        HEADER = file_contents[count].startswith('%')
        #-- add 1 to counter
        count += 1
    #-- rewind 1 line
    count -= 1
    #-- iterate over number of stations
    AntTG = {}
    constituents = ['q1','o1','p1','k1','n2','m2','s2','k2']
    antarctic_stations = (len(file_contents) - count)//10
    stations = [None]*antarctic_stations
    shortname = [None]*antarctic_stations
    station_lon = np.zeros((antarctic_stations))
    station_lat = np.zeros((antarctic_stations))
    station_amp = np.ma.zeros((antarctic_stations,len(constituents)))
    station_ph = np.ma.zeros((antarctic_stations,len(constituents)))
    for s in range(antarctic_stations):
        i = count + s*10
        stations[s] = file_contents[i + 1].strip()
        shortname[s] = file_contents[i + 3].strip()
        lon,lat,aux1,aux2 = file_contents[i + 4].split()
        station_lon[s] = np.float(lon)
        station_lat[s] = np.float(lat)
        amp = file_contents[i + 7].split()
        ph = file_contents[i + 8].split()
        station_amp.data[s,:] = np.array(amp,dtype=np.float)
        station_ph.data[s,:] = np.array(ph,dtype=np.float)
    #-- update masks where NaN
    station_amp.mask = np.isnan(station_amp.data) | (station_amp.data == 0.0)
    station_ph.mask = np.isnan(station_ph.data)
    #-- replace nans with fill values
    station_amp.data[station_amp.mask] = station_amp.fill_value
    station_ph.data[station_ph.mask] = station_ph.fill_value

    #-- extract amplitude and phase from tide model
    amp,ph,D,cons = pyTMD.read_tide_model.extract_tidal_constants(station_lon,
        station_lat, grid_file, model_file, EPSG, TYPE=TYPE, METHOD='spline',
        GRID=GRID)
    #-- reorder constituents of model and convert amplitudes to cm
    model_amp = np.ma.zeros((antarctic_stations,len(constituents)))
    model_ph = np.ma.zeros((antarctic_stations,len(constituents)))
    for i,c in enumerate(constituents):
        j, = [j for j,val in enumerate(cons) if (val == c)]
        model_amp[:,i] = 100.0*amp[:,j]
        model_ph[:,i] = ph[:,j]
    #-- calculate complex constituent oscillations
    station_z = station_amp*np.exp(-1j*station_ph*np.pi/180.0)
    model_z = model_amp*np.exp(-1j*model_ph*np.pi/180.0)
    #-- valid stations for all constituents
    valid = np.all((~station_z.mask) & (~model_z.mask), axis=1)
    nv = np.count_nonzero(valid)
    #-- compare with RMS values from King et al. (2011)
    #-- https://doi.org/10.1029/2011JC006949
    RMS = np.array([1.4,2.7,1.7,3.5,2.9,7.3,5.0,1.7])
    rms = np.zeros((len(constituents)))
    for i,c in enumerate(constituents):
        #-- calculate difference and rms
        difference = np.abs(station_z[valid,i] - model_z[valid,i])
        rms[i] = np.sqrt(np.sum(difference**2))/(2.0*nv)
    #-- test RMS differences
    assert np.all(rms <= RMS)

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
