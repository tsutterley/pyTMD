#!/usr/bin/env python
u"""
test_download_and_read.py (08/2020)
Tests that CATS2008 data can be downloaded from the US Antarctic Program (USAP)
Tests that AOTIM-5-2018 data can be downloaded from the NSF ArcticData server
Tests the read program to verify that constituents are being extracted
Tests that interpolated results are comparable to Matlab TMD program
    https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    Oct2Py: Python to GNU Octave Bridge
        https://oct2py.readthedocs.io/en/latest/

UPDATE HISTORY:
    Updated 08/2020: Download Antarctic tide gauge database and compare with RMS
        directly call Matlab program (octave+oct2py) and compare outputs
        compare outputs for both Antarctic (CATS2008) and Arctic (AOTIM-5-2018)
        will install octave and oct2py in development requirements
    Written 08/2020
"""
import os
import re
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
from oct2py import octave

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
    #-- find model files within zip file
    rx = re.compile(r'(grid|h[0f]?|UV[0]?|Model|xy)[_\.](.*?)',re.IGNORECASE)
    m = [m for m in zfile.filelist if rx.match(posixpath.basename(m.filename))]
    #-- verify that model files are within downloaded zip file
    assert all(m)
    #-- extract each member (model and configuration files)
    for member in m:
        #-- strip directories from member filename
        member.filename = posixpath.basename(member.filename)
        print('\t{0}\n'.format(os.path.join(filepath,member.filename)))
        zfile.extract(member, path=filepath)
    #-- close the zipfile object
    zfile.close()

#-- PURPOSE: Download AOTIM-5-2018 from NSF ArcticData server
def test_download_AOTIM5_2018():
    #-- build host url for model
    resource_map_doi = 'resource_map_doi:{0}'.format('10.18739/A21R6N14K')
    HOST = ['https://arcticdata.io','metacat','d1','mn','v2','packages',
        pyTMD.utilities.quote_plus(posixpath.join('application','bagit-097')),
        pyTMD.utilities.quote_plus(resource_map_doi)]
    #-- download zipfile from host
    FILE = pyTMD.utilities.from_http(HOST)
    zfile = zipfile.ZipFile(FILE)
    print('{0} -->\n'.format(posixpath.join(*HOST)))
    #-- find model files within zip file
    rx = re.compile(r'(grid|h[0f]?|UV[0]?|Model|xy)[_\.](.*?)',re.IGNORECASE)
    m = [m for m in zfile.filelist if rx.match(posixpath.basename(m.filename))]
    #-- verify that model files are within downloaded zip file
    assert all(m)
    #-- extract each member (model and configuration files)
    for member in m:
        #-- strip directories from member filename
        member.filename = posixpath.basename(member.filename)
        print('\t{0}\n'.format(os.path.join(filepath,member.filename)))
        #-- extract file
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
    assert os.access(local, os.F_OK)

#-- PURPOSE: Download Arctic Tidal Current Atlas list of records
def test_download_Arctic_Tide_Atlas():
    HOST = ['https://arcticdata.io','metacat','d1','mn','v2','object',
        'urn%3Auuid%3Ae3abe2cc-f903-44de-9758-0c6bfc5b66c9']
    local = os.path.join(filepath,'List_of_records.txt')
    pyTMD.utilities.from_http(HOST,local=local)
    assert os.access(local, os.F_OK)

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
        lat,lon,aux1,aux2 = file_contents[i + 4].split()
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
    invalid_list = ['Ablation Lake','Amery','Bahia Esperanza','Beaver Lake',
        'Cape Roberts','Casey','Doake Ice Rumples','EE4A','EE4B',
        'Eklund Islands','Gerlache C','Groussac','Gurrachaga','Half Moon Is.',
        'Heard Island','Hobbs Pool','Mawson','McMurdo','Mikkelsen','Palmer',
        'Primavera','Rutford GL','Rutford GPS','Rothera','Scott Base',
        'Seymour Is','Terra Nova Bay']
    #-- remove coastal stations from the list
    invalid_stations = [i for i,s in enumerate(shortname) if s in invalid_list]
    valid[invalid_stations] = False
    nv = np.count_nonzero(valid)
    #-- compare with RMS values from King et al. (2011)
    #-- https://doi.org/10.1029/2011JC006949
    RMS = np.array([1.4,2.7,1.7,3.5,2.9,7.3,5.0,1.7])
    rms = np.zeros((len(constituents)))
    for i,c in enumerate(constituents):
        #-- calculate difference and rms
        difference = np.abs(station_z[valid,i] - model_z[valid,i])
        #-- round to precision of King et al. (2011)
        rms[i] = np.round(np.sqrt(np.sum(difference**2)/(2.0*nv)),decimals=1)
    #-- test RMS differences
    assert np.all(rms <= RMS)

#-- PURPOSE: calculate the matlab serial date from calendar date
#-- http://scienceworld.wolfram.com/astronomy/JulianDate.html
def convert_calendar_serial(year, month, day, hour=0.0, minute=0.0, second=0.0):
    #-- return the date in days since serial epoch 0000-01-01T00:00:00
    sd = 367.0*year - np.floor(7.0*(year + np.floor((month+9.0)/12.0))/4.0) - \
        np.floor(3.0*(np.floor((year + (month - 9.0)/7.0)/100.0) + 1.0)/4.0) + \
        np.floor(275.0*month/9.0) + day + hour/24.0 + minute/1440.0 + \
        second/86400.0 - 30.0
    return sd

#-- PURPOSE: Tests that interpolated results are comparable to Matlab program
def test_verify_CATS2008():
    #-- model parameters for CATS2008
    grid_file = os.path.join(filepath,'grid_CATS2008')
    elevation_file = os.path.join(filepath,'hf.CATS2008.out')
    transport_file = os.path.join(filepath,'uv.CATS2008.out')
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
    antarctic_stations = (len(file_contents) - count)//10
    stations = [None]*antarctic_stations
    shortname = [None]*antarctic_stations
    station_type = [None]*antarctic_stations
    station_lon = np.zeros((antarctic_stations))
    station_lat = np.zeros((antarctic_stations))
    for s in range(antarctic_stations):
        i = count + s*10
        stations[s] = file_contents[i + 1].strip()
        shortname[s] = file_contents[i + 3].strip()
        lat,lon,aux1,aux2 = file_contents[i + 4].split()
        station_type[s] = file_contents[i + 6].strip()
        station_lon[s] = np.float(lon)
        station_lat[s] = np.float(lat)

    #-- calculate daily results for a time period
    #-- convert time to days since 1992-01-01T00:00:00
    tide_time = np.arange(pyTMD.time.convert_calendar_dates(2000,1,1),
        pyTMD.time.convert_calendar_dates(2000,12,31)+1)
    #-- serial dates for matlab program (days since 0000-01-01T00:00:00)
    SDtime = np.arange(convert_calendar_serial(2000,1,1),
        convert_calendar_serial(2000,12,31)+1)
    #-- presently not converting times to dynamic times for model comparisons
    deltat = np.zeros_like(tide_time)
    #-- number of days
    ndays = len(tide_time)

    #-- extract amplitude and phase from tide model
    amp,ph,D,c = pyTMD.read_tide_model.extract_tidal_constants(station_lon,
        station_lat, grid_file, elevation_file, EPSG, TYPE=TYPE,
        METHOD='spline', GRID=GRID)
    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- will verify differences between model outputs are within tolerance
    eps = np.finfo(np.float16).eps

    #-- compare daily outputs at each station point
    invalid_list = ['Ablation Lake','Amery','Bahia Esperanza','Beaver Lake',
        'Cape Roberts','Casey','Doake Ice Rumples','EE4A','EE4B',
        'Eklund Islands','Gerlache C','Groussac','Gurrachaga','Half Moon Is.',
        'Heard Island','Hobbs Pool','Mawson','McMurdo','Mikkelsen','Palmer',
        'Primavera','Rutford GL','Rutford GPS','Rothera','Scott Base',
        'Seymour Is','Terra Nova Bay']
    #-- remove coastal stations from the list
    valid_stations=[i for i,s in enumerate(shortname) if s not in invalid_list]
    for i in valid_stations:
        #-- calculate constituent oscillation for station
        hc = amp[i,None,:]*np.exp(cph[i,None,:])
        #-- allocate for out tides at point
        tide = np.ma.zeros((ndays))
        tide.mask = np.zeros((ndays),dtype=np.bool)
        #-- predict tidal elevations at time and infer minor corrections
        tide.mask[:] = np.any(hc.mask)
        tide.data[:] = pyTMD.predict_tidal_ts(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=GRID)
        minor = pyTMD.infer_minor_corrections(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=GRID)
        tide.data[:] += minor.data[:]

        #-- compute validation data from Matlab TMD program using octave
        #-- https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5
        TMDpath = os.path.join(filepath,'..','TMD_Matlab_Toolbox','TMD')
        octave.addpath(octave.genpath(os.path.normpath(TMDpath)))
        octave.addpath(filepath)
        octave.warning('off', 'all')
        CFname = os.path.join(filepath,'Model_CATS2008')
        validation,cons = octave.tmd_tide_pred(CFname,SDtime,
            station_lat[i],station_lon[i],'z',nout=2)

        #-- calculate differences between matlab and python version
        difference = np.ma.zeros((ndays))
        difference.data[:] = tide.data - validation.T
        difference.mask = (tide.mask | np.isnan(validation))
        if not np.all(difference.mask):
            assert np.all(np.abs(difference) < eps)

#-- PURPOSE: Tests that interpolated results are comparable to Matlab program
def test_verify_AOTIM5_2018():
    #-- model parameters for AOTIM-5-2018
    grid_file = os.path.join(filepath,'grid_Arc5km2018')
    elevation_file = os.path.join(filepath,'h_Arc5km2018')
    transport_file = os.path.join(filepath,'UV_Arc5km2018')
    GRID = 'OTIS'
    EPSG = 'PSNorth'
    TYPE = 'z'

    #-- open Arctic Tidal Current Atlas list of records
    with open(os.path.join(filepath,'List_of_records.txt'),'r') as f:
        file_contents = f.read().splitlines()
    #-- skip 2 header rows
    count = 2
    #-- iterate over number of stations
    arctic_stations = len(file_contents) - count
    stations = [None]*arctic_stations
    shortname = [None]*arctic_stations
    station_lon = np.zeros((arctic_stations))
    station_lat = np.zeros((arctic_stations))
    for s in range(arctic_stations):
        line_contents = file_contents[count+s].split()
        stations[s] = line_contents[1]
        shortname[s] = line_contents[2]
        station_lat[s] = np.float(line_contents[10])
        station_lon[s] = np.float(line_contents[11])

    #-- calculate daily results for a time period
    #-- convert time to days since 1992-01-01T00:00:00
    tide_time = np.arange(pyTMD.time.convert_calendar_dates(2000,1,1),
        pyTMD.time.convert_calendar_dates(2000,12,31)+1)
    #-- serial dates for matlab program (days since 0000-01-01T00:00:00)
    SDtime = np.arange(convert_calendar_serial(2000,1,1),
        convert_calendar_serial(2000,12,31)+1)
    #-- presently not converting times to dynamic times for model comparisons
    deltat = np.zeros_like(tide_time)
    #-- number of days
    ndays = len(tide_time)

    #-- extract amplitude and phase from tide model
    amp,ph,D,c = pyTMD.read_tide_model.extract_tidal_constants(station_lon,
        station_lat, grid_file, elevation_file, EPSG, TYPE=TYPE,
        METHOD='spline', GRID=GRID)
    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- will verify differences between model outputs are within tolerance
    eps = np.finfo(np.float16).eps

    #-- compare daily outputs at each station point
    invalid_list = ['KS14']
    #-- remove coastal stations from the list
    valid_stations=[i for i,s in enumerate(shortname) if s not in invalid_list]
    for i in valid_stations:
        #-- calculate constituent oscillation for station
        hc = amp[i,None,:]*np.exp(cph[i,None,:])
        #-- allocate for out tides at point
        tide = np.ma.zeros((ndays))
        tide.mask = np.zeros((ndays),dtype=np.bool)
        #-- predict tidal elevations at time and infer minor corrections
        tide.mask[:] = np.any(hc.mask)
        tide.data[:] = pyTMD.predict_tidal_ts(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=GRID)
        minor = pyTMD.infer_minor_corrections(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=GRID)
        tide.data[:] += minor.data[:]

        #-- compute validation data from Matlab TMD program using octave
        #-- https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5
        TMDpath = os.path.join(filepath,'..','TMD_Matlab_Toolbox','TMD')
        octave.addpath(octave.genpath(os.path.normpath(TMDpath)))
        octave.addpath(filepath)
        octave.warning('off', 'all')
        CFname = os.path.join(filepath,'Model_Arc5km2018')
        validation,cons = octave.tmd_tide_pred(CFname,SDtime,
            station_lat[i],station_lon[i],'z',nout=2)

        #-- calculate differences between matlab and python version
        difference = np.ma.zeros((ndays))
        difference.data[:] = tide.data - validation.T
        difference.mask = (tide.mask | np.isnan(validation))
        if not np.all(difference.mask):
            assert np.all(np.abs(difference) < eps)
