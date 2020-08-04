#!/usr/bin/env python
u"""
compute_tides_icebridge_data.py
Written by Tyler Sutterley (07/2020)
Calculates tidal elevations for correcting Operation IceBridge elevation data

Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
or Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

INPUTS:
    ATM1B, ATM icessn or LVIS file from NSIDC

COMMAND LINE OPTIONS:
    -D X, --directory=X: Working data directory
    -T X, --tide=X: Tide model to use in correction
        CATS0201
        CATS2008
        CATS2008_load
        TPXO9-atlas-v2
        TPXO9-atlas
        TPXO9.1
        TPXO8-atlas
        TPXO7.2
        TPXO7.2_load
        AODTM-5
        AOTIM-5
        AOTIM-5-2018
        GOT4.7
        GOT4.7_load
        GOT4.8
        GOT4.8_load
        GOT4.10
        GOT4.10_load
        FES2014
        FES2014_load
    -M X, --mode=X: Permission mode of directories and files created
    -V, --verbose: Output information about each created file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    utilities: download and management utilities for syncing files
    convert_julian.py: returns the calendar date and time given a Julian date
    calc_astrol_longitudes.py: computes the basic astronomical mean longitudes
    calc_delta_time.py: calculates difference between universal and dynamic time
    convert_ll_xy.py: convert lat/lon points to and from projected coordinates
    infer_minor_corrections.py: return corrections for minor constituents
    load_constituent.py: loads parameters for a given tidal constituent
    load_nodal_corrections.py: load the nodal corrections for tidal constituents
    predict_tide_drift.py: predict tidal elevations using harmonic constants
    read_tide_model.py: extract tidal harmonic constants from OTIS tide models
    read_netcdf_model.py: extract tidal harmonic constants from netcdf models
    read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
    read_FES_model.py: extract tidal harmonic constants from FES tide models
    read_ATM1b_QFIT_binary.py: read ATM1b QFIT binary files (NSIDC version 1)

UPDATE HISTORY:
    Updated 08/2020: using builtin time operations
    Updated 07/2020: added FES2014 and FES2014_load.  use merged delta times
    Updated 06/2020: added version 2 of TPX09-atlas (TPX09-atlas-v2)
    Updated 03/2020: use read_ATM1b_QFIT_binary from repository
    Updated 02/2020: changed CATS2008 grid to match version on U.S. Antarctic
        Program Data Center http://www.usap-dc.org/view/dataset/601235
    Updated 11/2019: added AOTIM-5-2018 tide model (2018 update to 2004 model)
    Updated 09/2019: added TPXO9_atlas reading from netcdf4 tide files
    Updated 05/2019: added option interpolate to choose the interpolation method
    Updated 02/2019: using range for python3 compatibility
    Updated 10/2018: updated GPS time calculation for calculating leap seconds
    Updated 07/2018: added GSFC Global Ocean Tides (GOT) models
    Written 06/2018
"""
from __future__ import print_function

import sys
import os
import re
import time
import h5py
import getopt
import numpy as np
import pyTMD.time
import pyTMD.utilities
from pyTMD.convert_julian import convert_julian
from read_ATM1b_QFIT_binary.read_ATM1b_QFIT_binary import read_ATM1b_QFIT_binary
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide_drift import predict_tide_drift
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants
from pyTMD.read_FES_model import extract_FES_constants

#-- PURPOSE: reading the number of file lines removing commented lines
def file_length(input_file, input_subsetter, HDF5=False, QFIT=False):
    #-- subset the data to indices if specified
    if input_subsetter:
        file_lines = len(input_subsetter)
    elif HDF5:
        #-- read the size of an input variable within a HDF5 file
        with h5py.File(input_file,'r') as fileID:
            file_lines, = fileID[HDF5].shape
    elif QFIT:
        #-- read the size of a QFIT binary file
        file_lines = read_ATM1b_QFIT_binary.ATM1b_QFIT_shape(input_file)
    else:
        #-- read the input file, split at lines and remove all commented lines
        with open(input_file,'r') as f:
            i = [i for i in f.read().splitlines() if re.match('^(?!#)',i)]
        file_lines = len(i)
    #-- return the number of lines
    return file_lines

##-- PURPOSE: read the ATM Level-1b data file for variables of interest
def read_ATM_qfit_file(input_file, input_subsetter):
    #-- regular expression pattern for extracting parameters
    mission_flag = '(BLATM1B|ILATM1B|ILNSA1B)'
    regex_pattern = '{0}_(\d+)_(\d+)(.*?).(qi|TXT|h5)'.format(mission_flag)
    #-- extract mission and other parameters from filename
    MISSION,YYMMDD,HHMMSS,AUX,SFX = re.findall(regex_pattern,input_file).pop()
    #-- early date strings omitted century and millenia (e.g. 93 for 1993)
    if (len(YYMMDD) == 6):
        ypre,month,day = np.array([YYMMDD[:2],YYMMDD[2:4],YYMMDD[4:]],dtype='i')
        year = (ypre + 1900.0) if (ypre >= 90) else (ypre + 2000.0)
    elif (len(YYMMDD) == 8):
        year,month,day = np.array([YYMMDD[:4],YYMMDD[4:6],YYMMDD[6:]],dtype='i')
    #-- output python dictionary with variables
    ATM_L1b_input = {}
    #-- Version 1 of ATM QFIT files (ascii)
    #-- output text file from qi2txt with proper filename format
    #-- do not use the shortened output format from qi2txt
    if (SFX == 'TXT'):
        #-- compile regular expression operator for reading lines
        regex_pattern = '[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)
        #-- read the input file, split at lines and remove all commented lines
        f = open(input_file,'r')
        file_contents=[i for i in f.read().splitlines() if re.match('^(?!#)',i)]
        f.close()
        #-- number of lines of data within file
        file_lines = file_length(input_file,input_subsetter)
        #-- create output variables with length equal to the number of lines
        ATM_L1b_input['lat'] = np.zeros_like(file_contents,dtype=np.float)
        ATM_L1b_input['lon'] = np.zeros_like(file_contents,dtype=np.float)
        ATM_L1b_input['data'] = np.zeros_like(file_contents,dtype=np.float)
        hour = np.zeros_like(file_contents,dtype=np.float)
        minute = np.zeros_like(file_contents,dtype=np.float)
        second = np.zeros_like(file_contents,dtype=np.float)
        #-- for each line within the file
        for i,line in enumerate(file_contents):
            #-- find numerical instances within the line
            line_contents = rx.findall(line)
            ATM_L1b_input['lat'][i] = np.float(line_contents[1])
            ATM_L1b_input['lon'][i] = np.float(line_contents[2])
            ATM_L1b_input['data'][i] = np.float(line_contents[3])
            hour[i] = np.float(line_contents[-1][:2])
            minute[i] = np.float(line_contents[-1][2:4])
            second[i] = np.float(line_contents[-1][4:])
    #-- Version 1 of ATM QFIT files (binary)
    elif (SFX == 'qi'):
        #-- read input QFIT data file and subset if specified
        fid,h = read_ATM1b_QFIT_binary(input_file)
        #-- number of lines of data within file
        file_lines = file_length(input_file,input_subsetter,QFIT=True)
        ATM_L1b_input['lat'] = fid['latitude'][:]
        ATM_L1b_input['lon'] = fid['longitude'][:]
        ATM_L1b_input['data'] = fid['elevation'][:]
        time_hhmmss = fid['time_hhmmss'][:]
        #-- extract hour, minute and second from time_hhmmss
        hour = np.zeros_like(time_hhmmss,dtype=np.float)
        minute = np.zeros_like(time_hhmmss,dtype=np.float)
        second = np.zeros_like(time_hhmmss,dtype=np.float)
        #-- for each line within the file
        for i,packed_time in enumerate(time_hhmmss):
            #-- convert to zero-padded string with 3 decimal points
            line_contents = '{0:010.3f}'.format(packed_time)
            hour[i] = np.float(line_contents[:2])
            minute[i] = np.float(line_contents[2:4])
            second[i] = np.float(line_contents[4:])
    #-- Version 2 of ATM QFIT files (HDF5)
    elif (SFX == 'h5'):
        #-- Open the HDF5 file for reading
        fileID = h5py.File(os.path.expanduser(input_file), 'r')
        #-- number of lines of data within file
        file_lines = file_length(input_file,input_subsetter,HDF5='elevation')
        #-- create output variables with length equal to input elevation
        ATM_L1b_input['lat'] = fileID['latitude'][:]
        ATM_L1b_input['lon'] = fileID['longitude'][:]
        ATM_L1b_input['data'] = fileID['elevation'][:]
        time_hhmmss = fileID['instrument_parameters']['time_hhmmss'][:]
        #-- extract hour, minute and second from time_hhmmss
        hour = np.zeros_like(time_hhmmss,dtype=np.float)
        minute = np.zeros_like(time_hhmmss,dtype=np.float)
        second = np.zeros_like(time_hhmmss,dtype=np.float)
        #-- for each line within the file
        for i,packed_time in enumerate(time_hhmmss):
            #-- convert to zero-padded string with 3 decimal points
            line_contents = '{0:010.3f}'.format(packed_time)
            hour[i] = np.float(line_contents[:2])
            minute[i] = np.float(line_contents[2:4])
            second[i] = np.float(line_contents[4:])
        #-- close the input HDF5 file
        fileID.close()
    #-- calculate the number of leap seconds between GPS time (seconds
    #-- since Jan 6, 1980 00:00:00) and UTC
    gps_seconds = pyTMD.time.convert_calendar_dates(year,month,day,
        hour=hour,minute=minute,second=second,
        epoch=(1980,1,6,0,0,0),scale=86400.0)
    leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)
    #-- calculation of Julian day taking into account leap seconds
    #-- converting to J2000 seconds
    ATM_L1b_input['time'] = pyTMD.time.convert_calendar_dates(year,month,day,
        hour=hour,minute=minute,second=second-leap_seconds,
        epoch=(2000,1,1,12,0,0,0),scale=86400.0)
    #-- subset the data to indices if specified
    if input_subsetter:
        for key,val in ATM_L1b_input.items():
            ATM_L1b_input[key] = val[input_subsetter]
    #-- hemispheric shot count
    count = {}
    count['N'] = np.count_nonzero(ATM_L1b_input['lat'] >= 0.0)
    count['S'] = np.count_nonzero(ATM_L1b_input['lat'] < 0.0)
    #-- determine hemisphere with containing shots in file
    HEM, = [key for key, val in count.items() if val]
    #-- return the output variables
    return ATM_L1b_input,file_lines,HEM

#-- PURPOSE: read the ATM Level-2 data file for variables of interest
def read_ATM_icessn_file(input_file, input_subsetter):
    #-- regular expression pattern for extracting parameters
    regex_pattern = '(BLATM2|ILATM2)_(\d+)_(\d+)_smooth_nadir(.*?)(csv|seg|pt)$'
    #-- extract mission and other parameters from filename
    MISSION,YYMMDD,HHMMSS,AUX,SFX = re.findall(regex_pattern,input_file).pop()
    #-- early date strings omitted century and millenia (e.g. 93 for 1993)
    if (len(YYMMDD) == 6):
        ypre,month,day = np.array([YYMMDD[:2],YYMMDD[2:4],YYMMDD[4:]],dtype='i')
        year = (ypre + 1900.0) if (ypre >= 90) else (ypre + 2000.0)
    elif (len(YYMMDD) == 8):
        year,month,day = np.array([YYMMDD[:4],YYMMDD[4:6],YYMMDD[6:]],dtype='i')
    #-- input file column names for variables of interest with column indices
    #-- variables not used: (SNslope:4, WEslope:5, npt_used:7, npt_edit:8, d:9)
    file_dtype = {'seconds':0, 'lat':1, 'lon':2, 'data':3, 'RMS':6, 'track':-1}
    #-- compile regular expression operator for reading lines (extracts numbers)
    regex_pattern = '[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)
    #-- read the input file, split at lines and remove all commented lines
    with open(input_file,'r') as f:
        file_contents=[i for i in f.read().splitlines() if re.match('^(?!#)',i)]
    #-- number of lines of data within file
    file_lines = file_length(input_file,input_subsetter)
    #-- output python dictionary with variables
    ATM_L2_input = {}
    #-- create output variables with length equal to the number of file lines
    for key in file_dtype.keys():
        ATM_L2_input[key] = np.zeros_like(file_contents, dtype=np.float)
    #-- for each line within the file
    for line_number,line_entries in enumerate(file_contents):
        #-- find numerical instances within the line
        line_contents = rx.findall(line_entries)
        #-- for each variable of interest: save to dinput as float
        for key,val in file_dtype.items():
            ATM_L2_input[key][line_number] = np.float(line_contents[val])
    #-- convert shot time (seconds of day) to J2000
    hour = np.floor(ATM_L2_input['seconds']/3600.0)
    minute = np.floor((ATM_L2_input['seconds'] % 3600)/60.0)
    second = ATM_L2_input['seconds'] % 60.0
    #-- First column in Pre-IceBridge and ICESSN Version 1 files is GPS time
    if (MISSION == 'BLATM2') or (SFX != 'csv'):
        #-- calculate the number of leap seconds between GPS time (seconds
        #-- since Jan 6, 1980 00:00:00) and UTC
        gps_seconds = pyTMD.time.convert_calendar_dates(year,month,day,
            hour=hour,minute=minute,second=second,
            epoch=(1980,1,6,0,0,0),scale=86400.0)
        leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)
    else:
        leap_seconds = 0.0
    #-- calculation of Julian day
    #-- converting to J2000 seconds
    ATM_L2_input['time'] = pyTMD.time.convert_calendar_dates(year,month,day,
        hour=hour,minute=minute,second=second-leap_seconds,
        epoch=(2000,1,1,12,0,0,0),scale=86400.0)
    #-- convert RMS from centimeters to meters
    ATM_L2_input['error'] = ATM_L2_input['RMS']/100.0
    #-- subset the data to indices if specified
    if input_subsetter:
        for key,val in ATM_L2_input.items():
            ATM_L2_input[key] = val[input_subsetter]
    #-- hemispheric shot count
    count = {}
    count['N'] = np.count_nonzero(ATM_L2_input['lat'] >= 0.0)
    count['S'] = np.count_nonzero(ATM_L2_input['lat'] < 0.0)
    #-- determine hemisphere with containing shots in file
    HEM, = [key for key, val in count.items() if val]
    #-- return the output variables
    return ATM_L2_input,file_lines,HEM

#-- PURPOSE: read the LVIS Level-2 data file for variables of interest
def read_LVIS_HDF5_file(input_file, input_subsetter):
    #-- LVIS region flags: GL for Greenland and AQ for Antarctica
    lvis_flag = {'GL':'N','AQ':'S'}
    #-- regular expression pattern for extracting parameters from HDF5 files
    #-- computed in read_icebridge_lvis.py
    mission_flag = '(BLVIS2|BVLIS2|ILVIS2|ILVGH2)'
    regex_pattern = '{0}_(.*?)(\d+)_(\d+)_(R\d+)_(\d+).H5'.format(mission_flag)
    #-- extract mission, region and other parameters from filename
    MISSION,REGION,YY,MMDD,RLD,SS = re.findall(regex_pattern,input_file).pop()
    LDS_VERSION = '2.0.2' if (np.int(RLD[1:3]) >= 18) else '1.04'
    #-- input and output python dictionaries with variables
    file_input = {}
    LVIS_L2_input = {}
    fileID = h5py.File(input_file,'r')
    #-- create output variables with length equal to input shot number
    file_lines = file_length(input_file,input_subsetter,HDF5='Shot_Number')
    #-- https://lvis.gsfc.nasa.gov/Data/Data_Structure/DataStructure_LDS104.html
    #-- https://lvis.gsfc.nasa.gov/Data/Data_Structure/DataStructure_LDS202.html
    if (LDS_VERSION == '1.04'):
        #-- elevation surfaces
        file_input['elev'] = fileID['Elevation_Surfaces/Elevation_Centroid'][:]
        file_input['elev_low'] = fileID['Elevation_Surfaces/Elevation_Low'][:]
        file_input['elev_high'] = fileID['Elevation_Surfaces/Elevation_High'][:]
        #-- latitude
        file_input['lat'] = fileID['Geolocation/Latitude_Centroid'][:]
        file_input['lat_low'] = fileID['Geolocation/Latitude_Low'][:]
        #-- longitude
        file_input['lon'] = fileID['Geolocation/Longitude_Centroid'][:]
        file_input['lon_low'] = fileID['Geolocation/Longitude_Low'][:]
    elif (LDS_VERSION == '2.0.2'):
        #-- elevation surfaces
        file_input['elev_low'] = fileID['Elevation_Surfaces/Elevation_Low'][:]
        file_input['elev_high'] = fileID['Elevation_Surfaces/Elevation_High'][:]
        #-- heights above lowest detected mode
        file_input['RH50'] = fileID['Waveform/RH50'][:]
        file_input['RH100'] = fileID['Waveform/RH100'][:]
        #-- calculate centroidal elevation using 50% of waveform energy
        file_input['elev'] = file_input['elev_low'] + file_input['RH50']
        #-- latitude
        file_input['lat_top'] = fileID['Geolocation/Latitude_Top'][:]
        file_input['lat_low'] = fileID['Geolocation/Latitude_Low'][:]
        #-- longitude
        file_input['lon_top'] = fileID['Geolocation/Longitude_Top'][:]
        file_input['lon_low'] = fileID['Geolocation/Longitude_Low'][:]
        #-- linearly interpolate latitude and longitude to RH50
        file_input['lat'] = file_input['lat_low'] + file_input['RH50'] * \
            (file_input['lat_top'] - file_input['lat_low'])/file_input['RH100']
        file_input['lon'] = file_input['lon_low'] + file_input['RH50'] * \
            (file_input['lon_top'] - file_input['lon_low'])/file_input['RH100']
    #-- J2000 seconds
    LVIS_L2_input['time'] = fileID['Time/J2000'][:]
    #-- close the input HDF5 file
    fileID.close()
    #-- output combined variables
    LVIS_L2_input['data'] = np.zeros_like(file_input['elev'],dtype=np.float)
    LVIS_L2_input['lon'] = np.zeros_like(file_input['elev'],dtype=np.float)
    LVIS_L2_input['lat'] = np.zeros_like(file_input['elev'],dtype=np.float)
    LVIS_L2_input['error'] = np.zeros_like(file_input['elev'],dtype=np.float)
    #-- find where elev high is equal to elev low
    #-- see note about using LVIS centroid elevation product
    #-- http://lvis.gsfc.nasa.gov/OIBDataStructure.html
    ii = np.nonzero(file_input['elev_low'] == file_input['elev_high'])
    jj = np.nonzero(file_input['elev_low'] != file_input['elev_high'])
    #-- where lowest point of waveform is equal to highest point -->
    #-- using the elev_low elevation
    LVIS_L2_input['data'][ii] = file_input['elev_low'][ii]
    #-- for other locations use the centroid elevation
    #-- as the centroid is a useful product over rough terrain
    #-- when you are calculating ice volume change
    LVIS_L2_input['data'][jj] = file_input['elev'][jj]
    #-- latitude and longitude for each case
    #-- elevation low == elevation high
    LVIS_L2_input['lon'][ii] = file_input['lon_low'][ii]
    LVIS_L2_input['lat'][ii] = file_input['lat_low'][ii]
    #-- centroid elevations
    LVIS_L2_input['lon'][jj] = file_input['lon'][jj]
    LVIS_L2_input['lat'][jj] = file_input['lat'][jj]
    #-- estimated uncertainty for both cases
    LVIS_variance_low = (file_input['elev_low'] - file_input['elev'])**2
    LVIS_variance_high = (file_input['elev_high'] - file_input['elev'])**2
    LVIS_L2_input['error']=np.sqrt((LVIS_variance_low + LVIS_variance_high)/2.0)
    #-- subset the data to indices if specified
    if input_subsetter:
        for key,val in LVIS_L2_input.items():
            LVIS_L2_input[key] = val[input_subsetter]
    #-- return the output variables
    return LVIS_L2_input,file_lines,lvis_flag[REGION]

#-- PURPOSE: read Operation IceBridge data from NSIDC
#-- compute tides at points and times using tidal model driver algorithms
def compute_tides_icebridge_data(tide_dir, arg, MODEL, METHOD=None,
    VERBOSE=False, MODE=0o775):

    #-- extract file name and subsetter indices lists
    match_object = re.match('(.*?)(\[(.*?)\])?$',arg)
    input_file = os.path.expanduser(match_object.group(1))
    #-- subset input file to indices
    if match_object.group(2):
        #-- decompress ranges and add to list
        input_subsetter = []
        for i in re.findall('((\d+)-(\d+)|(\d+))',match_object.group(3)):
            input_subsetter.append(int(i[3])) if i[3] else \
                input_subsetter.extend(range(int(i[1]),int(i[2])+1))
    else:
        input_subsetter = None

    #-- output directory for input_file
    DIRECTORY = os.path.dirname(input_file)
    #-- calculate if input files are from ATM or LVIS (+GH)
    regex = {}
    regex['ATM'] = '(BLATM2|ILATM2)_(\d+)_(\d+)_smooth_nadir(.*?)(csv|seg|pt)$'
    regex['ATM1b'] = '(BLATM1b|ILATM1b)_(\d+)_(\d+)(.*?).(qi|TXT|h5)$'
    regex['LVIS'] = '(BLVIS2|BVLIS2|ILVIS2)_(.*?)(\d+)_(\d+)_(R\d+)_(\d+).H5$'
    regex['LVGH'] = '(ILVGH2)_(.*?)(\d+)_(\d+)_(R\d+)_(\d+).H5$'
    for key,val in regex.items():
        if re.match(val, os.path.basename(input_file)):
            OIB = key

    #-- HDF5 file attributes
    attrib = dict(lon={},lat={},tide={},day={})
    #-- select between tide models
    if (MODEL == 'CATS0201'):
        grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
        model_file = os.path.join(tide_dir,'cats0201_tmd','h0_CATS02_01')
        reference = 'https://mail.esr.org/polar_tide_models/Model_CATS0201.html'
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'CATS2008'):
        grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008')
        model_file = os.path.join(tide_dir,'CATS2008','hf.CATS2008.out')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = 'CATS2008'
        TYPE = 'z'
    elif (MODEL == 'CATS2008_load'):
        grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
        model_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        attrib['tide']['long_name'] = 'Load_Tide'
        model_format = 'OTIS'
        EPSG = 'CATS2008'
        TYPE = 'z'
    elif (MODEL == 'TPXO9-atlas'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas')
        grid_file = 'grid_tpxo9_atlas.nc.gz'
        model_files = ['h_q1_tpxo9_atlas_30.nc.gz','h_o1_tpxo9_atlas_30.nc.gz',
            'h_p1_tpxo9_atlas_30.nc.gz','h_k1_tpxo9_atlas_30.nc.gz',
            'h_n2_tpxo9_atlas_30.nc.gz','h_m2_tpxo9_atlas_30.nc.gz',
            'h_s2_tpxo9_atlas_30.nc.gz','h_k2_tpxo9_atlas_30.nc.gz',
            'h_m4_tpxo9_atlas_30.nc.gz','h_ms4_tpxo9_atlas_30.nc.gz',
            'h_mn4_tpxo9_atlas_30.nc.gz','h_2n2_tpxo9_atlas_30.nc.gz']
        reference = 'http://volkov.oce.orst.edu/tides/tpxo9_atlas.html'
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'netcdf'
        TYPE = 'z'
        SCALE = 1.0/1000.0
    elif (MODEL == 'TPXO9-atlas-v2'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas_v2')
        grid_file = 'grid_tpxo9_atlas_v2.nc.gz'
        model_files = ['h_q1_tpxo9_atlas_30_v2.nc.gz','h_o1_tpxo9_atlas_30_v2.nc.gz',
            'h_p1_tpxo9_atlas_30_v2.nc.gz','h_k1_tpxo9_atlas_30_v2.nc.gz',
            'h_n2_tpxo9_atlas_30_v2.nc.gz','h_m2_tpxo9_atlas_30_v2.nc.gz',
            'h_s2_tpxo9_atlas_30_v2.nc.gz','h_k2_tpxo9_atlas_30_v2.nc.gz',
            'h_m4_tpxo9_atlas_30_v2.nc.gz','h_ms4_tpxo9_atlas_30_v2.nc.gz',
            'h_mn4_tpxo9_atlas_30_v2.nc.gz','h_2n2_tpxo9_atlas_30_v2.nc.gz']
        reference = 'https://www.tpxo.net/global/tpxo9-atlas'
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'netcdf'
        TYPE = 'z'
        SCALE = 1.0/1000.0
    elif (MODEL == 'TPXO9.1'):
        grid_file = os.path.join(tide_dir,'TPXO9.1','DATA','grid_tpxo9')
        model_file = os.path.join(tide_dir,'TPXO9.1','DATA','h_tpxo9.v1')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_file = os.path.join(tide_dir,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
        reference = 'http://volkov.oce.orst.edu/tides/tpxo8_atlas.html'
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'ATLAS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'TPXO7.2'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_tmd','h_tpxo7.2')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'TPXO7.2_load'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_load','grid_tpxo6.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_load','h_tpxo7.2_load')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        attrib['tide']['long_name'] = 'Load_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'AODTM-5'):
        grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aodtm5_tmd','h0_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aodtm-5/')
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPE = 'z'
    elif (MODEL == 'AOTIM-5'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPE = 'z'
    elif (MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(tide_dir,'Arc5km2018','grid_Arc5km2018')
        model_file = os.path.join(tide_dir,'Arc5km2018','h_Arc5km2018')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPE = 'z'
    elif (MODEL == 'GOT4.7'):
        model_directory = os.path.join(tide_dir,'GOT4.7','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.7_load'):
        model_directory = os.path.join(tide_dir,'GOT4.7','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        attrib['tide']['long_name'] = 'Load_Tide'
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'GOT4.8'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.8_load'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        attrib['tide']['long_name'] = 'Load_Tide'
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'GOT4.10'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.10_load'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        attrib['tide']['long_name'] = 'Load_Tide'
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'FES2014'):
        model_directory = os.path.join(tide_dir,'fes2014','ocean_tide')
        model_files = ['2n2.nc.gz','eps2.nc.gz','j1.nc.gz','k1.nc.gz',
            'k2.nc.gz','l2.nc.gz','la2.nc.gz','m2.nc.gz','m3.nc.gz','m4.nc.gz',
            'm6.nc.gz','m8.nc.gz','mf.nc.gz','mks2.nc.gz','mm.nc.gz',
            'mn4.nc.gz','ms4.nc.gz','msf.nc.gz','msqm.nc.gz','mtm.nc.gz',
            'mu2.nc.gz','n2.nc.gz','n4.nc.gz','nu2.nc.gz','o1.nc.gz','p1.nc.gz',
            'q1.nc.gz','r2.nc.gz','s1.nc.gz','s2.nc.gz','s4.nc.gz','sa.nc.gz',
            'ssa.nc.gz','t2.nc.gz']
        c = ['2n2','eps2','j1','k1','k2','l2','lambda2','m2','m3','m4','m6',
            'm8','mf','mks2','mm','mn4','ms4','msf','msqm','mtm','mu2','n2',
            'n4','nu2','o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2']
        reference = ('https://www.aviso.altimetry.fr/data/products/'
            'auxiliary-products/global-tide-fes.html')
        attrib['tide']['long_name'] = 'Ocean_Tide'
        model_format = 'FES'
        TYPE = 'z'
        SCALE = 1.0/100.0
    elif (MODEL == 'FES2014_load'):
        model_directory = os.path.join(tide_dir,'fes2014','load_tide')
        model_files = ['2n2.nc.gz','eps2.nc.gz','j1.nc.gz','k1.nc.gz',
            'k2.nc.gz','l2.nc.gz','la2.nc.gz','m2.nc.gz','m3.nc.gz','m4.nc.gz',
            'm6.nc.gz','m8.nc.gz','mf.nc.gz','mks2.nc.gz','mm.nc.gz',
            'mn4.nc.gz','ms4.nc.gz','msf.nc.gz','msqm.nc.gz','mtm.nc.gz',
            'mu2.nc.gz','n2.nc.gz','n4.nc.gz','nu2.nc.gz','o1.nc.gz','p1.nc.gz',
            'q1.nc.gz','r2.nc.gz','s1.nc.gz','s2.nc.gz','s4.nc.gz','sa.nc.gz',
            'ssa.nc.gz','t2.nc.gz']
        c = ['2n2','eps2','j1','k1','k2','l2','lambda2','m2','m3','m4','m6',
            'm8','mf','mks2','mm','mn4','ms4','msf','msqm','mtm','mu2','n2',
            'n4','nu2','o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2']
        reference = ('https://www.aviso.altimetry.fr/data/products/'
            'auxiliary-products/global-tide-fes.html')
        attrib['tide']['long_name'] = 'Load_Tide'
        model_format = 'FES'
        TYPE = 'z'
        SCALE = 1.0/100.0

    #-- latitude
    attrib['lat']['long_name'] = 'Latitude_of_measurement'
    attrib['lat']['description'] = ('Corresponding_to_the_measurement_'
        'position_at_the_acquisition_time')
    attrib['lat']['units'] = 'Degrees_North'
    #-- longitude
    attrib['lon']['long_name'] = 'Longitude_of_measurement'
    attrib['lon']['description'] = ('Corresponding_to_the_measurement_'
        'position_at_the_acquisition_time')
    attrib['lon']['units'] = 'Degrees_East'
    #-- tides
    attrib['tide']['description'] = ('Tidal_elevation_from_harmonic_constants_'
        'at_the_measurement_position_at_the_acquisition_time')
    attrib['tide']['reference'] = reference
    attrib['tide']['model'] = MODEL
    attrib['tide']['units'] = 'meters'
    #-- days
    attrib['day']['long_name'] = 'Time'
    attrib['day']['description'] = 'Days relative to 1992-01-01 (MJD:48622)'
    attrib['day']['units'] = 'Days'

    #-- extract information from first input file
    #-- acquisition year, month and day
    #-- number of points
    #-- instrument (PRE-OIB ATM or LVIS, OIB ATM or LVIS)
    if OIB in ('ATM','ATM1b'):
        M1,YYMMDD1,HHMMSS1,AX1,SF1 = re.findall(regex[OIB], input_file).pop()
        #-- early date strings omitted century and millenia (e.g. 93 for 1993)
        if (len(YYMMDD1) == 6):
            ypre,MM1,DD1 = YYMMDD1[:2],YYMMDD1[2:4],YYMMDD1[4:]
            if (np.float(ypre) >= 90):
                YY1 = '{0:4.0f}'.format(np.float(ypre) + 1900.0)
            else:
                YY1 = '{0:4.0f}'.format(np.float(ypre) + 2000.0)
        elif (len(YYMMDD1) == 8):
            YY1,MM1,DD1 = YYMMDD1[:4],YYMMDD1[4:6],YYMMDD1[6:]
    elif OIB in ('LVIS','LVGH'):
        M1,RG1,YY1,MMDD1,RLD1,SS1 = re.findall(regex[OIB], input_file).pop()
        MM1,DD1 = MMDD1[:2],MMDD1[2:]

    #-- read data from input_file
    print('{0} -->'.format(input_file)) if VERBOSE else None
    if (OIB == 'ATM'):
        #-- load IceBridge ATM data from input_file
        dinput,file_lines,HEM = read_ATM_icessn_file(input_file,input_subsetter)
    elif (OIB == 'ATM1b'):
        #-- load IceBridge Level-1b ATM data from input_file
        dinput,file_lines,HEM = read_ATM_qfit_file(input_file,input_subsetter)
    elif OIB in ('LVIS','LVGH'):
        #-- load IceBridge LVIS data from input_file
        dinput,file_lines,HEM = read_LVIS_HDF5_file(input_file,input_subsetter)

    #-- convert time from J2000 to days relative to Jan 1, 1992 (48622mjd)
    #-- J2000: seconds since 2000-01-01 12:00:00 UTC
    t = pyTMD.time.convert_delta_time(dinput['time'],
        epoch1=(2000,1,1,12,0,0), epoch2=(1992,1,1,0,0,0),
        scale=1.0/86400.0)

    #-- read tidal constants and interpolate to grid points
    if model_format in ('OTIS','ATLAS'):
        amp,ph,D,c = extract_tidal_constants(dinput['lon'], dinput['lat'],
            grid_file, model_file, EPSG, TYPE=TYPE, METHOD=METHOD,
            GRID=model_format)
        deltat = np.zeros_like(t)
    elif model_format in ('netcdf'):
        amp,ph,D,c = extract_netcdf_constants(dinput['lon'], dinput['lat'],
            model_directory, grid_file, model_files, TYPE=TYPE, METHOD=METHOD,
            SCALE=SCALE)
        deltat = np.zeros_like(t)
    elif (model_format == 'GOT'):
        amp,ph = extract_GOT_constants(dinput['lon'], dinput['lat'],
            model_directory, model_files, METHOD=METHOD, SCALE=SCALE)
        #-- interpolate delta times from calendar dates to tide time
        delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
        deltat = calc_delta_time(delta_file, t)
    elif (model_format == 'FES'):
        amp,ph = extract_FES_constants(dinput['lon'], dinput['lat'],
            model_directory, model_files, TYPE=TYPE, VERSION=MODEL,
            METHOD=METHOD, SCALE=SCALE)
        deltat = np.zeros_like(t)

    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- calculate constituent oscillation
    hc = amp*np.exp(cph)

    #-- output tidal HDF5 file
    #-- form: rg_NASA_model_TIDES_WGS84_fl1yyyymmddjjjjj.H5
    #-- where rg is the hemisphere flag (GR or AN) for the region
    #-- model is the tidal MODEL flag (e.g. CATS0201)
    #-- fl1 and fl2 are the data flags (ATM, LVIS, GLAS)
    #-- yymmddjjjjj is the year, month, day and second of the input file
    #-- output region flags: GR for Greenland and AN for Antarctica
    hem_flag = {'N':'GR','S':'AN'}
    #-- use starting second to distinguish between files for the day
    JJ1 = np.min(dinput['time']) % 86400
    #-- output file format
    args = (hem_flag[HEM],MODEL,OIB,YY1,MM1,DD1,JJ1)
    FILENAME = '{0}_NASA_{1}_TIDES_WGS84_{2}{3}{4}{5}{6:05.0f}.H5'.format(*args)
    #-- print file information
    print('\t{0}'.format(FILENAME)) if VERBOSE else None

    #-- open output HDF5 file
    fid = h5py.File(os.path.join(DIRECTORY,FILENAME), 'w')

    #-- predict tidal elevations at time and infer minor corrections
    fill_value = -9999.0
    tide = np.ma.empty((file_lines),fill_value=fill_value)
    tide.mask = np.any(hc.mask,axis=1)
    tide.data[:] = predict_tide_drift(t, hc, c,
        DELTAT=deltat, CORRECTIONS=model_format)
    minor = infer_minor_corrections(t, hc, c,
        DELTAT=deltat, CORRECTIONS=model_format)
    tide.data[:] += minor.data[:]
    #-- replace invalid values with fill value
    tide.data[tide.mask] = tide.fill_value

    #-- add latitude and longitude to output file
    for key in ['lat','lon']:
        #-- Defining the HDF5 dataset variables for lat/lon
        h5 = fid.create_dataset(key, (file_lines,), data=dinput[key][:],
            dtype=dinput[key].dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in attrib[key].items():
            h5.attrs[att_name] = att_val
        #-- attach dimensions
        h5.dims[0].label = 'RECORD_SIZE'

    #-- output tides to HDF5 dataset
    h5 = fid.create_dataset('tide', (file_lines,), data=tide,
        dtype=tide.dtype, fillvalue=fill_value, compression='gzip')
    #-- add HDF5 variable attributes
    tide_count = np.count_nonzero(tide != fill_value)
    h5.attrs['tide_count'] = tide_count
    h5.attrs['_FillValue'] = fill_value
    for att_name,att_val in attrib['tide'].items():
        h5.attrs[att_name] = att_val
    #-- attach dimensions
    h5.dims[0].label = 'RECORD_SIZE'

    #-- output days to HDF5 dataset
    h5 = fid.create_dataset('day', (file_lines,), data=t,
        dtype=t.dtype, compression='gzip')
    #-- add HDF5 variable attributes
    for att_name,att_val in attrib['day'].items():
        h5.attrs[att_name] = att_val
    #-- attach dimensions
    h5.dims[0].label = 'RECORD_SIZE'

    #-- HDF5 file attributes
    fid.attrs['featureType'] = 'trajectory'
    fid.attrs['title'] = 'Tidal_correction_for_elevation_measurements'
    fid.attrs['summary'] = ('Tidal_correction_computed_at_elevation_'
        'measurements_using_a_tidal_model_driver.')
    fid.attrs['project'] = 'NASA_Operation_IceBridge'
    fid.attrs['processing_level'] = '4'
    fid.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())
    #-- add attributes for input files
    fid.attrs['elevation_file'] = os.path.basename(input_file)
    fid.attrs['tide_model'] = MODEL
    #-- add geospatial and temporal attributes
    fid.attrs['geospatial_lat_min'] = dinput['lat'].min()
    fid.attrs['geospatial_lat_max'] = dinput['lat'].max()
    fid.attrs['geospatial_lon_min'] = dinput['lon'].min()
    fid.attrs['geospatial_lon_max'] = dinput['lon'].max()
    fid.attrs['geospatial_lat_units'] = "degrees_north"
    fid.attrs['geospatial_lon_units'] = "degrees_east"
    fid.attrs['geospatial_ellipsoid'] = "WGS84"
    fid.attrs['time_type'] = 'UTC'
    #-- convert start/end time from days since 1992-01-01 into Julian days
    time_range = np.array([np.min(t),np.max(t)])
    time_julian = 2400000.5 + pyTMD.time.convert_delta_time(time_range,
        epoch1=(1992,1,1,0,0,0), epoch2=(1858,11,17,0,0,0), scale=1.0)
    #-- convert to calendar date with convert_julian.py
    cal = convert_julian(time_julian,ASTYPE=np.int)
    #-- add attributes with measurement date start, end and duration
    args = (cal['hour'][0],cal['minute'][0],cal['second'][0])
    fid.attrs['RangeBeginningTime'] = '{0:02d}:{1:02d}:{2:02d}'.format(*args)
    args = (cal['hour'][-1],cal['minute'][-1],cal['second'][-1])
    fid.attrs['RangeEndingTime'] = '{0:02d}:{1:02d}:{2:02d}'.format(*args)
    args = (cal['year'][0],cal['month'][0],cal['day'][0])
    fid.attrs['RangeBeginningDate'] = '{0:4d}-{1:02d}-{2:02d}'.format(*args)
    args = (cal['year'][-1],cal['month'][-1],cal['day'][-1])
    fid.attrs['RangeEndingDate'] = '{0:4d}-{1:02d}-{2:02d}'.format(*args)
    duration = np.round(time_julian[-1]*86400.0 - time_julian[0]*86400.0)
    fid.attrs['DurationTimeSeconds'] = '{0:0.0f}'.format(duration)
    #-- close the output HDF5 dataset
    fid.close()
    #-- change the permissions level to MODE
    os.chmod(os.path.join(DIRECTORY,FILENAME), MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -T X, --tide=X\t\tTide model to use in correction')
    print(' -I X, --interpolate\tInterpolation method (default spline)')
    print(' -M X, --mode=X\t\tPermission mode of directories and files created')
    print(' -V, --verbose\t\tOutput information about each created file\n')

#-- Main program that calls compute_tides_icebridge_data()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','tide=','interpolate=','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:T:I:VM:', long_options)

    #-- directory with tide data
    tide_dir = os.getcwd()
    #-- tide model to use
    MODEL = 'CATS2008'
    #-- interpolation method
    METHOD = 'spline'
    #-- verbosity settings
    VERBOSE = False
    #-- permissions mode of the local files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            tide_dir = os.path.expanduser(arg)
        elif opt in ("-T","--tide"):
            MODEL = arg
        elif opt in ("-I","--interpolate"):
            METHOD = arg.lower()
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- enter input file from NSIDC as system argument
    if not arglist:
        raise Exception('No System Arguments Listed')

    #-- verify model before running program
    model_list = ['CATS0201','CATS2008','CATS2008_load','TPXO9-atlas','TPXO9.1',
        'TPXO8-atlas','TPXO7.2','TPXO7.2_load','AODTM-5','AOTIM-5',
        'AOTIM-5-2018','GOT4.7','GOT4.7_load','GOT4.8','GOT4.8_load',
        'GOT4.10','GOT4.10_load','FES2014','FES2014_load']
    assert MODEL in model_list, 'Unlisted tide model'

    #-- run for each input file
    for arg in arglist:
        compute_tides_icebridge_data(tide_dir, os.path.expanduser(arg), MODEL,
            METHOD=METHOD, VERBOSE=VERBOSE,MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
