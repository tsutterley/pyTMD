#!/usr/bin/env python
u"""
compute_LPT_icebridge_data.py
Written by Tyler Sutterley (08/2020)
Calculates load pole tide displacements for correcting Operation IceBridge
    elevation data following IERS Convention (2010) guidelines
    http://maia.usno.navy.mil/conventions/2010officialinfo.php
    http://maia.usno.navy.mil/conventions/chapter7.php

INPUTS:
    ATM1B, ATM icessn or LVIS file from NSIDC

COMMAND LINE OPTIONS:
    -D X, --directory=X: Working data directory
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

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    utilities: download and management utilities for syncing files
    convert_julian.py: returns the calendar date and time given a Julian date
    convert_calendar_decimal.py: converts from calendar dates into decimal years
    iers_mean_pole.py: provides the angular coordinates of IERS Mean Pole
    read_iers_EOP.py: read daily earth orientation parameters from IERS
    read_ATM1b_QFIT_binary.py: read ATM1b QFIT binary files (NSIDC version 1)

UPDATE HISTORY:
    Updated 08/2020: using builtin time operations
    Updated 03/2020: use read_ATM1b_QFIT_binary from repository
    Updated 02/2019: using range for python3 compatibility
    Updated 10/2018: updated GPS time calculation for calculating leap seconds
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
import scipy.interpolate
import pyTMD.time
from pyTMD.convert_julian import convert_julian
from pyTMD.iers_mean_pole import iers_mean_pole
from pyTMD.read_iers_EOP import read_iers_EOP
from read_ATM1b_QFIT_binary.read_ATM1b_QFIT_binary import read_ATM1b_QFIT_binary

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
        file_lines = read_ATM1b_QFIT_binary(input_file)
    else:
        #-- read the input file, split at lines and remove all commented lines
        with open(input_file,'r') as f:
            i = [i for i in f.read().splitlines() if re.match('^(?!#)',i)]
        file_lines = len(i)
    #-- return the number of lines
    return file_lines

#-- PURPOSE: read the ATM Level-1b data file for variables of interest
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
#-- compute load pole tide radial displacements at data points and times
def compute_LPT_icebridge_data(tide_dir, arg, VERBOSE=False, MODE=0o775):

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
    attrib = {}
    #-- latitude
    attrib['lat'] = {}
    attrib['lat']['long_name'] = 'Latitude_of_measurement'
    attrib['lat']['description'] = ('Corresponding_to_the_measurement_'
        'position_at_the_acquisition_time')
    attrib['lat']['units'] = 'Degrees_North'
    #-- longitude
    attrib['lon'] = {}
    attrib['lon']['long_name'] = 'Longitude_of_measurement'
    attrib['lon']['description'] = ('Corresponding_to_the_measurement_'
        'position_at_the_acquisition_time')
    attrib['lon']['units'] = 'Degrees_East'
    #-- load pole tides
    attrib['tide_pole'] = {}
    attrib['tide_pole']['long_name'] = 'Solid_Earth_Pole_Tide'
    attrib['tide_pole']['description'] = ('Solid_Earth_pole_tide_radial_'
        'displacements_at_the_measurement_position_at_the_acquisition_'
        'time_due_to_polar_motion')
    attrib['tide_pole']['reference'] = ('ftp://tai.bipm.org/iers/conv2010/'
        'chapter7/opoleloadcoefcmcor.txt.gz')
    attrib['tide_pole']['units'] = 'meters'
    #-- Modified Julian Days
    attrib['MJD'] = {}
    attrib['MJD']['long_name'] = 'Time'
    attrib['MJD']['description'] = 'Modified Julian Days'
    attrib['MJD']['units'] = 'Days'

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

    #-- extract lat/lon
    lon = dinput['lon'][:]
    lat = dinput['lat'][:]
    #-- convert time from UTC time of day to modified julian days (MJD)
    #-- J2000: seconds since 2000-01-01 12:00:00 UTC
    t = dinput['time'][:]/86400.0 + 51544.5
    #-- convert from MJD to calendar dates
    YY,MM,DD,HH,MN,SS = convert_julian(t + 2400000.5,FORMAT='tuple')
    #-- convert calendar dates into year decimal
    tdec = convert_calendar_decimal(YY,MM,DAY=DD,HOUR=HH,MINUTE=MN,SECOND=SS)
    #-- elevation
    h1 = dinput['data'][:]

    #-- degrees to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    #-- earth and physical parameters (IERS and WGS84)
    G = 6.67428e-11#-- universal constant of gravitation [m^3/(kg*s^2)]
    GM = 3.98004418e14#-- geocentric gravitational constant [m^3/s^2]
    ge = 9.7803278#-- mean equatorial gravity [m/s^2]
    a_axis = 6378136.6#-- semimajor axis of the WGS84 ellipsoid [m]
    flat = 1.0/298.257223563#-- flattening of the WGS84 ellipsoid
    b_axis = (1.0 -flat)*a_axis#-- semiminor axis of the WGS84 ellipsoid [m]
    omega = 7.292115e-5#-- mean rotation rate of the Earth [radians/s]
    #-- tidal love number appropriate for the load tide
    hb2 = 0.6207
    #-- Linear eccentricity, first and second numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    ecc2 = lin_ecc/b_axis
    #-- m parameter [omega^2*a^2*b/(GM)]. p. 70, Eqn.(2-137)
    m = omega**2*((1 -flat)*a_axis**3)/GM
    #-- flattening components
    f_2 = -flat + (5.0/2.0)*m + (1.0/2.0)*flat**2.0 - (26.0/7.0)*flat*m + \
        (15.0/4.0)*m**2.0
    f_4 = -(1.0/2.0)*flat**2.0 + (5.0/2.0)*flat*m

    #-- convert from geodetic latitude to geocentric latitude
    #-- geodetic latitude in radians
    latitude_geodetic_rad = lat*dtr
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    #-- calculate X, Y and Z from geodetic latitude and longitude
    X = (N + h1) * np.cos(latitude_geodetic_rad) * np.cos(lon*dtr)
    Y = (N + h1) * np.cos(latitude_geodetic_rad) * np.sin(lon*dtr)
    Z = (N * (1.0 - ecc1**2.0) + h1) * np.sin(latitude_geodetic_rad)
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    #-- calculate geocentric latitude and convert to degrees
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/dtr
    #-- colatitude and longitude in radians
    theta = dtr*(90.0 - latitude_geocentric)
    phi = lon*dtr

    #-- compute normal gravity at spatial location and elevation of points.
    #-- normal gravity at the equator. p. 79, Eqn.(2-186)
    gamma_a = (GM/(a_axis*b_axis)) * (1.0-(3.0/2.0)*m - (3.0/14.0)*ecc2**2.0*m)
    #-- Normal gravity. p. 80, Eqn.(2-199)
    gamma_0 = gamma_a*(1.0 + f_2*np.cos(theta)**2.0 +
        f_4*np.sin(np.pi*latitude_geocentric/180.0)**4.0)
    #-- Normal gravity at height h. p. 82, Eqn.(2-215)
    gamma_h = gamma_0*(1.0-(2.0/a_axis)*(1.0+flat+m-2.0*flat*np.cos(theta)**2.0)
        *h1 + (3.0/a_axis**2.0)*h1**2.0)

    #-- pole tide files (mean and daily)
    # mean_pole_file = os.path.join(tide_dir,'mean-pole.tab')
    mean_pole_file = os.path.join(tide_dir,'mean_pole_2017-10-23.tab')
    pole_tide_file = os.path.join(tide_dir,'finals_all_2017-09-01.tab')
    #-- read IERS daily polar motion values
    EOP = read_iers_EOP(pole_tide_file)
    #-- create cubic spline interpolations of daily polar motion values
    xSPL = scipy.interpolate.UnivariateSpline(EOP['MJD'],EOP['x'],k=3,s=0)
    ySPL = scipy.interpolate.UnivariateSpline(EOP['MJD'],EOP['y'],k=3,s=0)
    #-- bad value
    fill_value = -9999.0

    #-- output load pole tide HDF5 file
    #-- form: rg_NASA_LOAD_POLE_TIDE_WGS84_fl1yyyymmddjjjjj.H5
    #-- where rg is the hemisphere flag (GR or AN) for the region
    #-- fl1 and fl2 are the data flags (ATM, LVIS, GLAS)
    #-- yymmddjjjjj is the year, month, day and second of the input file
    #-- output region flags: GR for Greenland and AN for Antarctica
    hem_flag = {'N':'GR','S':'AN'}
    #-- use starting second to distinguish between files for the day
    JJ1 = np.min(dinput['time']) % 86400
    #-- output file format
    args = (hem_flag[HEM],'LOAD_POLE_TIDE',OIB,YY1,MM1,DD1,JJ1)
    FILENAME = '{0}_NASA_{1}_WGS84_{2}{3}{4}{5}{6:05.0f}.H5'.format(*args)
    #-- print file information
    print('\t{0}'.format(FILENAME)) if VERBOSE else None

    #-- open output HDF5 file
    fid = h5py.File(os.path.join(DIRECTORY,FILENAME), 'w')

    #-- calculate angular coordinates of mean pole at time tdec
    mpx,mpy,fl = iers_mean_pole(mean_pole_file,tdec,'2015')
    #-- interpolate daily polar motion values to time using cubic splines
    px = xSPL(t)
    py = ySPL(t)
    #-- calculate differentials from mean pole positions
    mx = px - mpx
    my = -(py - mpy)
    #-- calculate radial displacement at time
    dfactor = -hb2*atr*(omega**2*rr**2)/(2.0*gamma_h)
    Sr = dfactor*np.sin(2.0*theta)*(mx*np.cos(phi) + my*np.sin(phi))

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
    h5 = fid.create_dataset('tide_pole', (file_lines,), data=Sr,
        dtype=Sr.dtype, compression='gzip')
    #-- add HDF5 variable attributes
    h5.attrs['_FillValue'] = fill_value
    for att_name,att_val in attrib['tide_pole'].items():
        h5.attrs[att_name] = att_val
    #-- attach dimensions
    h5.dims[0].label = 'RECORD_SIZE'

    #-- output days to HDF5 dataset
    h5 = fid.create_dataset('MJD', (file_lines,), data=t, dtype=t.dtype,
        compression='gzip')
    #-- add HDF5 variable attributes
    for att_name,att_val in attrib['MJD'].items():
        h5.attrs[att_name] = att_val
    #-- attach dimensions
    h5.dims[0].label = 'RECORD_SIZE'

    #-- HDF5 file attributes
    fid.attrs['featureType'] = 'trajectory'
    fid.attrs['title'] = 'Load_Pole_Tide_correction_for_elevation_measurements'
    fid.attrs['summary'] = ('Solid_Earth_pole_tide_radial_displacements_'
        'computed_at_elevation_measurements.')
    fid.attrs['project'] = 'NASA_Operation_IceBridge'
    fid.attrs['processing_level'] = '4'
    fid.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())
    #-- add attributes for input files
    fid.attrs['elevation_file'] = os.path.basename(input_file)
    #-- add geospatial and temporal attributes
    fid.attrs['geospatial_lat_min'] = dinput['lat'].min()
    fid.attrs['geospatial_lat_max'] = dinput['lat'].max()
    fid.attrs['geospatial_lon_min'] = dinput['lon'].min()
    fid.attrs['geospatial_lon_max'] = dinput['lon'].max()
    fid.attrs['geospatial_lat_units'] = "degrees_north"
    fid.attrs['geospatial_lon_units'] = "degrees_east"
    fid.attrs['geospatial_ellipsoid'] = "WGS84"
    fid.attrs['time_type'] = 'UTC'

    #-- convert start/end time from MJD into Julian days
    JD_start = np.min(t) + 2400000.5
    JD_end = np.max(t) + 2400000.5
    #-- convert to calendar date with convert_julian.py
    cal = convert_julian(np.array([JD_start,JD_end]),ASTYPE=np.int)
    #-- add attributes with measurement date start, end and duration
    args = (cal['hour'][0],cal['minute'][0],cal['second'][0])
    fid.attrs['RangeBeginningTime'] = '{0:02d}:{1:02d}:{2:02d}'.format(*args)
    args = (cal['hour'][-1],cal['minute'][-1],cal['second'][-1])
    fid.attrs['RangeEndingTime'] = '{0:02d}:{1:02d}:{2:02d}'.format(*args)
    args = (cal['year'][0],cal['month'][0],cal['day'][0])
    fid.attrs['RangeBeginningDate'] = '{0:4d}-{1:02d}-{2:02d}'.format(*args)
    args = (cal['year'][-1],cal['month'][-1],cal['day'][-1])
    fid.attrs['RangeEndingDate'] = '{0:4d}-{1:02d}-{2:02d}'.format(*args)
    duration = np.round(JD_end*86400.0 - JD_start*86400.0)
    fid.attrs['DurationTimeSeconds'] ='{0:0.0f}'.format(duration)
    #-- close the output HDF5 dataset
    fid.close()
    #-- change the permissions level to MODE
    os.chmod(os.path.join(DIRECTORY,FILENAME), MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -M X, --mode=X\t\tPermission mode of directories and files created')
    print(' -V, --verbose\t\tOutput information about each created file\n')

#-- Main program that calls compute_LPT_icebridge_data()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:VM:', long_options)

    #-- directory with tide data
    tide_dir = os.getcwd()
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
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- enter input file from NSIDC as system argument
    if not arglist:
        raise Exception('No System Arguments Listed')

    #-- run for each input file
    for arg in arglist:
        compute_LPT_icebridge_data(tide_dir, arg, VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
