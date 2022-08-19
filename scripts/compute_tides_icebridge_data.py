#!/usr/bin/env python
u"""
compute_tides_icebridge_data.py
Written by Tyler Sutterley (07/2022)
Calculates tidal elevations for correcting Operation IceBridge elevation data

Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
or Finite Element Solution (FES) models provided by AVISO

INPUTS:
    ATM1B, ATM icessn or LVIS file from NSIDC

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -T X, --tide X: Tide model to use in correction
    --atlas-format X: ATLAS tide model format (OTIS, netcdf)
    --gzip, -G: Tide model files are gzip compressed
    --definition-file X: Model definition file for use as correction
    -I X, --interpolate X: Interpolation method
        spline
        linear
        nearest
        bilinear
    -E X, --extrapolate X: Extrapolate with nearest-neighbors
    -c X, --cutoff X: Extrapolation cutoff in kilometers
        set to inf to extrapolate for all points
    --apply-flexure: Apply ice flexure scaling factor to height constituents
        Only valid for models containing flexure fields
    -M X, --mode X: Permission mode of directories and files created
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
    model.py: retrieves tide model parameters for named tide models
    utilities.py: download and management utilities for syncing files
    calc_astrol_longitudes.py: computes the basic astronomical mean longitudes
    calc_delta_time.py: calculates difference between universal and dynamic time
    convert_ll_xy.py: convert lat/lon points to and from projected coordinates
    infer_minor_corrections.py: return corrections for minor constituents
    load_constituent.py: loads parameters for a given tidal constituent
    load_nodal_corrections.py: load the nodal corrections for tidal constituents
    read_tide_model.py: extract tidal harmonic constants from OTIS tide models
    read_netcdf_model.py: extract tidal harmonic constants from netcdf models
    read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
    read_FES_model.py: extract tidal harmonic constants from FES tide models
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates
    predict_tide_drift.py: predict tidal elevations using harmonic constants
    read_ATM1b_QFIT_binary.py: read ATM1b QFIT binary files (NSIDC version 1)

UPDATE HISTORY:
    Updated 07/2022: update imports of ATM1b QFIT functions to released version
        place some imports within try/except statements
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        updated keyword arguments to read tide model programs
        added command line option to apply flexure for applicable models
    Updated 04/2022: include utf-8 encoding in reads to be windows compliant
        use argparse descriptions within sphinx documentation
    Updated 03/2022: using static decorators to define available models
    Updated 02/2022: added Arctic 2km model (Arc2kmTM) to list of models
    Updated 12/2021: added TPXO9-atlas-v5 to list of available tide models
    Updated 10/2021: using python logging for handling verbose output
        using collections to store attributes in order of creation
    Updated 09/2021: refactor to use model class for files and attributes
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 06/2021: added new Gr1km-v2 1km Greenland model from ESR
    Updated 05/2021: added option for extrapolation cutoff in kilometers
        modified import of ATM1b QFIT reader
    Updated 03/2021: added TPXO9-atlas-v4 in binary OTIS format
        simplified netcdf inputs to be similar to binary OTIS read program
        replaced numpy bool/int to prevent deprecation warnings
    Updated 12/2020: added valid data extrapolation with nearest_extrap
         merged time conversion routines into module
    Updated 11/2020: added model constituents from TPXO9-atlas-v3
    Updated 10/2020: using argparse to set command line parameters
    Updated 09/2020: output ocean and load tide as tide_ocean and tide_load
    Updated 08/2020: using builtin time operations.  python3 regular expressions
    Updated 07/2020: added FES2014 and FES2014_load.  use merged delta times
    Updated 06/2020: added version 2 of TPXO9-atlas (TPXO9-atlas-v2)
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
import logging
import argparse
import warnings
import collections
import numpy as np
import pyTMD.time
import pyTMD.model
import pyTMD.utilities
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide_drift import predict_tide_drift
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants
from pyTMD.read_FES_model import extract_FES_constants
# attempt imports
try:
    import h5py
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("h5py not available")
try:
    import ATM1b_QFIT.read_ATM1b_QFIT_binary
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("ATM1b_QFIT not available")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: reading the number of file lines removing commented lines
def file_length(input_file, input_subsetter, HDF5=False, QFIT=False):
    # subset the data to indices if specified
    if input_subsetter:
        file_lines = len(input_subsetter)
    elif HDF5:
        # read the size of an input variable within a HDF5 file
        with h5py.File(input_file,'r') as fileID:
            file_lines, = fileID[HDF5].shape
    elif QFIT:
        # read the size of a QFIT binary file
        file_lines = ATM1b_QFIT.ATM1b_QFIT_shape(input_file)
    else:
        # read the input file, split at lines and remove all commented lines
        with open(input_file, mode='r', encoding='utf8') as f:
            i = [i for i in f.readlines() if re.match(r'^(?!\#|\n)',i)]
        file_lines = len(i)
    # return the number of lines
    return file_lines

## PURPOSE: read the ATM Level-1b data file for variables of interest
def read_ATM_qfit_file(input_file, input_subsetter):
    # regular expression pattern for extracting parameters
    mission_flag = r'(BLATM1B|ILATM1B|ILNSA1B)'
    regex_pattern = r'{0}_(\d+)_(\d+)(.*?).(qi|TXT|h5)'.format(mission_flag)
    # extract mission and other parameters from filename
    MISSION,YYMMDD,HHMMSS,AUX,SFX = re.findall(regex_pattern,input_file).pop()
    # early date strings omitted century and millenia (e.g. 93 for 1993)
    if (len(YYMMDD) == 6):
        ypre,month,day = np.array([YYMMDD[:2],YYMMDD[2:4],YYMMDD[4:]],dtype='i')
        year = (ypre + 1900.0) if (ypre >= 90) else (ypre + 2000.0)
    elif (len(YYMMDD) == 8):
        year,month,day = np.array([YYMMDD[:4],YYMMDD[4:6],YYMMDD[6:]],dtype='i')
    # output python dictionary with variables
    ATM_L1b_input = {}
    # Version 1 of ATM QFIT files (ascii)
    # output text file from qi2txt with proper filename format
    # do not use the shortened output format from qi2txt
    if (SFX == 'TXT'):
        # compile regular expression operator for reading lines
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)
        # read the input file, split at lines and remove all commented lines
        with open(input_file, mode='r', encoding='utf8') as f:
            file_contents = [i for i in f.read().splitlines() if
                re.match(r'^(?!\#|\n)',i)]
        # number of lines of data within file
        file_lines = file_length(input_file,input_subsetter)
        # create output variables with length equal to the number of lines
        ATM_L1b_input['lat'] = np.zeros_like(file_contents,dtype=np.float64)
        ATM_L1b_input['lon'] = np.zeros_like(file_contents,dtype=np.float64)
        ATM_L1b_input['data'] = np.zeros_like(file_contents,dtype=np.float64)
        hour = np.zeros_like(file_contents,dtype=np.float64)
        minute = np.zeros_like(file_contents,dtype=np.float64)
        second = np.zeros_like(file_contents,dtype=np.float64)
        # for each line within the file
        for i,line in enumerate(file_contents):
            # find numerical instances within the line
            line_contents = rx.findall(line)
            ATM_L1b_input['lat'][i] = np.float64(line_contents[1])
            ATM_L1b_input['lon'][i] = np.float64(line_contents[2])
            ATM_L1b_input['data'][i] = np.float64(line_contents[3])
            hour[i] = np.float64(line_contents[-1][:2])
            minute[i] = np.float64(line_contents[-1][2:4])
            second[i] = np.float64(line_contents[-1][4:])
    # Version 1 of ATM QFIT files (binary)
    elif (SFX == 'qi'):
        # read input QFIT data file and subset if specified
        fid,h = ATM1b_QFIT.read_ATM1b_QFIT_binary(input_file)
        # number of lines of data within file
        file_lines = file_length(input_file,input_subsetter,QFIT=True)
        ATM_L1b_input['lat'] = fid['latitude'][:]
        ATM_L1b_input['lon'] = fid['longitude'][:]
        ATM_L1b_input['data'] = fid['elevation'][:]
        time_hhmmss = fid['time_hhmmss'][:]
        # extract hour, minute and second from time_hhmmss
        hour = np.zeros_like(time_hhmmss,dtype=np.float64)
        minute = np.zeros_like(time_hhmmss,dtype=np.float64)
        second = np.zeros_like(time_hhmmss,dtype=np.float64)
        # for each line within the file
        for i,packed_time in enumerate(time_hhmmss):
            # convert to zero-padded string with 3 decimal points
            line_contents = '{0:010.3f}'.format(packed_time)
            hour[i] = np.float64(line_contents[:2])
            minute[i] = np.float64(line_contents[2:4])
            second[i] = np.float64(line_contents[4:])
    # Version 2 of ATM QFIT files (HDF5)
    elif (SFX == 'h5'):
        # Open the HDF5 file for reading
        fileID = h5py.File(os.path.expanduser(input_file), 'r')
        # number of lines of data within file
        file_lines = file_length(input_file,input_subsetter,HDF5='elevation')
        # create output variables with length equal to input elevation
        ATM_L1b_input['lat'] = fileID['latitude'][:]
        ATM_L1b_input['lon'] = fileID['longitude'][:]
        ATM_L1b_input['data'] = fileID['elevation'][:]
        time_hhmmss = fileID['instrument_parameters']['time_hhmmss'][:]
        # extract hour, minute and second from time_hhmmss
        hour = np.zeros_like(time_hhmmss,dtype=np.float64)
        minute = np.zeros_like(time_hhmmss,dtype=np.float64)
        second = np.zeros_like(time_hhmmss,dtype=np.float64)
        # for each line within the file
        for i,packed_time in enumerate(time_hhmmss):
            # convert to zero-padded string with 3 decimal points
            line_contents = '{0:010.3f}'.format(packed_time)
            hour[i] = np.float64(line_contents[:2])
            minute[i] = np.float64(line_contents[2:4])
            second[i] = np.float64(line_contents[4:])
        # close the input HDF5 file
        fileID.close()
    # calculate the number of leap seconds between GPS time (seconds
    # since Jan 6, 1980 00:00:00) and UTC
    gps_seconds = pyTMD.time.convert_calendar_dates(year,month,day,
        hour=hour,minute=minute,second=second,
        epoch=(1980,1,6,0,0,0),scale=86400.0)
    leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)
    # calculation of Julian day taking into account leap seconds
    # converting to J2000 seconds
    ATM_L1b_input['time'] = pyTMD.time.convert_calendar_dates(year,month,day,
        hour=hour,minute=minute,second=second-leap_seconds,
        epoch=(2000,1,1,12,0,0,0),scale=86400.0)
    # subset the data to indices if specified
    if input_subsetter:
        for key,val in ATM_L1b_input.items():
            ATM_L1b_input[key] = val[input_subsetter]
    # hemispheric shot count
    count = {}
    count['N'] = np.count_nonzero(ATM_L1b_input['lat'] >= 0.0)
    count['S'] = np.count_nonzero(ATM_L1b_input['lat'] < 0.0)
    # determine hemisphere with containing shots in file
    HEM, = [key for key, val in count.items() if val]
    # return the output variables
    return ATM_L1b_input,file_lines,HEM

# PURPOSE: read the ATM Level-2 data file for variables of interest
def read_ATM_icessn_file(input_file, input_subsetter):
    # regular expression pattern for extracting parameters
    regex_pattern=r'(BLATM2|ILATM2)_(\d+)_(\d+)_smooth_nadir(.*?)(csv|seg|pt)$'
    # extract mission and other parameters from filename
    MISSION,YYMMDD,HHMMSS,AUX,SFX = re.findall(regex_pattern,input_file).pop()
    # early date strings omitted century and millenia (e.g. 93 for 1993)
    if (len(YYMMDD) == 6):
        ypre,month,day = np.array([YYMMDD[:2],YYMMDD[2:4],YYMMDD[4:]],dtype='i')
        year = (ypre + 1900.0) if (ypre >= 90) else (ypre + 2000.0)
    elif (len(YYMMDD) == 8):
        year,month,day = np.array([YYMMDD[:4],YYMMDD[4:6],YYMMDD[6:]],dtype='i')
    # input file column names for variables of interest with column indices
    # variables not used: (SNslope:4, WEslope:5, npt_used:7, npt_edit:8, d:9)
    file_dtype = {'seconds':0, 'lat':1, 'lon':2, 'data':3, 'RMS':6, 'track':-1}
    # compile regular expression operator for reading lines (extracts numbers)
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)
    # read the input file, split at lines and remove all commented lines
    with open(input_file, mode='r', encoding='utf8') as f:
        file_contents = [i for i in f.read().splitlines() if
            re.match(r'^(?!\#|\n)',i)]
    # number of lines of data within file
    file_lines = file_length(input_file,input_subsetter)
    # output python dictionary with variables
    ATM_L2_input = {}
    # create output variables with length equal to the number of file lines
    for key in file_dtype.keys():
        ATM_L2_input[key] = np.zeros_like(file_contents, dtype=np.float64)
    # for each line within the file
    for line_number,line_entries in enumerate(file_contents):
        # find numerical instances within the line
        line_contents = rx.findall(line_entries)
        # for each variable of interest: save to dinput as float
        for key,val in file_dtype.items():
            ATM_L2_input[key][line_number] = np.float64(line_contents[val])
    # convert shot time (seconds of day) to J2000
    hour = np.floor(ATM_L2_input['seconds']/3600.0)
    minute = np.floor((ATM_L2_input['seconds'] % 3600)/60.0)
    second = ATM_L2_input['seconds'] % 60.0
    # First column in Pre-IceBridge and ICESSN Version 1 files is GPS time
    if (MISSION == 'BLATM2') or (SFX != 'csv'):
        # calculate the number of leap seconds between GPS time (seconds
        # since Jan 6, 1980 00:00:00) and UTC
        gps_seconds = pyTMD.time.convert_calendar_dates(year,month,day,
            hour=hour,minute=minute,second=second,
            epoch=(1980,1,6,0,0,0),scale=86400.0)
        leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)
    else:
        leap_seconds = 0.0
    # calculation of Julian day
    # converting to J2000 seconds
    ATM_L2_input['time'] = pyTMD.time.convert_calendar_dates(year,month,day,
        hour=hour,minute=minute,second=second-leap_seconds,
        epoch=(2000,1,1,12,0,0,0),scale=86400.0)
    # convert RMS from centimeters to meters
    ATM_L2_input['error'] = ATM_L2_input['RMS']/100.0
    # subset the data to indices if specified
    if input_subsetter:
        for key,val in ATM_L2_input.items():
            ATM_L2_input[key] = val[input_subsetter]
    # hemispheric shot count
    count = {}
    count['N'] = np.count_nonzero(ATM_L2_input['lat'] >= 0.0)
    count['S'] = np.count_nonzero(ATM_L2_input['lat'] < 0.0)
    # determine hemisphere with containing shots in file
    HEM, = [key for key, val in count.items() if val]
    # return the output variables
    return ATM_L2_input,file_lines,HEM

# PURPOSE: read the LVIS Level-2 data file for variables of interest
def read_LVIS_HDF5_file(input_file, input_subsetter):
    # LVIS region flags: GL for Greenland and AQ for Antarctica
    lvis_flag = {'GL':'N','AQ':'S'}
    # regular expression pattern for extracting parameters from HDF5 files
    # computed in read_icebridge_lvis.py
    mission_flag = '(BLVIS2|BVLIS2|ILVIS2|ILVGH2)'
    regex_pattern = r'{0}_(.*?)(\d+)_(\d+)_(R\d+)_(\d+).H5'.format(mission_flag)
    # extract mission, region and other parameters from filename
    MISSION,REGION,YY,MMDD,RLD,SS = re.findall(regex_pattern,input_file).pop()
    LDS_VERSION = '2.0.2' if (int(RLD[1:3]) >= 18) else '1.04'
    # input and output python dictionaries with variables
    file_input = {}
    LVIS_L2_input = {}
    fileID = h5py.File(input_file,'r')
    # create output variables with length equal to input shot number
    file_lines = file_length(input_file,input_subsetter,HDF5='Shot_Number')
    # https://lvis.gsfc.nasa.gov/Data/Data_Structure/DataStructure_LDS104.html
    # https://lvis.gsfc.nasa.gov/Data/Data_Structure/DataStructure_LDS202.html
    if (LDS_VERSION == '1.04'):
        # elevation surfaces
        file_input['elev'] = fileID['Elevation_Surfaces/Elevation_Centroid'][:]
        file_input['elev_low'] = fileID['Elevation_Surfaces/Elevation_Low'][:]
        file_input['elev_high'] = fileID['Elevation_Surfaces/Elevation_High'][:]
        # latitude
        file_input['lat'] = fileID['Geolocation/Latitude_Centroid'][:]
        file_input['lat_low'] = fileID['Geolocation/Latitude_Low'][:]
        # longitude
        file_input['lon'] = fileID['Geolocation/Longitude_Centroid'][:]
        file_input['lon_low'] = fileID['Geolocation/Longitude_Low'][:]
    elif (LDS_VERSION == '2.0.2'):
        # elevation surfaces
        file_input['elev_low'] = fileID['Elevation_Surfaces/Elevation_Low'][:]
        file_input['elev_high'] = fileID['Elevation_Surfaces/Elevation_High'][:]
        # heights above lowest detected mode
        file_input['RH50'] = fileID['Waveform/RH50'][:]
        file_input['RH100'] = fileID['Waveform/RH100'][:]
        # calculate centroidal elevation using 50% of waveform energy
        file_input['elev'] = file_input['elev_low'] + file_input['RH50']
        # latitude
        file_input['lat_top'] = fileID['Geolocation/Latitude_Top'][:]
        file_input['lat_low'] = fileID['Geolocation/Latitude_Low'][:]
        # longitude
        file_input['lon_top'] = fileID['Geolocation/Longitude_Top'][:]
        file_input['lon_low'] = fileID['Geolocation/Longitude_Low'][:]
        # linearly interpolate latitude and longitude to RH50
        file_input['lat'] = file_input['lat_low'] + file_input['RH50'] * \
            (file_input['lat_top'] - file_input['lat_low'])/file_input['RH100']
        file_input['lon'] = file_input['lon_low'] + file_input['RH50'] * \
            (file_input['lon_top'] - file_input['lon_low'])/file_input['RH100']
    # J2000 seconds
    LVIS_L2_input['time'] = fileID['Time/J2000'][:]
    # close the input HDF5 file
    fileID.close()
    # output combined variables
    LVIS_L2_input['data'] = np.zeros_like(file_input['elev'],dtype=np.float64)
    LVIS_L2_input['lon'] = np.zeros_like(file_input['elev'],dtype=np.float64)
    LVIS_L2_input['lat'] = np.zeros_like(file_input['elev'],dtype=np.float64)
    LVIS_L2_input['error'] = np.zeros_like(file_input['elev'],dtype=np.float64)
    # find where elev high is equal to elev low
    # see note about using LVIS centroid elevation product
    # http://lvis.gsfc.nasa.gov/OIBDataStructure.html
    ii = np.nonzero(file_input['elev_low'] == file_input['elev_high'])
    jj = np.nonzero(file_input['elev_low'] != file_input['elev_high'])
    # where lowest point of waveform is equal to highest point -->
    # using the elev_low elevation
    LVIS_L2_input['data'][ii] = file_input['elev_low'][ii]
    # for other locations use the centroid elevation
    # as the centroid is a useful product over rough terrain
    # when you are calculating ice volume change
    LVIS_L2_input['data'][jj] = file_input['elev'][jj]
    # latitude and longitude for each case
    # elevation low == elevation high
    LVIS_L2_input['lon'][ii] = file_input['lon_low'][ii]
    LVIS_L2_input['lat'][ii] = file_input['lat_low'][ii]
    # centroid elevations
    LVIS_L2_input['lon'][jj] = file_input['lon'][jj]
    LVIS_L2_input['lat'][jj] = file_input['lat'][jj]
    # estimated uncertainty for both cases
    LVIS_variance_low = (file_input['elev_low'] - file_input['elev'])**2
    LVIS_variance_high = (file_input['elev_high'] - file_input['elev'])**2
    LVIS_L2_input['error']=np.sqrt((LVIS_variance_low + LVIS_variance_high)/2.0)
    # subset the data to indices if specified
    if input_subsetter:
        for key,val in LVIS_L2_input.items():
            LVIS_L2_input[key] = val[input_subsetter]
    # return the output variables
    return LVIS_L2_input,file_lines,lvis_flag[REGION]

# PURPOSE: read Operation IceBridge data from NSIDC
# compute tides at points and times using tidal model driver algorithms
def compute_tides_icebridge_data(tide_dir, arg, TIDE_MODEL,
    ATLAS_FORMAT=None,
    GZIP=True,
    DEFINITION_FILE=None,
    METHOD='spline',
    EXTRAPOLATE=False,
    CUTOFF=None,
    APPLY_FLEXURE=False,
    VERBOSE=False,
    MODE=0o775):

    # create logger for verbosity level
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logger = pyTMD.utilities.build_logger('pytmd',level=loglevel)

    # get parameters for tide model
    if DEFINITION_FILE is not None:
        model = pyTMD.model(tide_dir).from_file(DEFINITION_FILE)
    else:
        model = pyTMD.model(tide_dir, format=ATLAS_FORMAT,
            compressed=GZIP).elevation(TIDE_MODEL)

    # extract file name and subsetter indices lists
    match_object = re.match(r'(.*?)(\[(.*?)\])?$',arg)
    input_file = os.path.expanduser(match_object.group(1))
    # subset input file to indices
    if match_object.group(2):
        # decompress ranges and add to list
        input_subsetter = []
        for i in re.findall(r'((\d+)-(\d+)|(\d+))',match_object.group(3)):
            input_subsetter.append(int(i[3])) if i[3] else \
                input_subsetter.extend(range(int(i[1]),int(i[2])+1))
    else:
        input_subsetter = None

    # output directory for input_file
    DIRECTORY = os.path.dirname(input_file)
    # calculate if input files are from ATM or LVIS (+GH)
    regex = {}
    regex['ATM'] = r'(BLATM2|ILATM2)_(\d+)_(\d+)_smooth_nadir(.*?)(csv|seg|pt)$'
    regex['ATM1b'] = r'(BLATM1b|ILATM1b)_(\d+)_(\d+)(.*?).(qi|TXT|h5)$'
    regex['LVIS'] = r'(BLVIS2|BVLIS2|ILVIS2)_(.*?)(\d+)_(\d+)_(R\d+)_(\d+).H5$'
    regex['LVGH'] = r'(ILVGH2)_(.*?)(\d+)_(\d+)_(R\d+)_(\d+).H5$'
    for key,val in regex.items():
        if re.match(val, os.path.basename(input_file)):
            OIB = key

    # HDF5 file attributes
    attrib = collections.OrderedDict()
    # time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['description'] = ('Time_corresponding_to_the_measurement_'
        'position')
    attrib['time']['units'] = 'Days since 1992-01-01T00:00:00'
    attrib['time']['standard_name'] = 'time'
    attrib['time']['calendar'] = 'standard'
    # latitude
    attrib['lat'] = {}
    attrib['lat']['long_name'] = 'Latitude_of_measurement'
    attrib['lat']['description'] = ('Corresponding_to_the_measurement_'
        'position_at_the_acquisition_time')
    attrib['lat']['units'] = 'Degrees_North'
    # longitude
    attrib['lon'] = {}
    attrib['lon']['long_name'] = 'Longitude_of_measurement'
    attrib['lon']['description'] = ('Corresponding_to_the_measurement_'
        'position_at_the_acquisition_time')
    attrib['lon']['units'] = 'Degrees_East'
    # tides
    attrib[model.variable] = {}
    attrib[model.variable]['description'] = model.description
    attrib[model.variable]['reference'] = model.reference
    attrib[model.variable]['model'] = model.name
    attrib[model.variable]['units'] = 'meters'
    attrib[model.variable]['long_name'] = model.long_name

    # extract information from first input file
    # acquisition year, month and day
    # number of points
    # instrument (PRE-OIB ATM or LVIS, OIB ATM or LVIS)
    if OIB in ('ATM','ATM1b'):
        M1,YYMMDD1,HHMMSS1,AX1,SF1 = re.findall(regex[OIB], input_file).pop()
        # early date strings omitted century and millenia (e.g. 93 for 1993)
        if (len(YYMMDD1) == 6):
            ypre,MM1,DD1 = YYMMDD1[:2],YYMMDD1[2:4],YYMMDD1[4:]
            if (np.float64(ypre) >= 90):
                YY1 = '{0:4.0f}'.format(np.float64(ypre) + 1900.0)
            else:
                YY1 = '{0:4.0f}'.format(np.float64(ypre) + 2000.0)
        elif (len(YYMMDD1) == 8):
            YY1,MM1,DD1 = YYMMDD1[:4],YYMMDD1[4:6],YYMMDD1[6:]
    elif OIB in ('LVIS','LVGH'):
        M1,RG1,YY1,MMDD1,RLD1,SS1 = re.findall(regex[OIB], input_file).pop()
        MM1,DD1 = MMDD1[:2],MMDD1[2:]

    # read data from input_file
    logger.info('{0} -->'.format(input_file))
    if (OIB == 'ATM'):
        # load IceBridge ATM data from input_file
        dinput,file_lines,HEM = read_ATM_icessn_file(input_file,input_subsetter)
    elif (OIB == 'ATM1b'):
        # load IceBridge Level-1b ATM data from input_file
        dinput,file_lines,HEM = read_ATM_qfit_file(input_file,input_subsetter)
    elif OIB in ('LVIS','LVGH'):
        # load IceBridge LVIS data from input_file
        dinput,file_lines,HEM = read_LVIS_HDF5_file(input_file,input_subsetter)

    # convert time from J2000 to days relative to Jan 1, 1992 (48622mjd)
    # J2000: seconds since 2000-01-01 12:00:00 UTC
    t = pyTMD.time.convert_delta_time(dinput['time'],
        epoch1=(2000,1,1,12,0,0), epoch2=(1992,1,1,0,0,0),
        scale=1.0/86400.0)
    # delta time (TT - UT1) file
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])

    # read tidal constants and interpolate to grid points
    if model.format in ('OTIS','ATLAS','ESR'):
        amp,ph,D,c = extract_tidal_constants(dinput['lon'], dinput['lat'],
            model.grid_file, model.model_file, model.projection,
            type=model.type, method=METHOD, extrapolate=EXTRAPOLATE,
            cutoff=CUTOFF, grid=model.format, apply_flexure=APPLY_FLEXURE)
        deltat = np.zeros_like(t)
    elif model.format in ('netcdf'):
        amp,ph,D,c = extract_netcdf_constants(dinput['lon'], dinput['lat'],
            model.grid_file, model.model_file, type=model.type, method=METHOD,
            extrapolate=EXTRAPOLATE, cutoff=CUTOFF, scale=model.scale,
            compressed=model.compressed)
        deltat = np.zeros_like(t)
    elif (model.format == 'GOT'):
        amp,ph,c = extract_GOT_constants(dinput['lon'], dinput['lat'],
            model.model_file, method=METHOD, extrapolate=EXTRAPOLATE,
            cutoff=CUTOFF, scale=model.scale, compressed=model.compressed)
        # interpolate delta times from calendar dates to tide time
        deltat = calc_delta_time(delta_file, t)
    elif (model.format == 'FES'):
        amp,ph = extract_FES_constants(dinput['lon'], dinput['lat'],
            model.model_file, type=model.type, version=model.version,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            scale=model.scale, compressed=model.compressed)
        # available model constituents
        c = model.constituents
        # interpolate delta times from calendar dates to tide time
        deltat = calc_delta_time(delta_file, t)

    # calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    # calculate constituent oscillation
    hc = amp*np.exp(cph)

    # output tidal HDF5 file
    # form: rg_NASA_model_TIDES_WGS84_fl1yyyymmddjjjjj.H5
    # where rg is the hemisphere flag (GR or AN) for the region
    # model is the tidal model name flag (e.g. CATS0201)
    # fl1 and fl2 are the data flags (ATM, LVIS, GLAS)
    # yymmddjjjjj is the year, month, day and second of the input file
    # output region flags: GR for Greenland and AN for Antarctica
    hem_flag = {'N':'GR','S':'AN'}
    # use starting second to distinguish between files for the day
    JJ1 = np.min(dinput['time']) % 86400
    # flexure flag if being applied
    flexure_flag = '_FLEXURE' if APPLY_FLEXURE and model.flexure else ''
    # output file format
    args = (hem_flag[HEM],model.name,flexure_flag,OIB,YY1,MM1,DD1,JJ1)
    FILENAME = '{0}_NASA_{1}{2}_TIDES_WGS84_{3}{4}{5}{6}{7:05.0f}.H5'.format(*args)
    # print file information
    logger.info('\t{0}'.format(FILENAME))

    # open output HDF5 file
    fid = h5py.File(os.path.join(DIRECTORY,FILENAME), 'w')

    # predict tidal elevations at time and infer minor corrections
    fill_value = -9999.0
    tide = np.ma.empty((file_lines),fill_value=fill_value)
    tide.mask = np.any(hc.mask,axis=1)
    tide.data[:] = predict_tide_drift(t, hc, c,
        deltat=deltat, corrections=model.format)
    minor = infer_minor_corrections(t, hc, c,
        deltat=deltat, corrections=model.format)
    tide.data[:] += minor.data[:]
    # replace invalid values with fill value
    tide.data[tide.mask] = tide.fill_value
    # copy tide to output variable
    dinput[model.variable] = tide.copy()

    # output dictionary with HDF5 variables
    h5 = {}
    # add variables to output file
    for key,attributes in attrib.items():
        # Defining the HDF5 dataset variables for lat/lon
        h5[key] = fid.create_dataset(key, (file_lines,),
            data=dinput[key][:], dtype=dinput[key].dtype,
            compression='gzip')
        # add HDF5 variable attributes
        for att_name,att_val in attributes.items():
            h5[key].attrs[att_name] = att_val
        # attach dimensions
        if key not in ('time',):
            for i,dim in enumerate(['time']):
                h5[key].dims[i].label = 'RECORD_SIZE'
                h5[key].dims[i].attach_scale(h5[dim])

    # HDF5 file attributes
    fid.attrs['featureType'] = 'trajectory'
    fid.attrs['title'] = 'Tidal_correction_for_elevation_measurements'
    fid.attrs['summary'] = ('Tidal_correction_computed_at_elevation_'
        'measurements_using_a_tidal_model_driver.')
    fid.attrs['project'] = 'NASA_Operation_IceBridge'
    fid.attrs['processing_level'] = '4'
    fid.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())
    # add attributes for input file
    fid.attrs['elevation_file'] = os.path.basename(input_file)
    fid.attrs['tide_model'] = model.name
    # add geospatial and temporal attributes
    fid.attrs['geospatial_lat_min'] = dinput['lat'].min()
    fid.attrs['geospatial_lat_max'] = dinput['lat'].max()
    fid.attrs['geospatial_lon_min'] = dinput['lon'].min()
    fid.attrs['geospatial_lon_max'] = dinput['lon'].max()
    fid.attrs['geospatial_lat_units'] = "degrees_north"
    fid.attrs['geospatial_lon_units'] = "degrees_east"
    fid.attrs['geospatial_ellipsoid'] = "WGS84"
    fid.attrs['time_type'] = 'UTC'
    # convert start/end time from days since 1992-01-01 into Julian days
    time_range = np.array([np.min(t),np.max(t)])
    time_julian = 2400000.5 + pyTMD.time.convert_delta_time(time_range,
        epoch1=(1992,1,1,0,0,0), epoch2=(1858,11,17,0,0,0), scale=1.0)
    # convert to calendar date
    cal = pyTMD.time.convert_julian(time_julian,astype=int)
    # add attributes with measurement date start, end and duration
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
    # close the output HDF5 dataset
    fid.close()
    # change the permissions level to MODE
    os.chmod(os.path.join(DIRECTORY,FILENAME), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates tidal elevations for correcting Operation
            IceBridge elevation data
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    # command line parameters
    group = parser.add_mutually_exclusive_group(required=True)
    # input operation icebridge files
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Input Operation IceBridge file to run')
    # directory with tide data
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # tide model to use
    choices = sorted(pyTMD.model.ocean_elevation() + pyTMD.model.load_elevation())
    group.add_argument('--tide','-T',
        metavar='TIDE', type=str,
        choices=choices,
        help='Tide model to use in correction')
    parser.add_argument('--atlas-format',
        type=str, choices=('OTIS','netcdf'), default='netcdf',
        help='ATLAS tide model format')
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='Tide model files are gzip compressed')
    # tide model definition file to set an undefined model
    group.add_argument('--definition-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Tide model definition file for use as correction')
    # interpolation method
    parser.add_argument('--interpolate','-I',
        metavar='METHOD', type=str, default='spline',
        choices=('spline','linear','nearest','bilinear'),
        help='Spatial interpolation method')
    # extrapolate with nearest-neighbors
    parser.add_argument('--extrapolate','-E',
        default=False, action='store_true',
        help='Extrapolate with nearest-neighbors')
    # extrapolation cutoff in kilometers
    # set to inf to extrapolate over all points
    parser.add_argument('--cutoff','-c',
        type=np.float64, default=10.0,
        help='Extrapolation cutoff in kilometers')
    # apply flexure scaling factors to height constituents
    parser.add_argument('--apply-flexure',
        default=False, action='store_true',
        help='Apply ice flexure scaling factor to height constituents')
    # verbosity settings
    # verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Output information about each created file')
    # permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files created')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run for each input Operation IceBridge file
    for arg in args.infile:
        compute_tides_icebridge_data(args.directory, arg,
            TIDE_MODEL=args.tide,
            ATLAS_FORMAT=args.atlas_format,
            GZIP=args.gzip,
            DEFINITION_FILE=args.definition_file,
            METHOD=args.interpolate,
            EXTRAPOLATE=args.extrapolate,
            CUTOFF=args.cutoff,
            APPLY_FLEXURE=args.apply_flexure,
            VERBOSE=args.verbose,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
