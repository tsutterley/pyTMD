#!/usr/bin/env python
u"""
compute_tides_ICESat_GLA12.py
Written by Tyler Sutterley (07/2022)
Calculates tidal elevations for correcting ICESat/GLAS L2 GLA12
    Antarctic and Greenland Ice Sheet elevation data

Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
or Finite Element Solution (FES) models provided by AVISO

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
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    model.py: retrieves tide model parameters for named tide models
    spatial: utilities for reading, writing and operating on spatial data
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

UPDATE HISTORY:
    Updated 07/2022: place some imports within try/except statements
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        updated keyword arguments to read tide model programs
        added command line option to apply flexure for applicable models
    Updated 04/2022: use argparse descriptions within documentation
    Updated 03/2022: using static decorators to define available models
    Updated 02/2022: save ICESat campaign attribute to output file
        added Arctic 2km model (Arc2kmTM) to list of available tide models
    Updated 12/2021: added TPXO9-atlas-v5 to list of available tide models
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: refactor to use model class for files and attributes
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 06/2021: added new Gr1km-v2 1km Greenland model from ESR
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 04/2021: can use a generically named GLA12 file as input
    Updated 03/2021: added TPXO9-atlas-v4 in binary OTIS format
        simplified netcdf inputs to be similar to binary OTIS read program
    Updated 12/2020: updated for public release
        H5py deprecation warning change to use make_scale and not create_scale
        added valid data extrapolation with nearest_extrap
    Updated 11/2020: added model constituents from TPXO9-atlas-v3
    Updated 10/2020: using argparse to set command line parameters
    Updated 08/2020: using builtin time operations.  python3 regular expressions
    Updated 07/2020: added FES2014 and FES2014_load.  use merged delta times
    Updated 06/2020: added version 2 of TPXO9-atlas (TPXO9-atlas-v2)
    Updated 02/2020: changed CATS2008 grid to match version on U.S. Antarctic
        Program Data Center http://www.usap-dc.org/view/dataset/601235
    Updated 11/2019: calculate minor constituents as separate variable
        compute tide values at all segments and then mask to valid
        added AOTIM-5-2018 tide model (2018 update to 2004 model)
    Updated 10/2019: external read functions.  adjust regex for processed files
        changing Y/N flags to True/False
    Updated 09/2019: using date functions paralleling public repository
        add option for TPXO9-atlas.  add OTIS netcdf tide option
    Written 12/2018
"""
from __future__ import print_function

import sys
import os
import re
import logging
import argparse
import warnings
import numpy as np
import pyTMD.time
import pyTMD.model
import pyTMD.spatial
import pyTMD.utilities
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants
from pyTMD.read_FES_model import extract_FES_constants
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide_drift import predict_tide_drift
# attempt imports
try:
    import h5py
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("h5py not available")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: read ICESat ice sheet HDF5 elevation data (GLAH12) from NSIDC
# compute tides at points and times using tidal model driver algorithms
def compute_tides_ICESat(tide_dir, INPUT_FILE,
    TIDE_MODEL=None,
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

    # get directory from INPUT_FILE
    logger.info('{0} -->'.format(INPUT_FILE))
    DIRECTORY = os.path.dirname(INPUT_FILE)
    # flexure flag if being applied
    flexure_flag = '_FLEXURE' if APPLY_FLEXURE and model.flexure else ''
    # compile regular expression operator for extracting information from file
    rx = re.compile((r'GLAH(\d{2})_(\d{3})_(\d{1})(\d{1})(\d{2})_(\d{3})_'
        r'(\d{4})_(\d{1})_(\d{2})_(\d{4})\.H5'), re.VERBOSE)
    # extract parameters from ICESat/GLAS HDF5 file name
    # PRD:  Product number (01, 05, 06, 12, 13, 14, or 15)
    # RL:  Release number for process that created the product = 634
    # RGTP:  Repeat ground-track phase (1=8-day, 2=91-day, 3=transfer orbit)
    # ORB:   Reference orbit number (starts at 1 and increments each time a
    #           new reference orbit ground track file is obtained.)
    # INST:  Instance number (increments every time the satellite enters a
    #           different reference orbit)
    # CYCL:   Cycle of reference orbit for this phase
    # TRK: Track within reference orbit
    # SEG:   Segment of orbit
    # GRAN:  Granule version number
    # TYPE:  File type
    try:
        PRD,RL,RGTP,ORB,INST,CYCL,TRK,SEG,GRAN,TYPE = rx.findall(INPUT_FILE).pop()
    except:
        # output tide HDF5 file (generic)
        fileBasename,fileExtension = os.path.splitext(INPUT_FILE)
        args = (fileBasename,model.name,flexure_flag,fileExtension)
        OUTPUT_FILE = '{0}_{1}{2}_TIDES{3}'.format(*args)
    else:
        # output tide HDF5 file for NSIDC granules
        args = (PRD,RL,model.name,flexure_flag,RGTP,ORB,INST,CYCL,TRK,SEG,GRAN,TYPE)
        file_format = 'GLAH{0}_{1}_{2}{3}_TIDES_{4}{5}{6}_{7}_{8}_{9}_{10}_{11}.h5'
        OUTPUT_FILE = file_format.format(*args)

    # read GLAH12 HDF5 file
    fileID = h5py.File(INPUT_FILE,'r')
    n_40HZ, = fileID['Data_40HZ']['Time']['i_rec_ndx'].shape
    # get variables and attributes
    rec_ndx_40HZ = fileID['Data_40HZ']['Time']['i_rec_ndx'][:].copy()
    # seconds since 2000-01-01 12:00:00 UTC (J2000)
    DS_UTCTime_40HZ = fileID['Data_40HZ']['DS_UTCTime_40'][:].copy()
    # Latitude (degrees North)
    lat_TPX = fileID['Data_40HZ']['Geolocation']['d_lat'][:].copy()
    # Longitude (degrees East)
    lon_40HZ = fileID['Data_40HZ']['Geolocation']['d_lon'][:].copy()
    # Elevation (height above TOPEX/Poseidon ellipsoid in meters)
    elev_TPX = fileID['Data_40HZ']['Elevation_Surfaces']['d_elev'][:].copy()
    fv = fileID['Data_40HZ']['Elevation_Surfaces']['d_elev'].attrs['_FillValue']

    # semimajor axis (a) and flattening (f) for TP and WGS84 ellipsoids
    atop,ftop = (6378136.3,1.0/298.257)
    awgs,fwgs = (6378137.0,1.0/298.257223563)
    # convert from Topex/Poseidon to WGS84 Ellipsoids
    lat_40HZ,elev_40HZ = pyTMD.spatial.convert_ellipsoid(lat_TPX, elev_TPX,
        atop, ftop, awgs, fwgs, eps=1e-12, itmax=10)

    # convert time from J2000 to days relative to Jan 1, 1992 (48622mjd)
    # J2000: seconds since 2000-01-01 12:00:00 UTC
    tide_time = pyTMD.time.convert_delta_time(DS_UTCTime_40HZ,
        epoch1=(2000,1,1,12,0,0), epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
    # delta time (TT - UT1) file
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
    # read tidal constants and interpolate to grid points
    if model.format in ('OTIS','ATLAS','ESR'):
        amp,ph,D,c = extract_tidal_constants(lon_40HZ, lat_40HZ,
            model.grid_file, model.model_file, model.projection,
            type=model.type, method=METHOD, extrapolate=EXTRAPOLATE,
            cutoff=CUTOFF, grid=model.format, apply_flexure=APPLY_FLEXURE)
        deltat = np.zeros_like(tide_time)
    elif (model.format == 'netcdf'):
        amp,ph,D,c = extract_netcdf_constants(lon_40HZ, lat_40HZ,
            model.grid_file, model.model_file, type=model.type,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            scale=model.scale, compressed=model.compressed)
        deltat = np.zeros_like(tide_time)
    elif (model.format == 'GOT'):
        amp,ph,c = extract_GOT_constants(lon_40HZ, lat_40HZ,
            model.model_file, method=METHOD, extrapolate=EXTRAPOLATE,
            cutoff=CUTOFF, scale=model.scale, compressed=model.compressed)
        # interpolate delta times from calendar dates to tide time
        deltat = calc_delta_time(delta_file, tide_time)
    elif (model.format == 'FES'):
        amp,ph = extract_FES_constants(lon_40HZ, lat_40HZ,
            model.model_file, type=model.type, version=model.version,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            scale=model.scale, compressed=model.compressed)
        # available model constituents
        c = model.constituents
        # interpolate delta times from calendar dates to tide time
        deltat = calc_delta_time(delta_file, tide_time)

    # calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    # calculate constituent oscillation
    hc = amp*np.exp(cph)

    # predict tidal elevations at time and infer minor corrections
    tide = np.ma.empty((n_40HZ),fill_value=fv)
    tide.mask = np.any(hc.mask,axis=1)
    tide.data[:] = predict_tide_drift(tide_time, hc, c,
        deltat=deltat, corrections=model.format)
    minor = infer_minor_corrections(tide_time, hc, c,
        deltat=deltat, corrections=model.format)
    tide.data[:] += minor.data[:]
    # replace masked and nan values with fill value
    invalid, = np.nonzero(np.isnan(tide.data) | tide.mask)
    tide.data[invalid] = tide.fill_value
    tide.mask[invalid] = True

    # copy variables for outputting to HDF5 file
    IS_gla12_tide = dict(Data_40HZ={})
    IS_gla12_fill = dict(Data_40HZ={})
    IS_gla12_tide_attrs = dict(Data_40HZ={})

    # copy global file attributes of interest
    global_attribute_list = ['featureType','title','comment','summary','license',
        'references','AccessConstraints','CitationforExternalPublication',
        'contributor_role','contributor_name','creator_name','creator_email',
        'publisher_name','publisher_email','publisher_url','platform','instrument',
        'processing_level','date_created','spatial_coverage_type','history',
        'keywords','keywords_vocabulary','naming_authority','project','time_type',
        'date_type','time_coverage_start','time_coverage_end',
        'time_coverage_duration','source','HDFVersion','identifier_product_type',
        'identifier_product_format_version','Conventions','institution',
        'ReprocessingPlanned','ReprocessingActual','LocalGranuleID',
        'ProductionDateTime','LocalVersionID','PGEVersion','OrbitNumber',
        'StartOrbitNumber','StopOrbitNumber','EquatorCrossingLongitude',
        'EquatorCrossingTime','EquatorCrossingDate','ShortName','VersionID',
        'InputPointer','RangeBeginningTime','RangeEndingTime','RangeBeginningDate',
        'RangeEndingDate','PercentGroundHit','OrbitQuality','Cycle','Track',
        'Instrument_State','Timing_Bias','ReferenceOrbit','SP_ICE_PATH_NO',
        'SP_ICE_GLAS_StartBlock','SP_ICE_GLAS_EndBlock','Instance','Range_Bias',
        'Instrument_State_Date','Instrument_State_Time','Range_Bias_Date',
        'Range_Bias_Time','Timing_Bias_Date','Timing_Bias_Time',
        'identifier_product_doi','identifier_file_uuid',
        'identifier_product_doi_authority']
    for att in global_attribute_list:
        IS_gla12_tide_attrs[att] = fileID.attrs[att]
    # copy ICESat campaign name from ancillary data
    IS_gla12_tide_attrs['Campaign'] = fileID['ANCILLARY_DATA'].attrs['Campaign']

    # add attributes for input GLA12 file
    IS_gla12_tide_attrs['input_files'] = os.path.basename(INPUT_FILE)
    # update geospatial ranges for ellipsoid
    IS_gla12_tide_attrs['geospatial_lat_min'] = np.min(lat_40HZ)
    IS_gla12_tide_attrs['geospatial_lat_max'] = np.max(lat_40HZ)
    IS_gla12_tide_attrs['geospatial_lon_min'] = np.min(lon_40HZ)
    IS_gla12_tide_attrs['geospatial_lon_max'] = np.max(lon_40HZ)
    IS_gla12_tide_attrs['geospatial_lat_units'] = "degrees_north"
    IS_gla12_tide_attrs['geospatial_lon_units'] = "degrees_east"
    IS_gla12_tide_attrs['geospatial_ellipsoid'] = "WGS84"

    # copy 40Hz group attributes
    for att_name,att_val in fileID['Data_40HZ'].attrs.items():
        IS_gla12_tide_attrs['Data_40HZ'][att_name] = att_val
    # copy attributes for time, geolocation and geophysical groups
    for var in ['Time','Geolocation','Geophysical']:
        IS_gla12_tide['Data_40HZ'][var] = {}
        IS_gla12_fill['Data_40HZ'][var] = {}
        IS_gla12_tide_attrs['Data_40HZ'][var] = {}
        for att_name,att_val in fileID['Data_40HZ'][var].attrs.items():
            IS_gla12_tide_attrs['Data_40HZ'][var][att_name] = att_val

    # J2000 time
    IS_gla12_tide['Data_40HZ']['DS_UTCTime_40'] = DS_UTCTime_40HZ
    IS_gla12_fill['Data_40HZ']['DS_UTCTime_40'] = None
    IS_gla12_tide_attrs['Data_40HZ']['DS_UTCTime_40'] = {}
    for att_name,att_val in fileID['Data_40HZ']['DS_UTCTime_40'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['DS_UTCTime_40'][att_name] = att_val
    # record
    IS_gla12_tide['Data_40HZ']['Time']['i_rec_ndx'] = rec_ndx_40HZ
    IS_gla12_fill['Data_40HZ']['Time']['i_rec_ndx'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Time']['i_rec_ndx'] = {}
    IS_gla12_tide_attrs['Data_40HZ']['Time']['i_rec_ndx']['coordinates'] = \
        "../DS_UTCTime_40"
    for att_name,att_val in fileID['Data_40HZ']['Time']['i_rec_ndx'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Time']['i_rec_ndx'][att_name] = att_val
    # latitude
    IS_gla12_tide['Data_40HZ']['Geolocation']['d_lat'] = lat_40HZ
    IS_gla12_fill['Data_40HZ']['Geolocation']['d_lat'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lat'] = {}
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lat']['coordinates'] = \
        "../DS_UTCTime_40"
    for att_name,att_val in fileID['Data_40HZ']['Geolocation']['d_lat'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lat'][att_name] = att_val
    # longitude
    IS_gla12_tide['Data_40HZ']['Geolocation']['d_lon'] = lon_40HZ
    IS_gla12_fill['Data_40HZ']['Geolocation']['d_lon'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lon'] = {}
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lon']['coordinates'] = \
        "../DS_UTCTime_40"
    for att_name,att_val in fileID['Data_40HZ']['Geolocation']['d_lon'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lon'][att_name] = att_val

    # geophysical variables
    # computed tide
    IS_gla12_tide['Data_40HZ']['Geophysical'][model.gla12] = tide
    IS_gla12_fill['Data_40HZ']['Geophysical'][model.gla12] = tide.fill_value
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][model.gla12] = {}
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][model.gla12]['units'] = "meters"
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][model.gla12]['long_name'] = model.long_name
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][model.gla12]['description'] = model.description
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][model.gla12]['source'] = model.name
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][model.gla12]['reference'] = model.reference
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][model.gla12]['coordinates'] = \
        "../DS_UTCTime_40"

    # close the input HDF5 file
    fileID.close()

    # print file information
    logger.info('\t{0}'.format(OUTPUT_FILE))
    HDF5_GLA12_tide_write(IS_gla12_tide, IS_gla12_tide_attrs,
        FILENAME=os.path.join(DIRECTORY,OUTPUT_FILE),
        FILL_VALUE=IS_gla12_fill, CLOBBER=True)
    # change the permissions mode
    os.chmod(os.path.join(DIRECTORY,OUTPUT_FILE), MODE)

# PURPOSE: outputting the tide values for ICESat data to HDF5
def HDF5_GLA12_tide_write(IS_gla12_tide, IS_gla12_attrs,
    FILENAME='', FILL_VALUE=None, CLOBBER=False):
    # setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    # open output HDF5 file
    fileID = h5py.File(os.path.expanduser(FILENAME), clobber)
    # create 40HZ HDF5 records
    h5 = dict(Data_40HZ={})

    # add HDF5 file attributes
    attrs = {a:v for a,v in IS_gla12_attrs.items() if not isinstance(v,dict)}
    for att_name,att_val in attrs.items():
       fileID.attrs[att_name] = att_val

    # create Data_40HZ group
    fileID.create_group('Data_40HZ')
    # add HDF5 40HZ group attributes
    for att_name,att_val in IS_gla12_attrs['Data_40HZ'].items():
        if att_name not in ('DS_UTCTime_40',) and not isinstance(att_val,dict):
            fileID['Data_40HZ'].attrs[att_name] = att_val

    # add 40HZ time variable
    val = IS_gla12_tide['Data_40HZ']['DS_UTCTime_40']
    attrs = IS_gla12_attrs['Data_40HZ']['DS_UTCTime_40']
    # Defining the HDF5 dataset variables
    var = '{0}/{1}'.format('Data_40HZ','DS_UTCTime_40')
    h5['Data_40HZ']['DS_UTCTime_40'] = fileID.create_dataset(var,
        np.shape(val), data=val, dtype=val.dtype, compression='gzip')
    # make dimension
    h5['Data_40HZ']['DS_UTCTime_40'].make_scale('DS_UTCTime_40')
    # add HDF5 variable attributes
    for att_name,att_val in attrs.items():
        h5['Data_40HZ']['DS_UTCTime_40'].attrs[att_name] = att_val

    # for each variable group
    for group in ['Time','Geolocation','Geophysical']:
        # add group to dict
        h5['Data_40HZ'][group] = {}
        # create Data_40HZ group
        fileID.create_group('Data_40HZ/{0}'.format(group))
        # add HDF5 group attributes
        for att_name,att_val in IS_gla12_attrs['Data_40HZ'][group].items():
            if not isinstance(att_val,dict):
                fileID['Data_40HZ'][group].attrs[att_name] = att_val
        # for each variable in the group
        for key,val in IS_gla12_tide['Data_40HZ'][group].items():
            fillvalue = FILL_VALUE['Data_40HZ'][group][key]
            attrs = IS_gla12_attrs['Data_40HZ'][group][key]
            # Defining the HDF5 dataset variables
            var = '{0}/{1}/{2}'.format('Data_40HZ',group,key)
            # use variable compression if containing fill values
            if fillvalue:
                h5['Data_40HZ'][group][key] = fileID.create_dataset(var,
                    np.shape(val), data=val, dtype=val.dtype,
                    fillvalue=fillvalue, compression='gzip')
            else:
                h5['Data_40HZ'][group][key] = fileID.create_dataset(var,
                    np.shape(val), data=val, dtype=val.dtype,
                    compression='gzip')
            # attach dimensions
            for i,dim in enumerate(['DS_UTCTime_40']):
                h5['Data_40HZ'][group][key].dims[i].attach_scale(
                    h5['Data_40HZ'][dim])
            # add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5['Data_40HZ'][group][key].attrs[att_name] = att_val

    # Closing the HDF5 file
    fileID.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates tidal elevations for correcting ICESat/GLAS
            L2 GLA12 Antarctic and Greenland Ice Sheet elevation data
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    # command line parameters
    group = parser.add_mutually_exclusive_group(required=True)
    # input ICESat GLAS files
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='ICESat GLA12 file to run')
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

    # run for each input GLA12 file
    for FILE in args.infile:
        compute_tides_ICESat(args.directory, FILE,
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
