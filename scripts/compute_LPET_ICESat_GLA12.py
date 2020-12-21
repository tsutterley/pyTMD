#!/usr/bin/env python
u"""
compute_LPET_ICESat_GLA12.py
Written by Tyler Sutterley (12/2020)
Calculates long-period equilibrium tidal elevations for correcting
    ICESat/GLAS L2 GLA12 Antarctic and Greenland Ice Sheet elevation data
Will calculate the long-period tides for all GLAS elevations and not just
    ocean elevations defined by the ocean tide mask

COMMAND LINE OPTIONS:
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
    spatial.py: utilities for reading, writing and operating on spatial data
    utilities: download and management utilities for syncing files
    calc_delta_time.py: calculates difference between universal and dynamic time
    compute_equilibrium_tide.py: calculates long-period equilibrium ocean tides

UPDATE HISTORY:
    Updated 12/2020: H5py deprecation warning change to use make_scale
    Written 12/2020
"""
from __future__ import print_function

import sys
import os
import re
import h5py
import argparse
import numpy as np
import pyTMD.time
import pyTMD.spatial
import pyTMD.utilities
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.compute_equilibrium_tide import compute_equilibrium_tide

#-- PURPOSE: read ICESat ice sheet HDF5 elevation data (GLAH12) from NSIDC
#-- compute long-period equilibrium tides at points and times
def compute_LPET_ICESat(FILE, VERBOSE=False, MODE=0o775):

    #-- get directory from FILE
    print('{0} -->'.format(os.path.basename(FILE))) if VERBOSE else None
    DIRECTORY = os.path.dirname(FILE)

    #-- compile regular expression operator for extracting information from file
    rx = re.compile((r'GLAH(\d{2})_(\d{3})_(\d{1})(\d{1})(\d{2})_(\d{3})_'
        r'(\d{4})_(\d{1})_(\d{2})_(\d{4})\.H5'), re.VERBOSE)
    #-- extract parameters from ICESat/GLAS HDF5 file name
    #-- PRD:  Product number (01, 05, 06, 12, 13, 14, or 15)
    #-- RL:  Release number for process that created the product = 634
    #-- RGTP:  Repeat ground-track phase (1=8-day, 2=91-day, 3=transfer orbit)
    #-- ORB:   Reference orbit number (starts at 1 and increments each time a
    #--           new reference orbit ground track file is obtained.)
    #-- INST:  Instance number (increments every time the satellite enters a
    #--           different reference orbit)
    #-- CYCL:   Cycle of reference orbit for this phase
    #-- TRK: Track within reference orbit
    #-- SEG:   Segment of orbit
    #-- GRAN:  Granule version number
    #-- TYPE:  File type
    PRD,RL,RGTP,ORB,INST,CYCL,TRK,SEG,GRAN,TYPE = rx.findall(FILE).pop()

    #-- read GLAH12 HDF5 file
    fileID = h5py.File(FILE,'r')
    n_40HZ, = fileID['Data_40HZ']['Time']['i_rec_ndx'].shape
    #-- get variables and attributes
    rec_ndx_40HZ = fileID['Data_40HZ']['Time']['i_rec_ndx'][:].copy()
    #-- seconds since 2000-01-01 12:00:00 UTC (J2000)
    DS_UTCTime_40HZ = fileID['Data_40HZ']['DS_UTCTime_40'][:].copy()
    #-- Latitude (degrees North)
    lat_TPX = fileID['Data_40HZ']['Geolocation']['d_lat'][:].copy()
    #-- Longitude (degrees East)
    lon_40HZ = fileID['Data_40HZ']['Geolocation']['d_lon'][:].copy()
    #-- Elevation (height above TOPEX/Poseidon ellipsoid in meters)
    elev_TPX = fileID['Data_40HZ']['Elevation_Surfaces']['d_elev'][:].copy()
    fv = fileID['Data_40HZ']['Elevation_Surfaces']['d_elev'].attrs['_FillValue']

    #-- semimajor axis (a) and flattening (f) for TP and WGS84 ellipsoids
    atop,ftop = (6378136.3,1.0/298.257)
    awgs,fwgs = (6378137.0,1.0/298.257223563)
    #-- convert from Topex/Poseidon to WGS84 Ellipsoids
    lat_40HZ,elev_40HZ = pyTMD.spatial.convert_ellipsoid(lat_TPX, elev_TPX,
        atop, ftop, awgs, fwgs, eps=1e-12, itmax=10)

    #-- convert time from J2000 to days relative to Jan 1, 1992 (48622mjd)
    #-- J2000: seconds since 2000-01-01 12:00:00 UTC
    tide_time = pyTMD.time.convert_delta_time(DS_UTCTime_40HZ,
        epoch1=(2000,1,1,12,0,0), epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
    #-- interpolate delta times from calendar dates to tide time
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
    deltat = calc_delta_time(delta_file, tide_time)

    #-- predict long-period equilibrium tides at latitudes and time
    tide_lpe = compute_equilibrium_tide(tide_time + deltat, lat_40HZ)

    #-- copy variables for outputting to HDF5 file
    IS_gla12_tide = dict(Data_40HZ={})
    IS_gla12_fill = dict(Data_40HZ={})
    IS_gla12_tide_attrs = dict(Data_40HZ={})

    #-- copy global file attributes
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

    #-- add attributes for input GLA12 file
    IS_gla12_tide_attrs['input_files'] = os.path.basename(FILE)
    #-- update geospatial ranges for ellipsoid
    IS_gla12_tide_attrs['geospatial_lat_min'] = np.min(lat_40HZ)
    IS_gla12_tide_attrs['geospatial_lat_max'] = np.max(lat_40HZ)
    IS_gla12_tide_attrs['geospatial_lon_min'] = np.min(lon_40HZ)
    IS_gla12_tide_attrs['geospatial_lon_max'] = np.max(lon_40HZ)
    IS_gla12_tide_attrs['geospatial_lat_units'] = "degrees_north"
    IS_gla12_tide_attrs['geospatial_lon_units'] = "degrees_east"
    IS_gla12_tide_attrs['geospatial_ellipsoid'] = "WGS84"

    #-- copy 40Hz group attributes
    for att_name,att_val in fileID['Data_40HZ'].attrs.items():
        IS_gla12_tide_attrs['Data_40HZ'][att_name] = att_val
    #-- copy attributes for time, geolocation and geophysical groups
    for var in ['Time','Geolocation','Geophysical']:
        IS_gla12_tide['Data_40HZ'][var] = {}
        IS_gla12_fill['Data_40HZ'][var] = {}
        IS_gla12_tide_attrs['Data_40HZ'][var] = {}
        for att_name,att_val in fileID['Data_40HZ'][var].attrs.items():
            IS_gla12_tide_attrs['Data_40HZ'][var][att_name] = att_val

    #-- J2000 time
    IS_gla12_tide['Data_40HZ']['DS_UTCTime_40'] = DS_UTCTime_40HZ
    IS_gla12_fill['Data_40HZ']['DS_UTCTime_40'] = None
    IS_gla12_tide_attrs['Data_40HZ']['DS_UTCTime_40'] = {}
    for att_name,att_val in fileID['Data_40HZ']['DS_UTCTime_40'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['DS_UTCTime_40'][att_name] = att_val
    #-- record
    IS_gla12_tide['Data_40HZ']['Time']['i_rec_ndx'] = rec_ndx_40HZ
    IS_gla12_fill['Data_40HZ']['Time']['i_rec_ndx'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Time']['i_rec_ndx'] = {}
    for att_name,att_val in fileID['Data_40HZ']['Time']['i_rec_ndx'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Time']['i_rec_ndx'][att_name] = att_val
    #-- latitude
    IS_gla12_tide['Data_40HZ']['Geolocation']['d_lat'] = lat_40HZ
    IS_gla12_fill['Data_40HZ']['Geolocation']['d_lat'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lat'] = {}
    for att_name,att_val in fileID['Data_40HZ']['Geolocation']['d_lat'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lat'][att_name] = att_val
    #-- longitude
    IS_gla12_tide['Data_40HZ']['Geolocation']['d_lon'] = lon_40HZ
    IS_gla12_fill['Data_40HZ']['Geolocation']['d_lon'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lon'] = {}
    for att_name,att_val in fileID['Data_40HZ']['Geolocation']['d_lon'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lon'][att_name] = att_val

    #-- geophysical variables
    #-- computed long-period equilibrium tide
    IS_gla12_tide['Data_40HZ']['Geophysical']['d_eqElv'] = tide_lpe
    IS_gla12_fill['Data_40HZ']['Geophysical']['d_eqElv'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical']['d_eqElv'] = {}
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical']['d_eqElv']['units'] = "meters"
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical']['d_eqElv']['long_name'] = \
        "Long Period Equilibrium Tide"
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical']['d_eqElv']['description'] = ("Long-period "
            "equilibrium tidal elevation from the summation of fifteen tidal spectral lines")
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical']['d_eqElv']['reference'] = \
        "https://doi.org/10.1111/j.1365-246X.1973.tb03420.x"
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical']['d_eqElv']['coordinates'] = \
        "../DS_UTCTime_40"

    #-- close the input HDF5 file
    fileID.close()

    #-- output tidal HDF5 file
    args = (PRD,RL,RGTP,ORB,INST,CYCL,TRK,SEG,GRAN,TYPE)
    file_format = 'GLAH{0}_{1}_LPET_{2}{3}{4}_{5}_{6}_{7}_{8}_{9}.h5'
    #-- print file information
    print('\t{0}'.format(file_format.format(*args))) if VERBOSE else None
    HDF5_GLA12_tide_write(IS_gla12_tide, IS_gla12_tide_attrs,
        FILENAME=os.path.join(DIRECTORY,file_format.format(*args)),
        FILL_VALUE=IS_gla12_fill, CLOBBER=True)
    #-- change the permissions mode
    os.chmod(os.path.join(DIRECTORY,file_format.format(*args)), MODE)

#-- PURPOSE: outputting the tide values for ICESat data to HDF5
def HDF5_GLA12_tide_write(IS_gla12_tide, IS_gla12_attrs,
    FILENAME='', FILL_VALUE=None, CLOBBER=False):
    #-- setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    #-- open output HDF5 file
    fileID = h5py.File(os.path.expanduser(FILENAME), clobber)
    #-- create 40HZ HDF5 records
    h5 = dict(Data_40HZ={})

    #-- add HDF5 file attributes
    attrs = {a:v for a,v in IS_gla12_attrs.items() if not isinstance(v,dict)}
    for att_name,att_val in attrs.items():
       fileID.attrs[att_name] = att_val

    #-- create Data_40HZ group
    fileID.create_group('Data_40HZ')
    #-- add HDF5 40HZ group attributes
    for att_name,att_val in IS_gla12_attrs['Data_40HZ'].items():
        if att_name not in ('DS_UTCTime_40',) and not isinstance(att_val,dict):
            fileID['Data_40HZ'].attrs[att_name] = att_val

    #-- add 40HZ time variable
    val = IS_gla12_tide['Data_40HZ']['DS_UTCTime_40']
    attrs = IS_gla12_attrs['Data_40HZ']['DS_UTCTime_40']
    #-- Defining the HDF5 dataset variables
    var = '{0}/{1}'.format('Data_40HZ','DS_UTCTime_40')
    h5['Data_40HZ']['DS_UTCTime_40'] = fileID.create_dataset(var,
        np.shape(val), data=val, dtype=val.dtype, compression='gzip')
    #-- make dimension
    h5['Data_40HZ']['DS_UTCTime_40'].make_scale('DS_UTCTime_40')
    #-- add HDF5 variable attributes
    for att_name,att_val in attrs.items():
        h5['Data_40HZ']['DS_UTCTime_40'].attrs[att_name] = att_val

    #-- for each variable group
    for group in ['Time','Geolocation','Geophysical']:
        #-- add group to dict
        h5['Data_40HZ'][group] = {}
        #-- create Data_40HZ group
        fileID.create_group('Data_40HZ/{0}'.format(group))
        #-- add HDF5 group attributes
        for att_name,att_val in IS_gla12_attrs['Data_40HZ'][group].items():
            if not isinstance(att_val,dict):
                fileID['Data_40HZ'][group].attrs[att_name] = att_val
        #-- for each variable in the group
        for key,val in IS_gla12_tide['Data_40HZ'][group].items():
            fillvalue = FILL_VALUE['Data_40HZ'][group][key]
            attrs = IS_gla12_attrs['Data_40HZ'][group][key]
            #-- Defining the HDF5 dataset variables
            var = '{0}/{1}/{2}'.format('Data_40HZ',group,key)
            #-- use variable compression if containing fill values
            if fillvalue:
                h5['Data_40HZ'][group][key] = fileID.create_dataset(var,
                    np.shape(val), data=val, dtype=val.dtype,
                    fillvalue=fillvalue, compression='gzip')
            else:
                h5['Data_40HZ'][group][key] = fileID.create_dataset(var,
                    np.shape(val), data=val, dtype=val.dtype,
                    compression='gzip')
            #-- attach dimensions
            for i,dim in enumerate(['DS_UTCTime_40']):
                h5['Data_40HZ'][group][key].dims[i].attach_scale(
                    h5['Data_40HZ'][dim])
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5['Data_40HZ'][group][key].attrs[att_name] = att_val

    #-- Closing the HDF5 file
    fileID.close()

#-- Main program that calls compute_LPET_ICESat()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Calculates long-period equilibrium tidal elevations for
            correcting ICESat/GLAS L2 GLA12 Antarctic and Greenland Ice Sheet
            elevation data
            """
    )
    #-- command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='ICESat GLA12 file to run')
    #-- verbosity settings
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Output information about each created file')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files created')
    args = parser.parse_args()

    #-- run for each input GLA12 file
    for FILE in args.infile:
        compute_LPET_ICESat(FILE, VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
