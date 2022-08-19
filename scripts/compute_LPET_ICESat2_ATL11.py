#!/usr/bin/env python
u"""
compute_LPET_ICESat2_ATL11.py
Written by Tyler Sutterley (07/2022)
Calculates long-period equilibrium tidal elevations for correcting ICESat-2
    annual land ice height data
Will calculate the long-period tides for all ATL11 segments and not just ocean
    segments defined by the ocean tide mask

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
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
    read_ICESat2_ATL11.py: reads ICESat-2 annual land ice height data files
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for syncing files
    calc_delta_time.py: calculates difference between universal and dynamic time
    compute_equilibrium_tide.py: calculates long-period equilibrium ocean tides

UPDATE HISTORY:
    Updated 07/2022: place some imports within try/except statements
    Updated 04/2022: use argparse descriptions within documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 04/2021: can use a generically named ATL11 file as input
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 02/2021: additionally calculate tides for crossing track data
    Updated 01/2021: using standalone ATL11 reader
    Updated 12/2020: merged time conversion routines into module
    Written 12/2020
"""
from __future__ import print_function

import sys
import os
import re
import logging
import argparse
import datetime
import warnings
import numpy as np
import collections
import pyTMD.time
import pyTMD.utilities
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.compute_equilibrium_tide import compute_equilibrium_tide
# attempt imports
try:
    import h5py
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("h5py not available")
try:
    from icesat2_toolkit.read_ICESat2_ATL11 import read_HDF5_ATL11
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("icesat2_toolkit not available")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: read ICESat-2 annual land ice height data (ATL11) from NSIDC
# compute long-period equilibrium tides at points and times
def compute_LPET_ICESat2(INPUT_FILE, VERBOSE=False, MODE=0o775):

    # create logger for verbosity level
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logger = pyTMD.utilities.build_logger('pytmd',level=loglevel)

    # read data from input file
    logger.info('{0} -->'.format(INPUT_FILE))
    IS2_atl11_mds,IS2_atl11_attrs,IS2_atl11_pairs = read_HDF5_ATL11(INPUT_FILE,
        ATTRIBUTES=True, CROSSOVERS=True)
    DIRECTORY = os.path.dirname(INPUT_FILE)
    # extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_'
        r'(\d{3})_(\d{2})(.*?).h5$')
    try:
        SUB,PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX = rx.findall(INPUT_FILE).pop()
    except:
        # output long-period equilibrium tide HDF5 file (generic)
        fileBasename,fileExtension = os.path.splitext(INPUT_FILE)
        OUTPUT_FILE = '{0}_{1}{2}'.format(fileBasename,'LPET',fileExtension)
    else:
        # output long-period equilibrium tide HDF5 file for ASAS/NSIDC granules
        args = (PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX)
        file_format = '{0}_LPET_{1}{2}_{3}{4}_{5}_{6}{7}.h5'
        OUTPUT_FILE = file_format.format(*args)

    # number of GPS seconds between the GPS epoch
    # and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = IS2_atl11_mds['ancillary_data']['atlas_sdp_gps_epoch']
    # delta time (TT - UT1) file
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])

    # copy variables for outputting to HDF5 file
    IS2_atl11_tide = {}
    IS2_atl11_fill = {}
    IS2_atl11_dims = {}
    IS2_atl11_tide_attrs = {}
    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    IS2_atl11_tide['ancillary_data'] = {}
    IS2_atl11_tide_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        # get each HDF5 variable
        IS2_atl11_tide['ancillary_data'][key] = IS2_atl11_mds['ancillary_data'][key]
        # Getting attributes of group and included variables
        IS2_atl11_tide_attrs['ancillary_data'][key] = {}
        for att_name,att_val in IS2_atl11_attrs['ancillary_data'][key].items():
            IS2_atl11_tide_attrs['ancillary_data'][key][att_name] = att_val
    # HDF5 group name for across-track data
    XT = 'crossing_track_data'

    # for each input beam pair within the file
    for ptx in sorted(IS2_atl11_pairs):
        # output data dictionaries for beam
        IS2_atl11_tide[ptx] = dict(cycle_stats=collections.OrderedDict(),
            crossing_track_data=collections.OrderedDict())
        IS2_atl11_fill[ptx] = dict(cycle_stats={},crossing_track_data={})
        IS2_atl11_dims[ptx] = dict(cycle_stats={},crossing_track_data={})
        IS2_atl11_tide_attrs[ptx] = dict(cycle_stats={},crossing_track_data={})

        # extract along-track and across-track variables
        ref_pt = {}
        latitude = {}
        longitude = {}
        delta_time = {}
        # along-track (AT) reference point, latitude, longitude and time
        ref_pt['AT'] = IS2_atl11_mds[ptx]['ref_pt'].copy()
        latitude['AT'] = np.ma.array(IS2_atl11_mds[ptx]['latitude'],
            fill_value=IS2_atl11_attrs[ptx]['latitude']['_FillValue'])
        longitude['AT'] = np.ma.array(IS2_atl11_mds[ptx]['longitude'],
            fill_value=IS2_atl11_attrs[ptx]['longitude']['_FillValue'])
        delta_time['AT'] = np.ma.array(IS2_atl11_mds[ptx]['delta_time'],
            fill_value=IS2_atl11_attrs[ptx]['delta_time']['_FillValue'])
        # across-track (XT) reference point, latitude, longitude and time
        ref_pt['XT'] = IS2_atl11_mds[ptx][XT]['ref_pt'].copy()
        latitude['XT'] = np.ma.array(IS2_atl11_mds[ptx][XT]['latitude'],
            fill_value=IS2_atl11_attrs[ptx][XT]['latitude']['_FillValue'])
        longitude['XT'] = np.ma.array(IS2_atl11_mds[ptx][XT]['longitude'],
            fill_value=IS2_atl11_attrs[ptx][XT]['longitude']['_FillValue'])
        delta_time['XT'] = np.ma.array(IS2_atl11_mds[ptx][XT]['delta_time'],
            fill_value=IS2_atl11_attrs[ptx][XT]['delta_time']['_FillValue'])

        # number of average segments and number of included cycles
        # fill_value for invalid heights and corrections
        fv = IS2_atl11_attrs[ptx]['h_corr']['_FillValue']
        # shape of along-track and across-track data
        n_points,n_cycles = delta_time['AT'].shape
        n_cross, = delta_time['XT'].shape
        # allocate for output long-period equilibrium tide variables
        tide_lpe = {}
        # along-track (AT) tides
        tide_lpe['AT'] = np.ma.empty((n_points,n_cycles),fill_value=fv)
        tide_lpe['AT'].mask = (delta_time['AT'] == delta_time['AT'].fill_value)
        # across-track (XT) tides
        tide_lpe['XT'] = np.ma.empty((n_cross),fill_value=fv)
        tide_lpe['XT'].mask = (delta_time['XT'] == delta_time['XT'].fill_value)

        # calculate tides for along-track and across-track data
        for track in ['AT','XT']:
            # convert time from ATLAS SDP to days relative to Jan 1, 1992
            gps_seconds = atlas_sdp_gps_epoch + delta_time[track]
            leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)
            tide_time = pyTMD.time.convert_delta_time(gps_seconds-leap_seconds,
                epoch1=(1980,1,6,0,0,0), epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
            # interpolate delta times from calendar dates to tide time
            delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
            deltat = calc_delta_time(delta_file, tide_time)

            # calculate  long-period equilibrium tides for track type
            if (track == 'AT'):
                # calculate LPET for each cycle if along-track
                for cycle in range(n_cycles):
                    # find valid time and spatial points for cycle
                    valid, = np.nonzero(~tide_lpe[track].mask[:,cycle])
                    # predict long-period equilibrium tides at latitudes and time
                    t = tide_time[valid,cycle] + deltat[valid,cycle]
                    tide_lpe[track].data[valid,cycle] = compute_equilibrium_tide(t,
                        latitude[track][valid])
            elif (track == 'XT'):
                # find valid time and spatial points for cycle
                valid, = np.nonzero(~tide_lpe[track].mask[:])
                # predict long-period equilibrium tides at latitudes and time
                t = tide_time[valid] + deltat[valid]
                tide_lpe[track].data[valid] = compute_equilibrium_tide(t,
                    latitude[track][valid])

            # replace masked and nan values with fill value
            invalid = np.nonzero(np.isnan(tide_lpe[track].data) | tide_lpe[track].mask)
            tide_lpe[track].data[invalid] = tide_lpe[track].fill_value
            tide_lpe[track].mask[invalid] = True

        # group attributes for beam
        IS2_atl11_tide_attrs[ptx]['description'] = ('Contains the primary science parameters '
            'for this data set')
        IS2_atl11_tide_attrs[ptx]['beam_pair'] = IS2_atl11_attrs[ptx]['beam_pair']
        IS2_atl11_tide_attrs[ptx]['ReferenceGroundTrack'] = IS2_atl11_attrs[ptx]['ReferenceGroundTrack']
        IS2_atl11_tide_attrs[ptx]['first_cycle'] = IS2_atl11_attrs[ptx]['first_cycle']
        IS2_atl11_tide_attrs[ptx]['last_cycle'] = IS2_atl11_attrs[ptx]['last_cycle']
        IS2_atl11_tide_attrs[ptx]['equatorial_radius'] = IS2_atl11_attrs[ptx]['equatorial_radius']
        IS2_atl11_tide_attrs[ptx]['polar_radius'] = IS2_atl11_attrs[ptx]['polar_radius']

        # geolocation, time and reference point
        # reference point
        IS2_atl11_tide[ptx]['ref_pt'] = ref_pt['AT'].copy()
        IS2_atl11_fill[ptx]['ref_pt'] = None
        IS2_atl11_dims[ptx]['ref_pt'] = None
        IS2_atl11_tide_attrs[ptx]['ref_pt'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx]['ref_pt']['units'] = "1"
        IS2_atl11_tide_attrs[ptx]['ref_pt']['contentType'] = "referenceInformation"
        IS2_atl11_tide_attrs[ptx]['ref_pt']['long_name'] = "Reference point number"
        IS2_atl11_tide_attrs[ptx]['ref_pt']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx]['ref_pt']['description'] = ("The reference point is the "
            "7 digit segment_id number corresponding to the center of the ATL06 data used "
            "for each ATL11 point.  These are sequential, starting with 1 for the first "
            "segment after an ascending equatorial crossing node.")
        IS2_atl11_tide_attrs[ptx]['ref_pt']['coordinates'] = \
            "delta_time latitude longitude"
        # cycle_number
        IS2_atl11_tide[ptx]['cycle_number'] = IS2_atl11_mds[ptx]['cycle_number'].copy()
        IS2_atl11_fill[ptx]['cycle_number'] = None
        IS2_atl11_dims[ptx]['cycle_number'] = None
        IS2_atl11_tide_attrs[ptx]['cycle_number'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx]['cycle_number']['units'] = "1"
        IS2_atl11_tide_attrs[ptx]['cycle_number']['long_name'] = "Orbital cycle number"
        IS2_atl11_tide_attrs[ptx]['cycle_number']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx]['cycle_number']['description'] = ("Number of 91-day periods "
            "that have elapsed since ICESat-2 entered the science orbit. Each of the 1,387 "
            "reference ground track (RGTs) is targeted in the polar regions once "
            "every 91 days.")
        # delta time
        IS2_atl11_tide[ptx]['delta_time'] = delta_time['AT'].copy()
        IS2_atl11_fill[ptx]['delta_time'] = delta_time['AT'].fill_value
        IS2_atl11_dims[ptx]['delta_time'] = ['ref_pt','cycle_number']
        IS2_atl11_tide_attrs[ptx]['delta_time'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx]['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl11_tide_attrs[ptx]['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl11_tide_attrs[ptx]['delta_time']['standard_name'] = "time"
        IS2_atl11_tide_attrs[ptx]['delta_time']['calendar'] = "standard"
        IS2_atl11_tide_attrs[ptx]['delta_time']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx]['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl11_tide_attrs[ptx]['delta_time']['coordinates'] = \
            "ref_pt cycle_number latitude longitude"
        # latitude
        IS2_atl11_tide[ptx]['latitude'] = latitude['AT'].copy()
        IS2_atl11_fill[ptx]['latitude'] = latitude['AT'].fill_value
        IS2_atl11_dims[ptx]['latitude'] = ['ref_pt']
        IS2_atl11_tide_attrs[ptx]['latitude'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx]['latitude']['units'] = "degrees_north"
        IS2_atl11_tide_attrs[ptx]['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl11_tide_attrs[ptx]['latitude']['long_name'] = "Latitude"
        IS2_atl11_tide_attrs[ptx]['latitude']['standard_name'] = "latitude"
        IS2_atl11_tide_attrs[ptx]['latitude']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx]['latitude']['description'] = ("Center latitude of "
            "selected segments")
        IS2_atl11_tide_attrs[ptx]['latitude']['valid_min'] = -90.0
        IS2_atl11_tide_attrs[ptx]['latitude']['valid_max'] = 90.0
        IS2_atl11_tide_attrs[ptx]['latitude']['coordinates'] = \
            "ref_pt delta_time longitude"
        # longitude
        IS2_atl11_tide[ptx]['longitude'] = longitude['AT'].copy()
        IS2_atl11_fill[ptx]['longitude'] = longitude['AT'].fill_value
        IS2_atl11_dims[ptx]['longitude'] = ['ref_pt']
        IS2_atl11_tide_attrs[ptx]['longitude'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx]['longitude']['units'] = "degrees_east"
        IS2_atl11_tide_attrs[ptx]['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl11_tide_attrs[ptx]['longitude']['long_name'] = "Longitude"
        IS2_atl11_tide_attrs[ptx]['longitude']['standard_name'] = "longitude"
        IS2_atl11_tide_attrs[ptx]['longitude']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx]['longitude']['description'] = ("Center longitude of "
            "selected segments")
        IS2_atl11_tide_attrs[ptx]['longitude']['valid_min'] = -180.0
        IS2_atl11_tide_attrs[ptx]['longitude']['valid_max'] = 180.0
        IS2_atl11_tide_attrs[ptx]['longitude']['coordinates'] = \
            "ref_pt delta_time latitude"

        # cycle statistics variables
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['Description'] = ("The cycle_stats subgroup "
            "contains summary information about segments for each reference point, including "
            "the uncorrected mean heights for reference surfaces, blowing snow and cloud "
            "indicators, and geolocation and height misfit statistics.")
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['data_rate'] = ("Data within this group "
            "are stored at the average segment rate.")
        # computed long-period equilibrium tide
        IS2_atl11_tide[ptx]['cycle_stats']['tide_equilibrium'] = tide_lpe['AT'].copy()
        IS2_atl11_fill[ptx]['cycle_stats']['tide_equilibrium'] = tide_lpe['AT'].fill_value
        IS2_atl11_dims[ptx]['cycle_stats']['tide_equilibrium'] = ['ref_pt','cycle_number']
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['tide_equilibrium'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['tide_equilibrium']['units'] = "meters"
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['tide_equilibrium']['contentType'] = "referenceInformation"
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['tide_equilibrium']['long_name'] = \
            "Long Period Equilibrium Tide"
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['tide_equilibrium']['description'] = ("Long-period "
            "equilibrium tidal elevation from the summation of fifteen tidal spectral lines")
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['tide_equilibrium']['reference'] = \
            "https://doi.org/10.1111/j.1365-246X.1973.tb03420.x"
        IS2_atl11_tide_attrs[ptx]['cycle_stats']['tide_equilibrium']['coordinates'] = \
            "../ref_pt ../cycle_number ../delta_time ../latitude ../longitude"

        # crossing track variables
        IS2_atl11_tide_attrs[ptx][XT]['Description'] = ("The crossing_track_data "
            "subgroup contains elevation data at crossover locations. These are "
            "locations where two ICESat-2 pair tracks cross, so data are available "
            "from both the datum track, for which the granule was generated, and "
            "from the crossing track.")
        IS2_atl11_tide_attrs[ptx][XT]['data_rate'] = ("Data within this group are "
            "stored at the average segment rate.")

        # reference point
        IS2_atl11_tide[ptx][XT]['ref_pt'] = ref_pt['XT'].copy()
        IS2_atl11_fill[ptx][XT]['ref_pt'] = None
        IS2_atl11_dims[ptx][XT]['ref_pt'] = None
        IS2_atl11_tide_attrs[ptx][XT]['ref_pt'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx][XT]['ref_pt']['units'] = "1"
        IS2_atl11_tide_attrs[ptx][XT]['ref_pt']['contentType'] = "referenceInformation"
        IS2_atl11_tide_attrs[ptx][XT]['ref_pt']['long_name'] = ("fit center reference point number, "
            "segment_id")
        IS2_atl11_tide_attrs[ptx][XT]['ref_pt']['source'] = "derived, ATL11 algorithm"
        IS2_atl11_tide_attrs[ptx][XT]['ref_pt']['description'] = ("The reference-point number of the "
            "fit center for the datum track. The reference point is the 7 digit segment_id number "
            "corresponding to the center of the ATL06 data used for each ATL11 point.  These are "
            "sequential, starting with 1 for the first segment after an ascending equatorial "
            "crossing node.")
        IS2_atl11_tide_attrs[ptx][XT]['ref_pt']['coordinates'] = \
            "delta_time latitude longitude"

        # reference ground track of the crossing track
        IS2_atl11_tide[ptx][XT]['rgt'] = IS2_atl11_mds[ptx][XT]['rgt'].copy()
        IS2_atl11_fill[ptx][XT]['rgt'] = IS2_atl11_attrs[ptx][XT]['rgt']['_FillValue']
        IS2_atl11_dims[ptx][XT]['rgt'] = None
        IS2_atl11_tide_attrs[ptx][XT]['rgt'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx][XT]['rgt']['units'] = "1"
        IS2_atl11_tide_attrs[ptx][XT]['rgt']['contentType'] = "referenceInformation"
        IS2_atl11_tide_attrs[ptx][XT]['rgt']['long_name'] = "crossover reference ground track"
        IS2_atl11_tide_attrs[ptx][XT]['rgt']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx][XT]['rgt']['description'] = "The RGT number for the crossing data."
        IS2_atl11_tide_attrs[ptx][XT]['rgt']['coordinates'] = \
            "ref_pt delta_time latitude longitude"
        # cycle_number of the crossing track
        IS2_atl11_tide[ptx][XT]['cycle_number'] = IS2_atl11_mds[ptx][XT]['cycle_number'].copy()
        IS2_atl11_fill[ptx][XT]['cycle_number'] = IS2_atl11_attrs[ptx][XT]['cycle_number']['_FillValue']
        IS2_atl11_dims[ptx][XT]['cycle_number'] = None
        IS2_atl11_tide_attrs[ptx][XT]['cycle_number'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx][XT]['cycle_number']['units'] = "1"
        IS2_atl11_tide_attrs[ptx][XT]['cycle_number']['long_name'] = "crossover cycle number"
        IS2_atl11_tide_attrs[ptx][XT]['cycle_number']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx][XT]['cycle_number']['description'] = ("Cycle number for the "
            "crossing data. Number of 91-day periods that have elapsed since ICESat-2 entered "
            "the science orbit. Each of the 1,387 reference ground track (RGTs) is targeted "
            "in the polar regions once every 91 days.")
        # delta time of the crossing track
        IS2_atl11_tide[ptx][XT]['delta_time'] = delta_time['XT'].copy()
        IS2_atl11_fill[ptx][XT]['delta_time'] = delta_time['XT'].fill_value
        IS2_atl11_dims[ptx][XT]['delta_time'] = ['ref_pt']
        IS2_atl11_tide_attrs[ptx][XT]['delta_time'] = {}
        IS2_atl11_tide_attrs[ptx][XT]['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl11_tide_attrs[ptx][XT]['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl11_tide_attrs[ptx][XT]['delta_time']['standard_name'] = "time"
        IS2_atl11_tide_attrs[ptx][XT]['delta_time']['calendar'] = "standard"
        IS2_atl11_tide_attrs[ptx][XT]['delta_time']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx][XT]['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl11_tide_attrs[ptx]['delta_time']['coordinates'] = \
            "ref_pt latitude longitude"
        # latitude of the crossover measurement
        IS2_atl11_tide[ptx][XT]['latitude'] = latitude['XT'].copy()
        IS2_atl11_fill[ptx][XT]['latitude'] = latitude['XT'].fill_value
        IS2_atl11_dims[ptx][XT]['latitude'] = ['ref_pt']
        IS2_atl11_tide_attrs[ptx][XT]['latitude'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['units'] = "degrees_north"
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['long_name'] = "crossover latitude"
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['standard_name'] = "latitude"
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['description'] = ("Center latitude of "
            "selected segments")
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['valid_min'] = -90.0
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['valid_max'] = 90.0
        IS2_atl11_tide_attrs[ptx][XT]['latitude']['coordinates'] = \
            "ref_pt delta_time longitude"
        # longitude of the crossover measurement
        IS2_atl11_tide[ptx][XT]['longitude'] = longitude['XT'].copy()
        IS2_atl11_fill[ptx][XT]['longitude'] = longitude['XT'].fill_value
        IS2_atl11_dims[ptx][XT]['longitude'] = ['ref_pt']
        IS2_atl11_tide_attrs[ptx][XT]['longitude'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['units'] = "degrees_east"
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['long_name'] = "crossover longitude"
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['standard_name'] = "longitude"
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['source'] = "ATL06"
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['description'] = ("Center longitude of "
            "selected segments")
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['valid_min'] = -180.0
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['valid_max'] = 180.0
        IS2_atl11_tide_attrs[ptx][XT]['longitude']['coordinates'] = \
            "ref_pt delta_time latitude"
        # computed long-period equilibrium tide for the crossover measurement
        IS2_atl11_tide[ptx][XT]['tide_equilibrium'] = tide_lpe['XT'].copy()
        IS2_atl11_fill[ptx][XT]['tide_equilibrium'] = tide_lpe['XT'].fill_value
        IS2_atl11_dims[ptx][XT]['tide_equilibrium'] = ['ref_pt']
        IS2_atl11_tide_attrs[ptx][XT]['tide_equilibrium'] = collections.OrderedDict()
        IS2_atl11_tide_attrs[ptx][XT]['tide_equilibrium']['units'] = "meters"
        IS2_atl11_tide_attrs[ptx][XT]['tide_equilibrium']['contentType'] = "referenceInformation"
        IS2_atl11_tide_attrs[ptx][XT]['tide_equilibrium']['long_name'] = \
            "Long Period Equilibrium Tide"
        IS2_atl11_tide_attrs[ptx][XT]['tide_equilibrium']['description'] = ("Long-period "
            "equilibrium tidal elevation from the summation of fifteen tidal spectral lines")
        IS2_atl11_tide_attrs[ptx][XT]['tide_equilibrium']['reference'] = \
            "https://doi.org/10.1111/j.1365-246X.1973.tb03420.x"
        IS2_atl11_tide_attrs[ptx][XT]['tide_equilibrium']['coordinates'] = \
            "ref_pt delta_time latitude longitude"

    # print file information
    logger.info('\t{0}'.format(OUTPUT_FILE))
    HDF5_ATL11_tide_write(IS2_atl11_tide, IS2_atl11_tide_attrs,
        CLOBBER=True, INPUT=os.path.basename(INPUT_FILE),
        FILL_VALUE=IS2_atl11_fill, DIMENSIONS=IS2_atl11_dims,
        FILENAME=os.path.join(DIRECTORY,OUTPUT_FILE))
    # change the permissions mode
    os.chmod(os.path.join(DIRECTORY,OUTPUT_FILE), MODE)

# PURPOSE: outputting the tide values for ICESat-2 data to HDF5
def HDF5_ATL11_tide_write(IS2_atl11_tide, IS2_atl11_attrs, INPUT=None,
    FILENAME='', FILL_VALUE=None, DIMENSIONS=None, CLOBBER=False):
    # setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    # open output HDF5 file
    fileID = h5py.File(os.path.expanduser(FILENAME), clobber)

    # create HDF5 records
    h5 = {}

    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    h5['ancillary_data'] = {}
    for k,v in IS2_atl11_tide['ancillary_data'].items():
        # Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        # add HDF5 variable attributes
        for att_name,att_val in IS2_atl11_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    # write each output beam pair
    pairs = [k for k in IS2_atl11_tide.keys() if bool(re.match(r'pt\d',k))]
    for ptx in pairs:
        fileID.create_group(ptx)
        h5[ptx] = {}
        # add HDF5 group attributes for beam
        for att_name in ['description','beam_pair','ReferenceGroundTrack',
            'first_cycle','last_cycle','equatorial_radius','polar_radius']:
            fileID[ptx].attrs[att_name] = IS2_atl11_attrs[ptx][att_name]

        # ref_pt, cycle number, geolocation and delta_time variables
        for k in ['ref_pt','cycle_number','delta_time','latitude','longitude']:
            # values and attributes
            v = IS2_atl11_tide[ptx][k]
            attrs = IS2_atl11_attrs[ptx][k]
            fillvalue = FILL_VALUE[ptx][k]
            # Defining the HDF5 dataset variables
            val = '{0}/{1}'.format(ptx,k)
            if fillvalue:
                h5[ptx][k] = fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
            else:
                h5[ptx][k] = fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, compression='gzip')
            # create or attach dimensions for HDF5 variable
            if DIMENSIONS[ptx][k]:
                # attach dimensions
                for i,dim in enumerate(DIMENSIONS[ptx][k]):
                    h5[ptx][k].dims[i].attach_scale(h5[ptx][dim])
            else:
                # make dimension
                h5[ptx][k].make_scale(k)
            # add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[ptx][k].attrs[att_name] = att_val

        # add to cycle_stats and crossing_track_data variables
        for key in ['cycle_stats','crossing_track_data']:
            fileID[ptx].create_group(key)
            h5[ptx][key] = {}
            for att_name in ['Description','data_rate']:
                att_val=IS2_atl11_attrs[ptx][key][att_name]
                fileID[ptx][key].attrs[att_name] = att_val
            for k,v in IS2_atl11_tide[ptx][key].items():
                # attributes
                attrs = IS2_atl11_attrs[ptx][key][k]
                fillvalue = FILL_VALUE[ptx][key][k]
                # Defining the HDF5 dataset variables
                val = '{0}/{1}/{2}'.format(ptx,key,k)
                if fillvalue:
                    h5[ptx][key][k] = fileID.create_dataset(val, np.shape(v), data=v,
                        dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
                else:
                    h5[ptx][key][k] = fileID.create_dataset(val, np.shape(v), data=v,
                        dtype=v.dtype, compression='gzip')
                # create or attach dimensions for HDF5 variable
                if DIMENSIONS[ptx][key][k]:
                    # attach dimensions
                    for i,dim in enumerate(DIMENSIONS[ptx][key][k]):
                        if (key == 'cycle_stats'):
                            h5[ptx][key][k].dims[i].attach_scale(h5[ptx][dim])
                        else:
                            h5[ptx][key][k].dims[i].attach_scale(h5[ptx][key][dim])
                else:
                    # make dimension
                    h5[ptx][key][k].make_scale(k)
                # add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[ptx][key][k].attrs[att_name] = att_val

    # HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 Annual Land Ice Height'
    fileID.attrs['summary'] = ('The purpose of ATL11 is to provide an ICESat-2 '
        'satellite cycle summary of heights and height changes of land-based '
        'ice and will be provided as input to ATL15 and ATL16, gridded '
        'estimates of heights and height-changes.')
    fileID.attrs['description'] = ('Land ice parameters for each beam pair. '
        'All parameters are calculated for the same along-track increments '
        'for each beam pair and repeat.')
    date_created = datetime.datetime.today()
    fileID.attrs['date_created'] = date_created.isoformat()
    project = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = project
    platform = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = platform
    # add attribute for elevation instrument and designated processing level
    instrument = 'ATLAS > Advanced Topographic Laser Altimeter System'
    fileID.attrs['instrument'] = instrument
    fileID.attrs['source'] = 'Spacecraft'
    fileID.attrs['references'] = 'https://nsidc.org/data/icesat-2'
    fileID.attrs['processing_level'] = '4'
    # add attributes for input ATL11 files
    fileID.attrs['input_files'] = os.path.basename(INPUT)
    # find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for ptx in pairs:
        lon = IS2_atl11_tide[ptx]['longitude']
        lat = IS2_atl11_tide[ptx]['latitude']
        delta_time = IS2_atl11_tide[ptx]['delta_time']
        valid = np.nonzero(delta_time != FILL_VALUE[ptx]['delta_time'])
        # setting the geospatial and temporal ranges
        lnmn = lon.min() if (lon.min() < lnmn) else lnmn
        lnmx = lon.max() if (lon.max() > lnmx) else lnmx
        ltmn = lat.min() if (lat.min() < ltmn) else ltmn
        ltmx = lat.max() if (lat.max() > ltmx) else ltmx
        tmn = delta_time[valid].min() if (delta_time[valid].min() < tmn) else tmn
        tmx = delta_time[valid].max() if (delta_time[valid].max() > tmx) else tmx
    # add geospatial and temporal attributes
    fileID.attrs['geospatial_lat_min'] = ltmn
    fileID.attrs['geospatial_lat_max'] = ltmx
    fileID.attrs['geospatial_lon_min'] = lnmn
    fileID.attrs['geospatial_lon_max'] = lnmx
    fileID.attrs['geospatial_lat_units'] = "degrees_north"
    fileID.attrs['geospatial_lon_units'] = "degrees_east"
    fileID.attrs['geospatial_ellipsoid'] = "WGS84"
    fileID.attrs['date_type'] = 'UTC'
    fileID.attrs['time_type'] = 'CCSDS UTC-A'
    # convert start and end time from ATLAS SDP seconds into GPS seconds
    atlas_sdp_gps_epoch=IS2_atl11_tide['ancillary_data']['atlas_sdp_gps_epoch']
    gps_seconds = atlas_sdp_gps_epoch + np.array([tmn,tmx])
    # calculate leap seconds
    leaps = pyTMD.time.count_leap_seconds(gps_seconds)
    # convert from seconds since 1980-01-06T00:00:00 to Julian days
    time_julian = 2400000.5 + pyTMD.time.convert_delta_time(gps_seconds - leaps,
        epoch1=(1980,1,6,0,0,0), epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    # convert to calendar date
    YY,MM,DD,HH,MN,SS = pyTMD.time.convert_julian(time_julian,format='tuple')
    # add attributes with measurement date start, end and duration
    tcs = datetime.datetime(int(YY[0]), int(MM[0]), int(DD[0]),
        int(HH[0]), int(MN[0]), int(SS[0]), int(1e6*(SS[0] % 1)))
    fileID.attrs['time_coverage_start'] = tcs.isoformat()
    tce = datetime.datetime(int(YY[1]), int(MM[1]), int(DD[1]),
        int(HH[1]), int(MN[1]), int(SS[1]), int(1e6*(SS[1] % 1)))
    fileID.attrs['time_coverage_end'] = tce.isoformat()
    fileID.attrs['time_coverage_duration'] = '{0:0.0f}'.format(tmx-tmn)
    # Closing the HDF5 file
    fileID.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates long-period equilibrium tidal elevations for
            correcting ICESat-2 ATL11 annual land ice height data
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='ICESat-2 ATL11 file to run')
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

    # run for each input ATL11 file
    for FILE in args.infile:
        compute_LPET_ICESat2(FILE, VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
