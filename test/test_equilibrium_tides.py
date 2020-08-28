#!/usr/bin/env python
u"""
test_equilibrium_tides.py (08/2020)
Download an ATL03 and ATL07 file from NSIDC and compare equilibrium tides values
"""
import os
import pytest
import warnings
import numpy as np
import pyTMD.time
import pyTMD.utilities
import pyTMD.calc_delta_time
import pyTMD.compute_equilibrium_tide
import icesat2_toolkit.utilities
from icesat2_toolkit.read_ICESat2_ATL03 import read_HDF5_ATL03
from icesat2_toolkit.read_ICESat2_ATL07 import read_HDF5_ATL07

#-- PURPOSE: Download an ATL03 file from NSIDC and compare equilibrium tides
def test_ATL03_equilibrium_tides(username,password):
    #-- path to an ATL03 file from NSIDC
    HOST = ['https://n5eil01u.ecs.nsidc.org','ATLAS','ATL03.003','2018.10.14',
        'ATL03_20181014000347_02350101_003_01.h5']
    #-- only download ATL03 file if not currently existing
    if not os.access(HOST[-1], os.F_OK):
        #-- download an ATL03 file from NSIDC
        icesat2_toolkit.utilities.from_nsidc(HOST,username=username,
            password=password,local=HOST[-1],verbose=True)
    #-- read ATL03 file using HDF5 reader
    IS2_atl03_mds,IS2_atl03_attrs,IS2_atl03_beams = read_HDF5_ATL03(HOST[-1],
        ATTRIBUTES=True, VERBOSE=True)
    #-- verify that data is imported correctly
    assert all(gtx in IS2_atl03_mds.keys() for gtx in IS2_atl03_beams)
    #-- number of GPS seconds between the GPS epoch
    #-- and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = IS2_atl03_mds['ancillary_data']['atlas_sdp_gps_epoch']
    #-- for each beam
    for gtx in IS2_atl03_beams:
        #-- read ICESat-2 delta time and latitude
        nref = len(IS2_atl03_mds[gtx]['geolocation']['segment_id'])
        delta_time = IS2_atl03_mds[gtx]['geophys_corr']['delta_time']
        latitude = IS2_atl03_mds[gtx]['geolocation']['reference_photon_lat']
        #-- read ASAS predicted long-period equilibrium tides
        fv = IS2_atl03_attrs[gtx]['geophys_corr']['tide_equilibrium']['_FillValue']
        tide_equilibrium = IS2_atl03_mds[gtx]['geophys_corr']['tide_equilibrium']
        #-- calculate tide time for beam
        gps_seconds = atlas_sdp_gps_epoch + delta_time
        leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)
        tide_time = pyTMD.time.convert_delta_time(gps_seconds-leap_seconds,
            epoch1=(1980,1,6,0,0,0), epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
        #-- interpolate delta times from calendar dates to tide time
        delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
        deltat = pyTMD.calc_delta_time(delta_file, tide_time)
        #-- calculate long-period equilibrium tides
        lpet = pyTMD.compute_equilibrium_tide(tide_time+deltat, latitude)
        #-- calculate differences between computed and data versions
        difference = np.ma.zeros((nref))
        difference.data[:] = lpet - tide_equilibrium
        difference.mask = (tide_equilibrium == fv)
        print(np.max(np.abs(difference)))
        #-- will verify differences between outputs are within tolerance
        eps = np.finfo(np.float16).eps
        if not np.all(difference.mask):
            assert np.all(np.abs(difference) < eps)

#-- PURPOSE: Download an ATL07 file from NSIDC and compare equilibrium tides
def test_ATL07_equilibrium_tides(username,password):
    #-- path to an ATL07 file from NSIDC
    HOST = ['https://n5eil01u.ecs.nsidc.org','ATLAS','ATL07.003','2018.10.14',
        'ATL07-01_20181014000347_02350101_003_02.h5']
    #-- only download ATL07 file if not currently existing
    if not os.access(HOST[-1], os.F_OK):
        #-- download an ATL07 file from NSIDC
        icesat2_toolkit.utilities.from_nsidc(HOST,username=username,
            password=password,local=HOST[-1],verbose=True)
    #-- read ATL07 file using HDF5 reader
    IS2_atl07_mds,IS2_atl07_attrs,IS2_atl07_beams = read_HDF5_ATL07(HOST[-1],
        ATTRIBUTES=True, VERBOSE=True)
    #-- verify that data is imported correctly
    assert all(gtx in IS2_atl07_mds.keys() for gtx in IS2_atl07_beams)
    #-- number of GPS seconds between the GPS epoch
    #-- and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = IS2_atl07_mds['ancillary_data']['atlas_sdp_gps_epoch']
    #-- for each beam
    for gtx in IS2_atl07_beams:
        #-- read ICESat-2 sea ice delta time and latitude
        nseg = len(IS2_atl07_mds[gtx]['sea_ice_segments']['height_segment_id'])
        val = IS2_atl07_mds[gtx]['sea_ice_segments']
        attrs = IS2_atl07_attrs[gtx]['sea_ice_segments']
        delta_time = val['delta_time']
        latitude = val['latitude']
        #-- read ASAS predicted long-period equilibrium tides
        fv = attrs['geophysical']['height_segment_lpe']['_FillValue']
        tide_equilibrium = val['geophysical']['height_segment_lpe'][:]
        #-- calculate tide time for beam
        gps_seconds = atlas_sdp_gps_epoch + delta_time
        leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)
        tide_time = pyTMD.time.convert_delta_time(gps_seconds-leap_seconds,
            epoch1=(1980,1,6,0,0,0), epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
        #-- interpolate delta times from calendar dates to tide time
        delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
        deltat = pyTMD.calc_delta_time(delta_file, tide_time)
        #-- calculate long-period equilibrium tides
        lpet = pyTMD.compute_equilibrium_tide(tide_time+deltat, latitude)
        #-- calculate differences between computed and data versions
        difference = np.ma.zeros((nseg))
        difference.data[:] = lpet - tide_equilibrium
        difference.mask = (tide_equilibrium == fv)
        #-- will verify differences between outputs are within tolerance
        eps = np.finfo(np.float16).eps
        if not np.all(difference.mask):
            assert np.all(np.abs(difference) < eps)
