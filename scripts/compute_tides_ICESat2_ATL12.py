#!/usr/bin/env python
u"""
compute_tides_ICESat2_ATL12.py
Written by Tyler Sutterley (07/2020)
Calculates tidal elevations for correcting ICESat-2 ocean surface height data

Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
or Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

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
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    read_ICESat2_ATL12.py: reads ICESat-2 ocean surface height data files
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

UPDATE HISTORY:
    Updated 08/2020: using builtin time operations
    Updated 07/2020: added FES2014 and FES2014_load.  use merged delta times
    Updated 06/2020: added version 2 of TPX09-atlas (TPX09-atlas-v2)
    Updated 03/2020: use read_ICESat2_ATL12.py from read-ICESat-2 repository
    Updated 02/2020: changed CATS2008 grid to match version on U.S. Antarctic
        Program Data Center http://www.usap-dc.org/view/dataset/601235
    Forked 12/2019 from compute_tides_ICESat2_ATL07.py
    Updated 11/2019: added AOTIM-5-2018 tide model (2018 update to 2004 model)
    Forked 11/2019 from compute_tides_ICESat2_atl06.py
    Updated 10/2019: external read functions.  adjust regex for processed files
        changing Y/N flags to True/False
    Updated 09/2019: using date functions paralleling public repository
        add option for TPXO9-atlas.  add OTIS netcdf tide option
    Updated 05/2019: check if beam exists in a try except else clause
    Updated 04/2019: check if subsetted beam contains land ice data
    Written 04/2019
"""
from __future__ import print_function

import sys
import os
import re
import time
import h5py
import getopt
import datetime
import numpy as np
import pyTMD.time
import pyTMD.utilities
from pyTMD.convert_julian import convert_julian
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants
from pyTMD.read_FES_model import extract_FES_constants
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide_drift import predict_tide_drift
from icesat2_toolkit.read_ICESat2_ATL12 import read_HDF5_ATL12

#-- PURPOSE: read ICESat-2 ocean surface height (ATL12) from NSIDC
#-- compute tides at points and times using tidal model driver algorithms
def compute_tides_ICESat2(tide_dir,FILE,MODEL,VERBOSE=False,MODE=0o775):
    #-- select between tide models
    if (MODEL == 'CATS0201'):
        grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
        model_file = os.path.join(tide_dir,'cats0201_tmd','h0_CATS02_01')
        reference = 'https://mail.esr.org/polar_tide_models/Model_CATS0201.html'
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'CATS2008'):
        grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008')
        model_file = os.path.join(tide_dir,'CATS2008','hf.CATS2008.out')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = 'CATS2008'
        TYPE = 'z'
    elif (MODEL == 'CATS2008_load'):
        grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
        model_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        variable = 'tide_load_seg'
        long_name = "Load Tide"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
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
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
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
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'netcdf'
        TYPE = 'z'
        SCALE = 1.0/1000.0
    elif (MODEL == 'TPXO9.1'):
        grid_file = os.path.join(tide_dir,'TPXO9.1','DATA','grid_tpxo9')
        model_file = os.path.join(tide_dir,'TPXO9.1','DATA','h_tpxo9.v1')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_file = os.path.join(tide_dir,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
        reference = 'http://volkov.oce.orst.edu/tides/tpxo8_atlas.html'
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'ATLAS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'TPXO7.2'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_tmd','h_tpxo7.2')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'TPXO7.2_load'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_load','grid_tpxo6.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_load','h_tpxo7.2_load')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        variable = 'tide_load_seg'
        long_name = "Load Tide"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (MODEL == 'AODTM-5'):
        grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aodtm5_tmd','h0_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aodtm-5/')
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPE = 'z'
    elif (MODEL == 'AOTIM-5'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPE = 'z'
    elif (MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(tide_dir,'Arc5km2018','grid_Arc5km2018')
        model_file = os.path.join(tide_dir,'Arc5km2018','h_Arc5km2018')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
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
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
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
        variable = 'tide_load_seg'
        long_name = "Load Tide"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'GOT4.8'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
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
        variable = 'tide_load_seg'
        long_name = "Load Tide"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'GOT4.10'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
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
        variable = 'tide_load_seg'
        long_name = "Load Tide"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
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
        variable = 'tide_ocean_seg'
        long_name = "Ocean Tide"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
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
        variable = 'tide_load_seg'
        long_name = "Load Tide"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'FES'
        TYPE = 'z'
        SCALE = 1.0/100.0

    #-- read data from FILE
    print('{0} -->'.format(os.path.basename(FILE))) if VERBOSE else None
    IS2_atl12_mds,IS2_atl12_attrs,IS2_atl12_beams = read_HDF5_ATL12(FILE,
        ATTRIBUTES=True)
    DIRECTORY = os.path.dirname(FILE)
    #-- extract parameters from ICESat-2 ATLAS HDF5 ocean surface file name
    rx = re.compile('(processed_)?(ATL\d{2})_(\d{4})(\d{2})(\d{2})(\d{2})'
        '(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX = rx.findall(FILE).pop()

    #-- number of GPS seconds between the GPS epoch
    #-- and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = IS2_atl12_mds['ancillary_data']['atlas_sdp_gps_epoch']

    #-- copy variables for outputting to HDF5 file
    IS2_atl12_tide = {}
    IS2_atl12_fill = {}
    IS2_atl12_tide_attrs = {}
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    IS2_atl12_tide['ancillary_data'] = {}
    IS2_atl12_tide_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl12_tide['ancillary_data'][key] = IS2_atl12_mds['ancillary_data'][key]
        #-- Getting attributes of group and included variables
        IS2_atl12_tide_attrs['ancillary_data'][key] = {}
        for att_name,att_val in IS2_atl12_attrs['ancillary_data'][key].items():
            IS2_atl12_tide_attrs['ancillary_data'][key][att_name] = att_val

    #-- for each input beam within the file
    for gtx in sorted(IS2_atl12_beams):
        #-- output data dictionaries for beam
        IS2_atl12_tide[gtx] = dict(ssh_segments={})
        IS2_atl12_fill[gtx] = dict(ssh_segments={})
        IS2_atl12_tide_attrs[gtx] = dict(ssh_segments={})

        #-- number of segments
        val = IS2_atl12_mds[gtx]['ssh_segments']
        n_seg = len(val['delt_seg'])

        #-- convert time from ATLAS SDP to days relative to Jan 1, 1992
        gps_seconds = atlas_sdp_gps_epoch + val['delta_time']
        leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)
        tide_time = pyTMD.time.convert_delta_time(gps_seconds-leap_seconds,
            epoch1=(1980,1,6,0,0,0), epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
        #-- read tidal constants and interpolate to grid points
        if model_format in ('OTIS','ATLAS'):
            amp,ph,D,c = extract_tidal_constants(val['longitude'],
                val['latitude'], grid_file, model_file, EPSG, TYPE=TYPE,
                METHOD='spline', GRID=model_format)
            deltat = np.zeros_like(tide_time)
        elif (model_format == 'netcdf'):
            amp,ph,D,c = extract_netcdf_constants(val['longitude'],
                val['latitude'], model_directory, grid_file,
                model_files, TYPE=TYPE, METHOD='spline', SCALE=SCALE)
            deltat = np.zeros_like(tide_time)
        elif (model_format == 'GOT'):
            amp,ph = extract_GOT_constants(val['longitude'], val['latitude'],
                model_directory, model_files, METHOD='spline', SCALE=SCALE)
            #-- interpolate delta times from calendar dates to tide time
            delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
            deltat = calc_delta_time(delta_file, tide_time)
        elif (model_format == 'FES'):
            amp,ph = extract_FES_constants(val['longitude'], val['latitude'],
                model_directory, model_files, TYPE=TYPE, VERSION=MODEL,
                METHOD=METHOD, SCALE=SCALE)
            deltat = np.zeros_like(tide_time)

        #-- calculate complex phase in radians for Euler's
        cph = -1j*ph*np.pi/180.0
        #-- calculate constituent oscillation
        hc = amp*np.exp(cph)

        #-- predict tidal elevations at time and infer minor corrections
        tide = np.ma.empty((n_seg))
        tide.mask = np.any(hc.mask,axis=1)
        tide.data[:] = predict_tide_drift(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=model_format)
        minor = infer_minor_corrections(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=model_format)
        tide.data[:] += minor.data[:]
        #-- replace masked and nan values with fill value
        invalid, = np.nonzero(np.isnan(tide.data) | tide.mask)
        tide.data[invalid] = tide.fill_value
        tide.mask[invalid] = True

        #-- group attributes for beam
        IS2_atl12_tide_attrs[gtx]['Description'] = IS2_atl12_attrs[gtx]['Description']
        IS2_atl12_tide_attrs[gtx]['atlas_pce'] = IS2_atl12_attrs[gtx]['atlas_pce']
        IS2_atl12_tide_attrs[gtx]['atlas_beam_type'] = IS2_atl12_attrs[gtx]['atlas_beam_type']
        IS2_atl12_tide_attrs[gtx]['groundtrack_id'] = IS2_atl12_attrs[gtx]['groundtrack_id']
        IS2_atl12_tide_attrs[gtx]['atmosphere_profile'] = IS2_atl12_attrs[gtx]['atmosphere_profile']
        IS2_atl12_tide_attrs[gtx]['atlas_spot_number'] = IS2_atl12_attrs[gtx]['atlas_spot_number']
        IS2_atl12_tide_attrs[gtx]['sc_orientation'] = IS2_atl12_attrs[gtx]['sc_orientation']
        #-- group attributes for ssh_segments
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['Description'] = ("Contains "
            "parameters relating to the calculated surface height.")
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['data_rate'] = ("Data within "
            "this group are stored at the variable ocean processing segment rate.")


        #-- geolocation, time and segment ID
        #-- delta time
        IS2_atl12_tide[gtx]['ssh_segments']['delta_time'] = val['delta_time'].copy()
        IS2_atl12_fill[gtx]['ssh_segments']['delta_time'] = None
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delta_time'] = {}
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delta_time']['standard_name'] = "time"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delta_time']['source'] = "telemetry"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delta_time']['calendar'] = "standard"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delta_time']['description'] = ("Number of "
            "GPS seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch "
            "offset is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS "
            "seconds between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP "
            "epoch. By adding the offset contained within atlas_sdp_gps_epoch to delta time "
            "parameters, the time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delta_time']['coordinates'] = \
            "latitude longitude"

        #-- latitude
        IS2_atl12_tide[gtx]['ssh_segments']['latitude'] = val['latitude']
        IS2_atl12_fill[gtx]['ssh_segments']['latitude'] = None
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude'] = {}
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude']['units'] = "degrees_north"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude']['long_name'] = "Latitude"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude']['standard_name'] = "latitude"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude']['description'] = ("Latitude of "
            "segment center")
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude']['valid_min'] = -90.0
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude']['valid_max'] = 90.0
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['latitude']['coordinates'] = \
            "delta_time longitude"
        #-- longitude
        IS2_atl12_tide[gtx]['ssh_segments']['longitude'] = val['longitude'].copy()
        IS2_atl12_fill[gtx]['ssh_segments']['longitude'] = None
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude'] = {}
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude']['units'] = "degrees_east"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude']['long_name'] = "Longitude"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude']['standard_name'] = "longitude"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude']['description'] = ("Longitude of "
            "segment center")
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude']['valid_min'] = -180.0
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude']['valid_max'] = 180.0
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['longitude']['coordinates'] = \
            "delta_time latitude"
        #-- Ocean Segment Duration
        IS2_atl12_tide[gtx]['ssh_segments']['delt_seg'] = val['delt_seg']
        IS2_atl12_fill[gtx]['ssh_segments']['delt_seg'] = None
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delt_seg'] = {}
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delt_seg']['units'] = "seconds"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delt_seg']['contentType'] = \
            "referenceInformation"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delt_seg']['long_name'] = \
            "Ocean Segment Duration"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delt_seg']['description'] = \
            "Time duration segment"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['delt_seg']['coordinates'] = \
            "delta_time latitude longitude"

        #-- stats variables
        IS2_atl12_tide[gtx]['ssh_segments']['stats'] = {}
        IS2_atl12_fill[gtx]['ssh_segments']['stats'] = {}
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'] = {}
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats']['Description'] = ("Contains parameters "
            "related to quality and corrections on the sea surface height parameters.")
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats']['data_rate'] = ("Data within this group "
            "are stored at the variable ocean processing segment rate.")

        #-- computed tide
        IS2_atl12_tide[gtx]['ssh_segments']['stats'][variable] = tide.copy()
        IS2_atl12_fill[gtx]['ssh_segments']['stats'][variable] = tide.fill_value
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'][variable] = {}
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'][variable]['units'] = "meters"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'][variable]['contentType'] = "referenceInformation"
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'][variable]['long_name'] = long_name
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'][variable]['description'] = description
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'][variable]['source'] = MODEL
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'][variable]['reference'] = reference
        IS2_atl12_tide_attrs[gtx]['ssh_segments']['stats'][variable]['coordinates'] = \
            "../../delta_time ../latitude ../longitude"

    #-- output tidal HDF5 file
    args = (PRD,MODEL,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX)
    file_format = '{0}_{1}_TIDES_{2}{3}{4}{5}{6}{7}_{8}{9}{10}_{11}_{12}{13}.h5'
    #-- print file information
    print('\t{0}'.format(file_format.format(*args))) if VERBOSE else None
    HDF5_atl12_tide_write(IS2_atl12_tide, IS2_atl12_tide_attrs, CLOBBER=True,
        INPUT=os.path.basename(FILE), FILL_VALUE=IS2_atl12_fill,
        FILENAME=os.path.join(DIRECTORY,file_format.format(*args)))
    #-- change the permissions mode
    os.chmod(os.path.join(DIRECTORY,file_format.format(*args)), MODE)

#-- PURPOSE: outputting the tide values for ICESat-2 data to HDF5
def HDF5_atl12_tide_write(IS2_atl12_tide, IS2_atl12_attrs, INPUT=None,
    FILENAME='', FILL_VALUE=None, CLOBBER=False):
    #-- setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    #-- open output HDF5 file
    fileID = h5py.File(os.path.expanduser(FILENAME), clobber)

    #-- create HDF5 records
    h5 = {}

    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    h5['ancillary_data'] = {}
    for k,v in IS2_atl12_tide['ancillary_data'].items():
        #-- Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl12_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- write each output beam
    beams = [k for k in IS2_atl12_tide.keys() if bool(re.match(r'gt\d[lr]',k))]
    for gtx in beams:
        fileID.create_group(gtx)
        #-- add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl12_attrs[gtx][att_name]
        #-- create ssh_segments group
        fileID[gtx].create_group('ssh_segments')
        h5[gtx] = dict(ssh_segments={})
        for att_name in ['Description','data_rate']:
            att_val = IS2_atl12_attrs[gtx]['ssh_segments'][att_name]
            fileID[gtx]['ssh_segments'].attrs[att_name] = att_val

        #-- delta_time
        v = IS2_atl12_tide[gtx]['ssh_segments']['delta_time']
        attrs = IS2_atl12_attrs[gtx]['ssh_segments']['delta_time']
        #-- Defining the HDF5 dataset variables
        val = '{0}/{1}/{2}'.format(gtx,'ssh_segments','delta_time')
        h5[gtx]['ssh_segments']['delta_time'] = fileID.create_dataset(val,
            np.shape(v), data=v, dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in attrs.items():
            h5[gtx]['ssh_segments']['delta_time'].attrs[att_name] = att_val

        #-- geolocation and segment description variables
        for k in ['latitude','longitude','delt_seg']:
            #-- values and attributes
            v = IS2_atl12_tide[gtx]['ssh_segments'][k]
            attrs = IS2_atl12_attrs[gtx]['ssh_segments'][k]
            fillvalue = FILL_VALUE[gtx]['ssh_segments'][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}'.format(gtx,'ssh_segments',k)
            h5[gtx]['ssh_segments'][k] = fileID.create_dataset(val,
                np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                compression='gzip')
            #-- attach dimensions
            for dim in ['delta_time']:
                h5[gtx]['ssh_segments'][k].dims.create_scale(
                    h5[gtx]['ssh_segments'][dim], dim)
                h5[gtx]['ssh_segments'][k].dims[0].attach_scale(
                    h5[gtx]['ssh_segments'][dim])
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['ssh_segments'][k].attrs[att_name] = att_val

        #-- add to stats variables
        key = 'stats'
        fileID[gtx]['ssh_segments'].create_group(key)
        h5[gtx]['ssh_segments'][key] = {}
        for att_name in ['Description','data_rate']:
            att_val=IS2_atl12_attrs[gtx]['ssh_segments'][key][att_name]
            fileID[gtx]['ssh_segments'][key].attrs[att_name] = att_val
        for k,v in IS2_atl12_tide[gtx]['ssh_segments'][key].items():
            #-- attributes
            attrs = IS2_atl12_attrs[gtx]['ssh_segments'][key][k]
            fillvalue = FILL_VALUE[gtx]['ssh_segments'][key][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}/{3}'.format(gtx,'ssh_segments',key,k)
            if fillvalue:
                h5[gtx]['ssh_segments'][key][k] = \
                    fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
            else:
                h5[gtx]['ssh_segments'][key][k] = \
                    fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, compression='gzip')
            #-- attach dimensions
            for dim in ['delta_time']:
                h5[gtx]['ssh_segments'][key][k].dims.create_scale(
                    h5[gtx]['ssh_segments'][dim], dim)
                h5[gtx]['ssh_segments'][key][k].dims[0].attach_scale(
                    h5[gtx]['ssh_segments'][dim])
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['ssh_segments'][key][k].attrs[att_name] = att_val

    #-- HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 L3A Ocean Surface Height'
    fileID.attrs['summary'] = ('Estimates of the ocean surface tidal parameters '
        'needed to interpret and assess the quality of ocean height estimates.')
    fileID.attrs['description'] = ('Sea Surface Height (SSH) of the global '
        'open ocean including the ice-free seasonal ice zone (SIZ) and '
        'near-coast regions.')
    date_created = datetime.datetime.today()
    fileID.attrs['date_created'] = date_created.isoformat()
    project = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = project
    platform = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = platform
    #-- add attribute for elevation instrument and designated processing level
    instrument = 'ATLAS > Advanced Topographic Laser Altimeter System'
    fileID.attrs['instrument'] = instrument
    fileID.attrs['source'] = 'Spacecraft'
    fileID.attrs['references'] = 'http://nsidc.org/data/icesat2/data.html'
    fileID.attrs['processing_level'] = '4'
    #-- add attributes for input ATL12 file
    fileID.attrs['input_files'] = os.path.basename(INPUT)
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        lon = IS2_atl12_tide[gtx]['ssh_segments']['longitude']
        lat = IS2_atl12_tide[gtx]['ssh_segments']['latitude']
        delta_time = IS2_atl12_tide[gtx]['ssh_segments']['delta_time']
        #-- setting the geospatial and temporal ranges
        lnmn = lon.min() if (lon.min() < lnmn) else lnmn
        lnmx = lon.max() if (lon.max() > lnmx) else lnmx
        ltmn = lat.min() if (lat.min() < ltmn) else ltmn
        ltmx = lat.max() if (lat.max() > ltmx) else ltmx
        tmn = delta_time.min() if (delta_time.min() < tmn) else tmn
        tmx = delta_time.max() if (delta_time.max() > tmx) else tmx
    #-- add geospatial and temporal attributes
    fileID.attrs['geospatial_lat_min'] = ltmn
    fileID.attrs['geospatial_lat_max'] = ltmx
    fileID.attrs['geospatial_lon_min'] = lnmn
    fileID.attrs['geospatial_lon_max'] = lnmx
    fileID.attrs['geospatial_lat_units'] = "degrees_north"
    fileID.attrs['geospatial_lon_units'] = "degrees_east"
    fileID.attrs['geospatial_ellipsoid'] = "WGS84"
    fileID.attrs['date_type'] = 'UTC'
    fileID.attrs['time_type'] = 'CCSDS UTC-A'
    #-- convert start and end time from ATLAS SDP seconds into GPS seconds
    atlas_sdp_gps_epoch=IS2_atl12_tide['ancillary_data']['atlas_sdp_gps_epoch']
    gps_seconds = atlas_sdp_gps_epoch + np.array([tmn,tmx])
    #-- calculate leap seconds
    leaps = pyTMD.time.count_leap_seconds(gps_seconds)
    #-- convert from seconds since 1980-01-06T00:00:00 to Julian days
    time_julian = 2400000.5 + pyTMD.time.convert_delta_time(gps_seconds - leaps,
        epoch1=(1980,1,6,0,0,0), epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    #-- convert to calendar date with convert_julian.py
    YY,MM,DD,HH,MN,SS = convert_julian(time_julian,FORMAT='tuple')
    #-- add attributes with measurement date start, end and duration
    tcs = datetime.datetime(np.int(YY[0]), np.int(MM[0]), np.int(DD[0]),
        np.int(HH[0]), np.int(MN[0]), np.int(SS[0]), np.int(1e6*(SS[0] % 1)))
    fileID.attrs['time_coverage_start'] = tcs.isoformat()
    tce = datetime.datetime(np.int(YY[1]), np.int(MM[1]), np.int(DD[1]),
        np.int(HH[1]), np.int(MN[1]), np.int(SS[1]), np.int(1e6*(SS[1] % 1)))
    fileID.attrs['time_coverage_end'] = tce.isoformat()
    fileID.attrs['time_coverage_duration'] = '{0:0.0f}'.format(tmx-tmn)
    #-- Closing the HDF5 file
    fileID.close()

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -T X, --tide=X\t\tTide model to use in correction')
    print(' -M X, --mode=X\t\tPermission mode of directories and files created')
    print(' -V, --verbose\t\tOutput information about each created file\n')

#-- Main program that calls compute_tides_ICESat2()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','tide=','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:T:VM:', long_options)

    #-- directory with tide data
    tide_dir = os.getcwd()
    #-- tide model to use
    MODEL = 'CATS2008'
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
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- enter HDF5 file as system argument
    if not arglist:
        raise Exception('No System Arguments Listed')

    #-- verify model before running program
    model_list = ['CATS0201','CATS2008','CATS2008_load','TPXO9-atlas','TPXO9.1',
        'TPXO8-atlas','TPXO7.2','TPXO7.2_load','AODTM-5','AOTIM-5',
        'AOTIM-5-2018','GOT4.7','GOT4.7_load','GOT4.8','GOT4.8_load',
        'GOT4.10','GOT4.10_load','FES2014','FES2014_load']
    assert MODEL in model_list, 'Unlisted tide model'

    #-- run for each input file
    for FILE in arglist:
        compute_tides_ICESat2(tide_dir, os.path.expanduser(FILE), MODEL,
            VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
