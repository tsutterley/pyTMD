#!/usr/bin/env python
u"""
compute_tides_ICESat_GLA12.py
Written by Tyler Sutterley (12/2020)
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
        CATS0201
        CATS2008
        CATS2008_load
        TPXO9-atlas
        TPXO9-atlas-v2
        TPXO9-atlas-v3
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
    -I X, --interpolate X: Interpolation method
        spline
        linear
        nearest
        bilinear
    -E X, --extrapolate X: Extrapolate with nearest-neighbors
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
import h5py
import argparse
import numpy as np
import pyTMD.time
import pyTMD.spatial
import pyTMD.utilities
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants
from pyTMD.read_FES_model import extract_FES_constants
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide_drift import predict_tide_drift

#-- PURPOSE: read ICESat ice sheet HDF5 elevation data (GLAH12) from NSIDC
#-- compute tides at points and times using tidal model driver algorithms
def compute_tides_ICESat(tide_dir, FILE, TIDE_MODEL=None, METHOD='spline',
     EXTRAPOLATE=False, VERBOSE=False, MODE=0o775):
    #-- select between tide models
    if (TIDE_MODEL == 'CATS0201'):
        grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
        model_file = os.path.join(tide_dir,'cats0201_tmd','h0_CATS02_01')
        reference = 'https://mail.esr.org/polar_tide_models/Model_CATS0201.html'
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (TIDE_MODEL == 'CATS2008'):
        grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008')
        model_file = os.path.join(tide_dir,'CATS2008','hf.CATS2008.out')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = 'CATS2008'
        TYPE = 'z'
    elif (TIDE_MODEL == 'CATS2008_load'):
        grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
        model_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        variable = 'd_ldElv'
        long_name = "Load Tide Elevation"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'OTIS'
        EPSG = 'CATS2008'
        TYPE = 'z'
    elif (TIDE_MODEL == 'TPXO9-atlas'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas')
        grid_file = 'grid_tpxo9_atlas.nc.gz'
        model_files = ['h_q1_tpxo9_atlas_30.nc.gz','h_o1_tpxo9_atlas_30.nc.gz',
            'h_p1_tpxo9_atlas_30.nc.gz','h_k1_tpxo9_atlas_30.nc.gz',
            'h_n2_tpxo9_atlas_30.nc.gz','h_m2_tpxo9_atlas_30.nc.gz',
            'h_s2_tpxo9_atlas_30.nc.gz','h_k2_tpxo9_atlas_30.nc.gz',
            'h_m4_tpxo9_atlas_30.nc.gz','h_ms4_tpxo9_atlas_30.nc.gz',
            'h_mn4_tpxo9_atlas_30.nc.gz','h_2n2_tpxo9_atlas_30.nc.gz']
        reference = 'http://volkov.oce.orst.edu/tides/tpxo9_atlas.html'
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'netcdf'
        TYPE = 'z'
        SCALE = 1.0/1000.0
    elif (TIDE_MODEL == 'TPXO9-atlas-v2'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas_v2')
        grid_file = 'grid_tpxo9_atlas_30_v2.nc.gz'
        model_files = ['h_q1_tpxo9_atlas_30_v2.nc.gz','h_o1_tpxo9_atlas_30_v2.nc.gz',
            'h_p1_tpxo9_atlas_30_v2.nc.gz','h_k1_tpxo9_atlas_30_v2.nc.gz',
            'h_n2_tpxo9_atlas_30_v2.nc.gz','h_m2_tpxo9_atlas_30_v2.nc.gz',
            'h_s2_tpxo9_atlas_30_v2.nc.gz','h_k2_tpxo9_atlas_30_v2.nc.gz',
            'h_m4_tpxo9_atlas_30_v2.nc.gz','h_ms4_tpxo9_atlas_30_v2.nc.gz',
            'h_mn4_tpxo9_atlas_30_v2.nc.gz','h_2n2_tpxo9_atlas_30_v2.nc.gz']
        reference = 'https://www.tpxo.net/global/tpxo9-atlas'
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'netcdf'
        TYPE = 'z'
        SCALE = 1.0/1000.0
    elif (TIDE_MODEL == 'TPXO9-atlas-v3'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas_v3')
        grid_file = 'grid_tpxo9_atlas_30_v3.nc.gz'
        model_files = ['h_q1_tpxo9_atlas_30_v3.nc.gz','h_o1_tpxo9_atlas_30_v3.nc.gz',
            'h_p1_tpxo9_atlas_30_v3.nc.gz','h_k1_tpxo9_atlas_30_v3.nc.gz',
            'h_n2_tpxo9_atlas_30_v3.nc.gz','h_m2_tpxo9_atlas_30_v3.nc.gz',
            'h_s2_tpxo9_atlas_30_v3.nc.gz','h_k2_tpxo9_atlas_30_v3.nc.gz',
            'h_m4_tpxo9_atlas_30_v3.nc.gz','h_ms4_tpxo9_atlas_30_v3.nc.gz',
            'h_mn4_tpxo9_atlas_30_v3.nc.gz','h_2n2_tpxo9_atlas_30_v3.nc.gz',
            'h_mf_tpxo9_atlas_30_v3.nc.gz','h_mm_tpxo9_atlas_30_v3.nc.gz']
        reference = 'https://www.tpxo.net/global/tpxo9-atlas'
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'netcdf'
        TYPE = 'z'
        SCALE = 1.0/1000.0
    elif (TIDE_MODEL == 'TPXO9.1'):
        grid_file = os.path.join(tide_dir,'TPXO9.1','DATA','grid_tpxo9')
        model_file = os.path.join(tide_dir,'TPXO9.1','DATA','h_tpxo9.v1')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (TIDE_MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_file = os.path.join(tide_dir,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
        reference = 'http://volkov.oce.orst.edu/tides/tpxo8_atlas.html'
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'ATLAS'
        EPSG = '4326'
        TYPE = 'z'
    elif (TIDE_MODEL == 'TPXO7.2'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_tmd','h_tpxo7.2')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (TIDE_MODEL == 'TPXO7.2_load'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_load','grid_tpxo6.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_load','h_tpxo7.2_load')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        variable = 'd_ldElv'
        long_name = "Load Tide Elevation"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'OTIS'
        EPSG = '4326'
        TYPE = 'z'
    elif (TIDE_MODEL == 'AODTM-5'):
        grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aodtm5_tmd','h0_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aodtm-5/')
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPE = 'z'
    elif (TIDE_MODEL == 'AOTIM-5'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPE = 'z'
    elif (TIDE_MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(tide_dir,'Arc5km2018','grid_Arc5km2018')
        model_file = os.path.join(tide_dir,'Arc5km2018','h_Arc5km2018')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPE = 'z'
    elif (TIDE_MODEL == 'GOT4.7'):
        model_directory = os.path.join(tide_dir,'GOT4.7','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (TIDE_MODEL == 'GOT4.7_load'):
        model_directory = os.path.join(tide_dir,'GOT4.7','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        variable = 'd_ldElv'
        long_name = "Load Tide Elevation"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (TIDE_MODEL == 'GOT4.8'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (TIDE_MODEL == 'GOT4.8_load'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        variable = 'd_ldElv'
        long_name = "Load Tide Elevation"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (TIDE_MODEL == 'GOT4.10'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (TIDE_MODEL == 'GOT4.10_load'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        variable = 'd_ldElv'
        long_name = "Load Tide Elevation"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (TIDE_MODEL == 'FES2014'):
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
        reference = ('https://www.aviso.altimetry.fr/en/data/products'
            'auxiliary-products/global-tide-fes.html')
        variable = 'd_ocElv'
        long_name = "Ocean Tide Elevation"
        description = ("Ocean Tides including diurnal and semi-diurnal "
            "(harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        model_format = 'FES'
        TYPE = 'z'
        SCALE = 1.0/100.0
    elif (TIDE_MODEL == 'FES2014_load'):
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
        reference = ('https://www.aviso.altimetry.fr/en/data/products'
            'auxiliary-products/global-tide-fes.html')
        variable = 'd_ldElv'
        long_name = "Load Tide Elevation"
        description = "Local displacement due to Ocean Loading (-6 to 0 cm)"
        model_format = 'FES'
        TYPE = 'z'
        SCALE = 1.0/100.0

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
    #-- read tidal constants and interpolate to grid points
    if model_format in ('OTIS','ATLAS'):
        amp,ph,D,c = extract_tidal_constants(lon_40HZ, lat_40HZ, grid_file,
            model_file, EPSG, TYPE=TYPE, METHOD=METHOD, EXTRAPOLATE=EXTRAPOLATE,
            GRID=model_format)
        deltat = np.zeros_like(tide_time)
    elif (model_format == 'netcdf'):
        amp,ph,D,c = extract_netcdf_constants(lon_40HZ, lat_40HZ, model_directory,
            grid_file, model_files, TYPE=TYPE, METHOD=METHOD, 
            EXTRAPOLATE=EXTRAPOLATE, SCALE=SCALE)
        deltat = np.zeros_like(tide_time)
    elif (model_format == 'GOT'):
        amp,ph = extract_GOT_constants(lon_40HZ, lat_40HZ, model_directory,
            model_files, METHOD=METHOD, EXTRAPOLATE=EXTRAPOLATE, SCALE=SCALE)
        #-- interpolate delta times from calendar dates to tide time
        delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
        deltat = calc_delta_time(delta_file, tide_time)
    elif (model_format == 'FES'):
        amp,ph = extract_FES_constants(lon_40HZ, lat_40HZ, model_directory,
            model_files, TYPE=TYPE, VERSION=TIDE_MODEL, METHOD=METHOD,
            EXTRAPOLATE=EXTRAPOLATE, SCALE=SCALE)
        #-- interpolate delta times from calendar dates to tide time
        delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
        deltat = calc_delta_time(delta_file, tide_time)

    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- calculate constituent oscillation
    hc = amp*np.exp(cph)

    #-- predict tidal elevations at time and infer minor corrections
    tide = np.ma.empty((n_40HZ),fill_value=fv)
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
    IS_gla12_tide_attrs['Data_40HZ']['Time']['i_rec_ndx']['coordinates'] = \
        "../DS_UTCTime_40"
    for att_name,att_val in fileID['Data_40HZ']['Time']['i_rec_ndx'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Time']['i_rec_ndx'][att_name] = att_val
    #-- latitude
    IS_gla12_tide['Data_40HZ']['Geolocation']['d_lat'] = lat_40HZ
    IS_gla12_fill['Data_40HZ']['Geolocation']['d_lat'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lat'] = {}
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lat']['coordinates'] = \
        "../DS_UTCTime_40"
    for att_name,att_val in fileID['Data_40HZ']['Geolocation']['d_lat'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lat'][att_name] = att_val
    #-- longitude
    IS_gla12_tide['Data_40HZ']['Geolocation']['d_lon'] = lon_40HZ
    IS_gla12_fill['Data_40HZ']['Geolocation']['d_lon'] = None
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lon'] = {}
    IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lon']['coordinates'] = \
        "../DS_UTCTime_40"
    for att_name,att_val in fileID['Data_40HZ']['Geolocation']['d_lon'].attrs.items():
        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
            IS_gla12_tide_attrs['Data_40HZ']['Geolocation']['d_lon'][att_name] = att_val

    #-- geophysical variables
    #-- computed tide
    IS_gla12_tide['Data_40HZ']['Geophysical'][variable] = tide
    IS_gla12_fill['Data_40HZ']['Geophysical'][variable] = tide.fill_value
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][variable] = {}
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][variable]['units'] = "meters"
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][variable]['long_name'] = long_name
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][variable]['description'] = description
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][variable]['source'] = TIDE_MODEL
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][variable]['reference'] = reference
    IS_gla12_tide_attrs['Data_40HZ']['Geophysical'][variable]['coordinates'] = \
        "../DS_UTCTime_40"

    #-- close the input HDF5 file
    fileID.close()

    #-- output tidal HDF5 file
    args = (PRD,RL,TIDE_MODEL,RGTP,ORB,INST,CYCL,TRK,SEG,GRAN,TYPE)
    file_format = 'GLAH{0}_{1}_{2}_TIDES_{3}{4}{5}_{6}_{7}_{8}_{9}_{10}.h5'
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

#-- Main program that calls compute_tides_ICESat()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Calculates tidal elevations for correcting ICESat/GLAS
            L2 GLA12 Antarctic and Greenland Ice Sheet elevation data
            """
    )
    #-- command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='ICESat GLA12 file to run')
    #-- directory with tide data
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- tide model to use
    model_choices = ('CATS0201','CATS2008','CATS2008_load',
        'TPXO9-atlas','TPXO9-atlas-v2','TPXO9-atlas-v3',
        'TPXO9.1','TPXO8-atlas','TPXO7.2','TPXO7.2_load',
        'AODTM-5','AOTIM-5','AOTIM-5-2018',
        'GOT4.7','GOT4.7_load','GOT4.8','GOT4.8_load','GOT4.10','GOT4.10_load',
        'FES2014','FES2014_load')
    parser.add_argument('--tide','-T',
        metavar='TIDE', type=str, default='CATS2008',
        choices=model_choices,
        help='Tide model to use in correction')
    #-- interpolation method
    parser.add_argument('--interpolate','-I',
        metavar='METHOD', type=str, default='spline',
        choices=('spline','linear','nearest','bilinear'),
        help='Spatial interpolation method')
    #-- extrapolate with nearest-neighbors
    parser.add_argument('--extrapolate','-E',
        default=False, action='store_true',
        help='Extrapolate with nearest-neighbors')
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
        compute_tides_ICESat(args.directory, FILE, TIDE_MODEL=args.tide,
            METHOD=args.interpolate, EXTRAPOLATE=args.extrapolate,
            VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
