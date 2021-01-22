#!/usr/bin/env python
u"""
compute_tidal_elevations.py
Written by Tyler Sutterley (12/2020)
Calculates tidal elevations for an input file

Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
or Finite Element Solution (FES) models provided by AVISO

INPUTS:
    csv file with columns for spatial and temporal coordinates
    HDF5 file with variables for spatial and temporal coordinates
    netCDF4 file with variables for spatial and temporal coordinates
    geotiff file with bands in spatial coordinates

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
    --format X: input and output data format
        csv (default)
        netCDF4
        HDF5
        geotiff
    --variables X: variable names of data in csv, HDF5 or netCDF4 file
        for csv files: the order of the columns within the file
        for HDF5 and netCDF4 files: time, y, x and data variable names
    -H X, --header X: number of header lines for csv files
    -t X, --type X: input data type
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time for all data points)
    -e X, --epoch X: Reference epoch of input time (default Modified Julian Day)
        days since 1858-11-17T00:00:00
    -d X, --deltatime X: Input delta time for files without date information
        can be set to 0 to use exact calendar date from epoch
    -P X, --projection X: spatial projection as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
    -I X, --interpolate X: Interpolation method
        spline
        linear
        nearest
        bilinear
    -E X, --extrapolate X: Extrapolate with nearest-neighbors
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permission mode of output file

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
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    spatial.py: utilities for reading and writing spatial data
    utilities: download and management utilities for syncing files
    calc_astrol_longitudes.py: computes the basic astronomical mean longitudes
    calc_delta_time.py: calculates difference between universal and dynamic time
    convert_ll_xy.py: convert lat/lon points to and from projected coordinates
    load_constituent.py: loads parameters for a given tidal constituent
    load_nodal_corrections.py: load the nodal corrections for tidal constituents
    infer_minor_corrections.py: return corrections for minor constituents
    read_tide_model.py: extract tidal harmonic constants from OTIS tide models
    read_netcdf_model.py: extract tidal harmonic constants from netcdf models
    read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
    read_FES_model.py: extract tidal harmonic constants from FES tide models
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates
    predict_tide_drift.py: predict tidal elevations using harmonic constants

UPDATE HISTORY:
    Updated 12/2020: added valid data extrapolation with nearest_extrap
    Updated 11/2020: added model constituents from TPXO9-atlas-v3
        added options to read from and write to geotiff image files
    Updated 10/2020: using argparse to set command line parameters
    Updated 09/2020: can use HDF5 and netCDF4 as inputs and outputs
    Updated 08/2020: using builtin time operations
    Updated 07/2020: added FES2014 and FES2014_load.  use merged delta times
    Updated 06/2020: added version 2 of TPXO9-atlas (TPXO9-atlas-v2)
    Updated 02/2020: changed CATS2008 grid to match version on U.S. Antarctic
        Program Data Center http://www.usap-dc.org/view/dataset/601235
    Updated 11/2019: added AOTIM-5-2018 tide model (2018 update to 2004 model)
    Updated 09/2019: added TPXO9_atlas reading from netcdf4 tide files
    Updated 07/2018: added GSFC Global Ocean Tides (GOT) models
    Written 10/2017 for public release
"""
from __future__ import print_function

import sys
import os
import pyproj
import argparse
import numpy as np
import pyTMD.time
import pyTMD.spatial
from pyTMD.utilities import get_data_path
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide import predict_tide
from pyTMD.predict_tide_drift import predict_tide_drift
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants
from pyTMD.read_FES_model import extract_FES_constants

#-- PURPOSE: read csv, netCDF or HDF5 data
#-- compute tides at points and times using tidal model driver algorithms
def compute_tidal_elevations(tide_dir, input_file, output_file,
    TIDE_MODEL=None, FORMAT='csv', VARIABLES=['time','lat','lon','data'],
    HEADER=0, TYPE='drift', TIME_UNITS='days since 1858-11-17T00:00:00',
    TIME=None, PROJECTION='4326', METHOD='spline', EXTRAPOLATE=False,
    VERBOSE=False, MODE=0o775):

    #-- select between tide models
    if (TIDE_MODEL == 'CATS0201'):
        grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
        model_file = os.path.join(tide_dir,'cats0201_tmd','h0_CATS02_01')
        reference = 'https://mail.esr.org/polar_tide_models/Model_CATS0201.html'
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        model_type = 'z'
    elif (TIDE_MODEL == 'CATS2008'):
        grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008')
        model_file = os.path.join(tide_dir,'CATS2008','hf.CATS2008.out')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = 'CATS2008'
        model_type = 'z'
    elif (TIDE_MODEL == 'CATS2008_load'):
        grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
        model_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        output_variable = 'tide_load'
        variable_long_name = 'Load_Tide'
        model_format = 'OTIS'
        EPSG = 'CATS2008'
        model_type = 'z'
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
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'netcdf'
        model_type = 'z'
        model_scale = 1.0/1000.0
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
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'netcdf'
        model_type = 'z'
        model_scale = 1.0/1000.0
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
        variable = 'tide_ocean'
        long_name = "Ocean_Tide"
        model_format = 'netcdf'
        model_type = 'z'
        model_scale = 1.0/1000.0
    elif (TIDE_MODEL == 'TPXO9.1'):
        grid_file = os.path.join(tide_dir,'TPXO9.1','DATA','grid_tpxo9')
        model_file = os.path.join(tide_dir,'TPXO9.1','DATA','h_tpxo9.v1')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        model_type = 'z'
    elif (TIDE_MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_file = os.path.join(tide_dir,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
        reference = 'http://volkov.oce.orst.edu/tides/tpxo8_atlas.html'
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'ATLAS'
        EPSG = '4326'
        model_type = 'z'
    elif (TIDE_MODEL == 'TPXO7.2'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_tmd','h_tpxo7.2')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        model_type = 'z'
    elif (TIDE_MODEL == 'TPXO7.2_load'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_load','grid_tpxo6.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_load','h_tpxo7.2_load')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        output_variable = 'tide_load'
        variable_long_name = 'Load_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        model_type = 'z'
    elif (TIDE_MODEL == 'AODTM-5'):
        grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aodtm5_tmd','h0_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aodtm-5/')
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        model_type = 'z'
    elif (TIDE_MODEL == 'AOTIM-5'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        model_type = 'z'
    elif (TIDE_MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(tide_dir,'Arc5km2018','grid_Arc5km2018')
        model_file = os.path.join(tide_dir,'Arc5km2018','h_Arc5km2018')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        model_type = 'z'
    elif (TIDE_MODEL == 'GOT4.7'):
        model_directory = os.path.join(tide_dir,'GOT4.7','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'GOT'
        model_scale = 1.0/100.0
    elif (TIDE_MODEL == 'GOT4.7_load'):
        model_directory = os.path.join(tide_dir,'GOT4.7','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        output_variable = 'tide_load'
        variable_long_name = 'Load_Tide'
        model_format = 'GOT'
        model_scale = 1.0/1000.0
    elif (TIDE_MODEL == 'GOT4.8'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'GOT'
        model_scale = 1.0/100.0
    elif (TIDE_MODEL == 'GOT4.8_load'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        output_variable = 'tide_load'
        variable_long_name = 'Load_Tide'
        model_format = 'GOT'
        model_scale = 1.0/1000.0
    elif (TIDE_MODEL == 'GOT4.10'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'GOT'
        model_scale = 1.0/100.0
    elif (TIDE_MODEL == 'GOT4.10_load'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
            'MiscPubs/19990089548_1999150788.pdf')
        output_variable = 'tide_load'
        variable_long_name = 'Load_Tide'
        model_format = 'GOT'
        model_scale = 1.0/1000.0
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
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'FES'
        model_type = 'z'
        model_scale = 1.0/100.0
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
        output_variable = 'tide_load'
        variable_long_name = 'Load_Tide'
        model_format = 'FES'
        model_type = 'z'
        model_scale = 1.0/100.0

    #-- invalid value
    fill_value = -9999.0
    #-- output netCDF4 and HDF5 file attributes
    #-- will be added to YAML header in csv files
    attrib = {}
    #-- latitude
    attrib['lat'] = {}
    attrib['lat']['long_name'] = 'Latitude'
    attrib['lat']['units'] = 'Degrees_North'
    #-- longitude
    attrib['lon'] = {}
    attrib['lon']['long_name'] = 'Longitude'
    attrib['lon']['units'] = 'Degrees_East'
    #-- tides
    attrib[output_variable] = {}
    attrib[output_variable]['description'] = ('tidal_elevation_from_harmonic_'
        'constants')
    attrib[output_variable]['model'] = TIDE_MODEL
    attrib[output_variable]['units'] = 'meters'
    attrib[output_variable]['long_name'] = variable_long_name
    attrib[output_variable]['_FillValue'] = fill_value
    #-- time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'days since 1992-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    #-- read input file to extract time, spatial coordinates and data
    if (FORMAT == 'csv'):
        dinput = pyTMD.spatial.from_ascii(input_file, columns=VARIABLES,
            header=HEADER, verbose=VERBOSE)
    elif (FORMAT == 'netCDF4'):
        dinput = pyTMD.spatial.from_netCDF4(input_file, timename=VARIABLES[0],
            xname=VARIABLES[2], yname=VARIABLES[1], varname=VARIABLES[3],
            verbose=VERBOSE)
    elif (FORMAT == 'HDF5'):
        dinput = pyTMD.spatial.from_HDF5(input_file, timename=VARIABLES[0],
            xname=VARIABLES[2], yname=VARIABLES[1], varname=VARIABLES[3],
            verbose=VERBOSE)
    elif (FORMAT == 'geotiff'):
        dinput = pyTMD.spatial.from_geotiff(input_file, verbose=VERBOSE)
        #-- copy global geotiff attributes for projection and grid parameters
        for att_name in ['projection','wkt','spacing','extent']:
            attrib[att_name] = dinput['attributes'][att_name]
    #-- update time variable if entered as argument
    if TIME is not None:
        dinput['time'] = np.copy(TIME)

    #-- converting x,y from projection to latitude/longitude
    #-- could try to extract projection attributes from netCDF4 and HDF5 files
    try:
        crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(int(PROJECTION)))
    except (ValueError,pyproj.exceptions.CRSError):
        crs1 = pyproj.CRS.from_string(PROJECTION)
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    if (TYPE == 'grid'):
        ny,nx = (len(dinput['y']),len(dinput['x']))
        gridx,gridy = np.meshgrid(dinput['x'],dinput['y'])
        lon,lat = transformer.transform(gridx,gridy)
    elif (TYPE == 'drift'):
        lon,lat = transformer.transform(dinput['x'],dinput['y'])

    #-- extract time units from netCDF4 and HDF5 attributes or from TIME_UNITS
    try:
        time_string = dinput['attributes']['time']['units']
    except (TypeError, KeyError):
        epoch1,to_secs = pyTMD.time.parse_date_string(TIME_UNITS)
    else:
        epoch1,to_secs = pyTMD.time.parse_date_string(time_string)
    #-- convert time from units to days since 1992-01-01T00:00:00
    tide_time = pyTMD.time.convert_delta_time(to_secs*dinput['time'].flatten(),
        epoch1=epoch1, epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
    #-- number of time points
    nt = len(tide_time)

    #-- read tidal constants and interpolate to grid points
    if model_format in ('OTIS','ATLAS'):
        amp,ph,D,c = extract_tidal_constants(lon.flatten(), lat.flatten(),
            grid_file, model_file, EPSG, TYPE=model_type, METHOD=METHOD,
            EXTRAPOLATE=EXTRAPOLATE)
        deltat = np.zeros((nt))
    elif (model_format == 'netcdf'):
        amp,ph,D,c = extract_netcdf_constants(lon.flatten(), lat.flatten(),
            model_directory, grid_file, model_files, TYPE=model_type,
            METHOD=METHOD, EXTRAPOLATE=EXTRAPOLATE, SCALE=model_scale)
        deltat = np.zeros((nt))
    elif (model_format == 'GOT'):
        amp,ph = extract_GOT_constants(lon.flatten(), lat.flatten(),
            model_directory, model_files, METHOD=METHOD,
            EXTRAPOLATE=EXTRAPOLATE, SCALE=model_scale)
        #-- convert times from modified julian days to days since 1992-01-01
        #-- interpolate delta times from calendar dates to tide time
        delta_file = get_data_path(['data','merged_deltat.data'])
        deltat = calc_delta_time(delta_file,tide_time)
    elif (model_format == 'FES'):
        amp,ph = extract_FES_constants(lon.flatten(), lat.flatten(),
            model_directory, model_files, TYPE=model_type, VERSION=TIDE_MODEL,
            METHOD=METHOD, EXTRAPOLATE=EXTRAPOLATE, SCALE=model_scale)
        #-- convert times from modified julian days to days since 1992-01-01
        #-- interpolate delta times from calendar dates to tide time
        delta_file = get_data_path(['data','merged_deltat.data'])
        deltat = calc_delta_time(delta_file,tide_time)

    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- calculate constituent oscillation
    hc = amp*np.exp(cph)

    #-- predict tidal elevations at time and infer minor corrections
    if (TYPE == 'grid'):
        tide = np.ma.zeros((ny,nx,nt),fill_value=fill_value)
        tide.mask = np.zeros((ny,nx,nt),dtype=np.bool)
        for i in range(nt):
            TIDE = predict_tide(tide_time[i], hc, c,
                DELTAT=deltat[i], CORRECTIONS=model_format)
            MINOR = infer_minor_corrections(tide_time[i], hc, c,
                DELTAT=deltat[i], CORRECTIONS=model_format)
            #-- add major and minor components and reform grid
            tide[:,:,i] = np.reshape((TIDE+MINOR), (ny,nx))
            tide.mask[:,:,i] = np.reshape((TIDE.mask | MINOR.mask), (ny,nx))
    elif (TYPE == 'drift'):
        tide = np.ma.zeros((nt), fill_value=fill_value)
        tide.mask = np.any(hc.mask,axis=1)
        tide.data[:] = predict_tide_drift(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=model_format)
        minor = infer_minor_corrections(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=model_format)
        tide.data[:] += minor.data[:]
    #-- replace invalid values with fill value
    tide.data[tide.mask] = tide.fill_value

    #-- output to file
    output = {'time':tide_time,'lon':lon,'lat':lat,output_variable:tide}
    if (FORMAT == 'csv'):
        pyTMD.spatial.to_ascii(output, attrib, output_file, delimiter=',',
            columns=['time','lat','lon',output_variable], verbose=VERBOSE)
    elif (FORMAT == 'netCDF4'):
        pyTMD.spatial.to_netCDF4(output, attrib, output_file, verbose=VERBOSE)
    elif (FORMAT == 'HDF5'):
        pyTMD.spatial.to_HDF5(output, attrib, output_file, verbose=VERBOSE)
    elif (FORMAT == 'geotiff'):
        pyTMD.spatial.to_geotiff(output, attrib, output_file, verbose=VERBOSE,
            varname=output_variable)
    #-- change the permissions level to MODE
    os.chmod(output_file, MODE)

#-- Main program that calls compute_tidal_elevations()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Calculates tidal elevations for an input file
            """
    )
    #-- command line options
    #-- input and output file
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Input file')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Output file')
    #-- set data directory containing the tidal data
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
    #-- input and output data format
    parser.add_argument('--format','-F',
        type=str, default='csv', choices=('csv','netCDF4','HDF5','geotiff'),
        help='Input and output data format')
    #-- variable names (for csv names of columns)
    parser.add_argument('--variables','-v',
        type=str, nargs='+', default=['time','lat','lon','data'],
        help='Variable names of data in input file')
    #-- number of header lines for csv files
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of header lines for csv files')
    #-- input data type
    #-- drift: drift buoys or satellite/airborne altimetry (time per data point)
    #-- grid: spatial grids or images (single time for all data points)
    parser.add_argument('--type','-t',
        type=str, default='drift',
        choices=('drift','grid'),
        help='Input data type')
    #-- time epoch (default Modified Julian Days)
    #-- in form "time-units since yyyy-mm-dd hh:mm:ss"
    parser.add_argument('--epoch','-e',
        type=str, default='days since 1858-11-17T00:00:00',
        help='Reference epoch of input time')
    #-- input delta time for files without date information
    parser.add_argument('--deltatime','-d',
        type=float, nargs='+',
        help='Input delta time for files without date variables')
    #-- spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection','-P',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
    #-- interpolation method
    parser.add_argument('--interpolate','-I',
        metavar='METHOD', type=str, default='spline',
        choices=('spline','linear','nearest','bilinear'),
        help='Spatial interpolation method')
    #-- extrapolate with nearest-neighbors
    parser.add_argument('--extrapolate','-E',
        default=False, action='store_true',
        help='Extrapolate with nearest-neighbors')
    #-- verbose output of processing run
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of output file')
    args = parser.parse_args()

    #-- set output file from input filename if not entered
    if not args.outfile:
        fileBasename,fileExtension = os.path.splitext(args.infile)
        vars = (fileBasename,args.tide,fileExtension)
        args.outfile = '{0}_{1}{2}'.format(*vars)

    #-- run tidal elevation program for input file
    compute_tidal_elevations(args.directory, args.infile, args.outfile,
        FORMAT=args.format, TIDE_MODEL=args.tide, VARIABLES=args.variables,
        HEADER=args.header, TYPE=args.type, TIME_UNITS=args.epoch,
        TIME=args.deltatime, PROJECTION=args.projection,
        METHOD=args.interpolate, EXTRAPOLATE=args.extrapolate,
        VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
