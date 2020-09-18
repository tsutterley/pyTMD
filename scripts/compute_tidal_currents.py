#!/usr/bin/env python
u"""
compute_tidal_currents.py
Written by Tyler Sutterley (09/2020)
Calculates tidal currents for an input file

Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
or Finite Element Solution (FES) models provided by AVISO

INPUTS:
    csv file with columns for spatial and temporal coordinates
    HDF5 file with variables for spatial and temporal coordinates
    netCDF4 file with variables for spatial and temporal coordinates

COMMAND LINE OPTIONS:
    -D X, --directory=X: Working data directory
    -T X, --tide=X: Tide model to use in correction
        CATS0201
        CATS2008
        TPX09-atlas-v2
        TPXO9-atlas
        TPXO9.1
        TPXO8-atlas
        TPXO7.2
        TPXO7.2_load
        AODTM-5
        AOTIM-5
        AOTIM-5-2018
        FES2014
    --format=X: input and output data format
        csv (default)
        netCDF4
        HDF5
    --variables=X: variable names of data in csv, HDF5 or netCDF4 file
        for csv files: the order of the columns within the file
        for HDF5 and netCDF4 files: time, y, x and data variable names
    --epoch=X: Reference epoch of input time (default Modified Julian Day)
        days since 1858-11-17T00:00:00
    --projection=X: spatial projection as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
    -V, --verbose: Verbose output of processing run
    -M X, --mode=X: Permission mode of output file

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
    read_FES_model.py: extract tidal harmonic constants from FES tide models
    predict_tide_drift.py: predict tidal elevations using harmonic constants

UPDATE HISTORY:
    Forked 09/2020 from compute_tidal_elevations.py
    Updated 09/2020: can use HDF5 and netCDF4 as inputs and outputs
    Updated 08/2020: using builtin time operations
    Updated 07/2020: added FES2014 and FES2014_load.  use merged delta times
    Updated 06/2020: added version 2 of TPX09-atlas (TPX09-atlas-v2)
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
import getopt
import pyproj
import numpy as np
import pyTMD.time
import pyTMD.spatial
from pyTMD.utilities import get_data_path
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide_drift import predict_tide_drift
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_FES_model import extract_FES_constants

#-- PURPOSE: read csv, netCDF or HDF5 data
#-- compute tides at points and times using tidal model driver algorithms
def compute_tidal_currents(tide_dir, input_file, output_file,
    TIDE_MODEL='', FORMAT='csv', VARIABLES=['time','lat','lon','data'],
    TIME_UNITS='days since 1858-11-17T00:00:00', PROJECTION='4326',
    VERBOSE=False, MODE=0o775):

    #-- select between tide models
    if (TIDE_MODEL == 'CATS0201'):
        grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
        model_file = os.path.join(tide_dir,'cats0201_tmd','UV0_CATS02_01')
        reference = 'https://mail.esr.org/polar_tide_models/Model_CATS0201.html'
        model_format = 'OTIS'
        EPSG = '4326'
        TYPES = ['u','v']
    elif (TIDE_MODEL == 'CATS2008'):
        grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008')
        model_file = os.path.join(tide_dir,'CATS2008','uv.CATS2008.out')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        model_format = 'OTIS'
        EPSG = 'CATS2008'
        TYPES = ['u','v']
    elif (TIDE_MODEL == 'TPXO9-atlas'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas')
        grid_file = 'grid_tpxo9_atlas.nc.gz'
        model_files = {}
        model_files['u'] = ['u_q1_tpxo9_atlas_30.nc.gz','u_o1_tpxo9_atlas_30.nc.gz',
            'u_p1_tpxo9_atlas_30.nc.gz','u_k1_tpxo9_atlas_30.nc.gz',
            'u_n2_tpxo9_atlas_30.nc.gz','u_m2_tpxo9_atlas_30.nc.gz',
            'u_s2_tpxo9_atlas_30.nc.gz','u_k2_tpxo9_atlas_30.nc.gz',
            'u_m4_tpxo9_atlas_30.nc.gz','u_ms4_tpxo9_atlas_30.nc.gz',
            'u_mn4_tpxo9_atlas_30.nc.gz','u_2n2_tpxo9_atlas_30.nc.gz']
        model_files['v'] = ['v_q1_tpxo9_atlas_30.nc.gz','v_o1_tpxo9_atlas_30.nc.gz',
            'v_p1_tpxo9_atlas_30.nc.gz','v_k1_tpxo9_atlas_30.nc.gz',
            'v_n2_tpxo9_atlas_30.nc.gz','v_m2_tpxo9_atlas_30.nc.gz',
            'v_s2_tpxo9_atlas_30.nc.gz','v_k2_tpxo9_atlas_30.nc.gz',
            'v_m4_tpxo9_atlas_30.nc.gz','v_ms4_tpxo9_atlas_30.nc.gz',
            'v_mn4_tpxo9_atlas_30.nc.gz','v_2n2_tpxo9_atlas_30.nc.gz']
        reference = 'http://volkov.oce.orst.edu/tides/tpxo9_atlas.html'
        model_format = 'netcdf'
        TYPES = ['u','v']
        SCALE = 1.0/100.0
        GZIP = True
    elif (TIDE_MODEL == 'TPXO9-atlas-v2'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas_v2')
        grid_file = 'grid_tpxo9_atlas_30_v2.nc.gz'
        model_files = {}
        model_files['u'] = ['u_q1_tpxo9_atlas_30_v2.nc.gz','u_o1_tpxo9_atlas_30_v2.nc.gz',
            'u_p1_tpxo9_atlas_30_v2.nc.gz','u_k1_tpxo9_atlas_30_v2.nc.gz',
            'u_n2_tpxo9_atlas_30_v2.nc.gz','u_m2_tpxo9_atlas_30_v2.nc.gz',
            'u_s2_tpxo9_atlas_30_v2.nc.gz','u_k2_tpxo9_atlas_30_v2.nc.gz',
            'u_m4_tpxo9_atlas_30_v2.nc.gz','u_ms4_tpxo9_atlas_30_v2.nc.gz',
            'u_mn4_tpxo9_atlas_30_v2.nc.gz','u_2n2_tpxo9_atlas_30_v2.nc.gz']
        model_files['v'] = ['v_q1_tpxo9_atlas_30_v2.nc.gz','v_o1_tpxo9_atlas_30_v2.nc.gz',
            'v_p1_tpxo9_atlas_30_v2.nc.gz','v_k1_tpxo9_atlas_30_v2.nc.gz',
            'v_n2_tpxo9_atlas_30_v2.nc.gz','v_m2_tpxo9_atlas_30_v2.nc.gz',
            'v_s2_tpxo9_atlas_30_v2.nc.gz','v_k2_tpxo9_atlas_30_v2.nc.gz',
            'v_m4_tpxo9_atlas_30_v2.nc.gz','v_ms4_tpxo9_atlas_30_v2.nc.gz',
            'v_mn4_tpxo9_atlas_30_v2.nc.gz','v_2n2_tpxo9_atlas_30_v2.nc.gz']
        reference = 'https://www.tpxo.net/global/tpxo9-atlas'
        model_format = 'netcdf'
        TYPES = ['u','v']
        SCALE = 1.0/100.0
        GZIP = True
    elif (TIDE_MODEL == 'TPXO9.1'):
        grid_file = os.path.join(tide_dir,'TPXO9.1','DATA','grid_tpxo9')
        model_file = os.path.join(tide_dir,'TPXO9.1','DATA','u_tpxo9.v1')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        model_format = 'OTIS'
        EPSG = '4326'
        TYPES = ['u','v']
    elif (TIDE_MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_file = os.path.join(tide_dir,'tpxo8_atlas','uv.tpxo8_atlas_30_v1')
        reference = 'http://volkov.oce.orst.edu/tides/tpxo8_atlas.html'
        model_format = 'ATLAS'
        EPSG = '4326'
        TYPES = ['u','v']
    elif (TIDE_MODEL == 'TPXO7.2'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_tmd','u_tpxo7.2')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        output_variable = 'tide_ocean'
        variable_long_name = 'Ocean_Tide'
        model_format = 'OTIS'
        EPSG = '4326'
        TYPES = ['u','v']
    elif (TIDE_MODEL == 'AODTM-5'):
        grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aodtm5_tmd','UV0_Arc5km')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aodtm-5/')
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPES = ['u','v']
    elif (TIDE_MODEL == 'AOTIM-5'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aotim5_tmd','UV_Arc5km')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPES = ['u','v']
    elif (TIDE_MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(tide_dir,'Arc5km2018','grid_Arc5km2018')
        model_file = os.path.join(tide_dir,'Arc5km2018','UV_Arc5km2018')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        model_format = 'OTIS'
        EPSG = 'PSNorth'
        TYPES = ['u','v']
    elif (TIDE_MODEL == 'FES2014'):
        model_directory = {}
        model_directory['u'] = os.path.join(tide_dir,'fes2014','eastward_velocity')
        model_directory['v'] = os.path.join(tide_dir,'fes2014','northward_velocity')
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
        model_format = 'FES'
        TYPES = ['u','v']
        SCALE = 1.0/100.0

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
    #-- zonal tidal currents
    attrib['u'] = {}
    attrib['u']['description'] = ('depth_averaged_tidal_zonal_current_'
        'from_harmonic_constants')
    attrib['u']['model'] = TIDE_MODEL
    attrib['u']['units'] = 'cm/s'
    attrib['u']['long_name'] = 'zonal_tidal_current'
    attrib['u']['_FillValue'] = fill_value
    #-- meridional tidal currents
    attrib['v'] = {}
    attrib['v']['description'] = ('depth_averaged_tidal_meridional_current_'
        'from_harmonic_constants')
    attrib['v']['model'] = TIDE_MODEL
    attrib['v']['units'] = 'cm/s'
    attrib['v']['long_name'] = 'meridional_tidal_current'
    attrib['v']['_FillValue'] = fill_value
    #-- time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'days since 1992-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    #-- read input file to extract time, spatial coordinates and data
    if (FORMAT == 'csv'):
        dinput = pyTMD.spatial.from_ascii(input_file, columns=VARIABLES,
            header=0, verbose=VERBOSE)
    elif (FORMAT == 'netCDF4'):
        dinput = pyTMD.spatial.from_netCDF4(input_file, timename=VARIABLES[0],
            xname=VARIABLES[2], yname=VARIABLES[1], varname=VARIABLES[3],
            verbose=VERBOSE)
    elif (FORMAT == 'HDF5'):
        dinput = pyTMD.spatial.from_HDF5(input_file, timename=VARIABLES[0],
            xname=VARIABLES[2], yname=VARIABLES[1], varname=VARIABLES[3],
            verbose=VERBOSE)

    #-- converting x,y from projection to latitude/longitude
    #-- could try to extract projection attributes from netCDF4 and HDF5 files
    try:
        crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(int(PROJECTION)))
    except (ValueError,pyproj.exceptions.CRSError):
        crs1 = pyproj.CRS.from_string(PROJECTION)
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(dinput['x'].flatten(),dinput['y'].flatten())

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
    n_time = len(tide_time)

    #-- python dictionary with output data
    output = {'time':tide_time,'lon':lon,'lat':lat}
    #-- iterate over u and v currents
    for TYPE in TYPES:
        #-- read tidal constants and interpolate to grid points
        if model_format in ('OTIS','ATLAS'):
            amp,ph,D,c = extract_tidal_constants(lon, lat, grid_file,
                model_file, EPSG, TYPE=TYPE, METHOD='spline', GRID=model_format)
            deltat = np.zeros_like(tide_time)
        elif (model_format == 'netcdf'):
            amp,ph,D,c = extract_netcdf_constants(lon, lat, model_directory,
                grid_file, model_files[TYPE], TYPE=TYPE, METHOD='spline',
                SCALE=SCALE, GZIP=GZIP)
            deltat = np.zeros_like(tide_time)
        elif (model_format == 'GOT'):
            amp,ph = extract_GOT_constants(lon, lat, model_directory,
                model_files, METHOD='spline', SCALE=SCALE)
            #-- interpolate delta times from calendar dates to tide time
            delta_file = get_data_path(['data','merged_deltat.data'])
            deltat = calc_delta_time(delta_file, tide_time)
        elif (model_format == 'FES'):
            amp,ph = extract_FES_constants(lon, lat, model_directory[TYPE],
                model_files, TYPE=TYPE, VERSION=MODEL, METHOD=METHOD,
                SCALE=SCALE)
            #-- interpolate delta times from calendar dates to tide time
            delta_file = get_data_path(['data','merged_deltat.data'])
            deltat = calc_delta_time(delta_file, tide_time)

        #-- calculate complex phase in radians for Euler's
        cph = -1j*ph*np.pi/180.0
        #-- calculate constituent oscillation
        hc = amp*np.exp(cph)

        #-- predict tidal currents at time and infer minor corrections
        output[TYPE] = np.ma.zeros((n_time), fill_value=fill_value)
        output[TYPE].mask = np.any(hc.mask,axis=1)
        output[TYPE].data[:] = predict_tide_drift(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=model_format)
        minor = infer_minor_corrections(tide_time, hc, c,
            DELTAT=deltat, CORRECTIONS=model_format)
        output[TYPE].data[:] += minor.data[:]
        #-- replace invalid values with fill value
        output[TYPE].data[output[TYPE].mask] = output[TYPE].fill_value

    #-- output to file
    if (FORMAT == 'csv'):
        pyTMD.spatial.to_ascii(output, attrib, output_file, delimiter=',',
            columns=['time','lat','lon',output_variable], verbose=VERBOSE)
    elif (FORMAT == 'netCDF4'):
        pyTMD.spatial.to_netCDF4(output, attrib, output_file, verbose=VERBOSE)
    elif (FORMAT == 'HDF5'):
        pyTMD.spatial.to_HDF5(output, attrib, output_file, verbose=VERBOSE)
    #-- change the permissions level to MODE
    os.chmod(output_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -T X, --tide=X\t\tTide model to use in correction')
    print(' --format=X\t\tInput and output data format')
    print('\tcsv\n\tnetCDF4\n\tHDF5')
    print(' --variables=X\t\tVariable names of data in input file')
    print(' --epoch=X\t\tReference epoch of input time')
    print('\tin form "time-units since yyyy-mm-dd hh:mm:ss"')
    print(' --projection=X\t\tSpatial projection as EPSG code or PROJ4 string')
    print(' -V, --verbose\t\tVerbose output of processing run')
    print(' -M X, --mode=X\t\tPermission mode of output file\n')

#-- Main program that calls compute_tidal_currents()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','tide=','variables=','epoch=',
        'projection=','format=','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:T:VM:', long_options)

    #-- command line options
    #-- set data directory containing the tidal data
    tide_dir = None
    #-- tide model to use
    TIDE_MODEL = 'CATS2008'
    #-- input and output data format
    FORMAT = 'csv'
    #-- variable names (for csv names of columns)
    VARIABLES = ['time','lat','lon','data']
    #-- time epoch (default Modified Julian Days)
    TIME_UNITS = 'days since 1858-11-17T00:00:00'
    #-- spatial projection (EPSG code or PROJ4 string)
    PROJECTION = '4326'
    #-- verbose output of processing run
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
            TIDE_MODEL = arg
        elif opt in ("--format",):
            FORMAT = arg
        elif opt in ("--variables",):
            VARIABLES = arg.split(',')
        elif opt in ("--epoch",):
            TIME_UNITS = arg
        elif opt in ("--projection",):
            PROJECTION = arg
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- enter input and output files as system argument
    if not arglist:
        raise Exception('No System Arguments Listed')

    #-- verify model before running program
    model_list=['CATS0201','CATS2008','TPXO9-atlas','TPXO9-atlas-v2','TPXO9.1',
        'TPXO8-atlas','TPXO7.2','AODTM-5','AOTIM-5','AOTIM-5-2018','FES2014']
    assert TIDE_MODEL in model_list, 'Unlisted tide model'

    #-- tilde-expand input and output files
    #-- set output file from input filename if not entered
    input_file=os.path.expanduser(arglist[0])
    try:
        output_file=os.path.expanduser(arglist[1])
    except IndexError:
        fileBasename,fileExtension = os.path.splitext(input_file)
        args = (fileBasename,TIDE_MODEL,fileExtension)
        output_file = '{0}_{1}_currents{2}'.format(*args)
    #-- set base directory from the input file if not currently set
    tide_dir = os.path.dirname(input_file) if (tide_dir is None) else tide_dir

    #-- run tidal current program for input file
    compute_tidal_currents(tide_dir, input_file, output_file, FORMAT=FORMAT,
        TIDE_MODEL=TIDE_MODEL, VARIABLES=VARIABLES, TIME_UNITS=TIME_UNITS,
        PROJECTION=PROJECTION, VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
