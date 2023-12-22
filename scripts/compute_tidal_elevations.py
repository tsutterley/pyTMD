#!/usr/bin/env python
u"""
compute_tidal_elevations.py
Written by Tyler Sutterley (12/2023)
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
    parquet file with variables for spatial and temporal coordinates
    geotiff file with bands in spatial coordinates

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -T X, --tide X: Tide model to use in correction
    --atlas-format X: ATLAS tide model format (OTIS, netcdf)
    --gzip, -G: Tide model files are gzip compressed
    --definition-file X: Model definition file for use as correction
    --format X: input and output data format
        csv (default)
        netCDF4
        HDF5
        geotiff
    --variables X: variable names of data in csv, HDF5 or netCDF4 file
        for csv files: the order of the columns within the file
        for HDF5, netCDF4 and parquet files: time, y, x and data variable names
    -H X, --header X: number of header lines for csv files
    --delimiter X: Delimiter for csv or ascii files
    -t X, --type X: input data type
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time for all data points)
    -e X, --epoch X: Reference epoch of input time (default Modified Julian Day)
        days since 1858-11-17T00:00:00
    -d X, --deltatime X: Input delta time for files without date information
        can be set to 0 to use exact calendar date from epoch
    -s X, --standard X: Input time standard for delta times or input time type
        UTC: Coordinate Universal Time
        GPS: GPS Time
        LORAN: Long Range Navigator Time
        TAI: International Atomic Time
        datetime: formatted datetime string in UTC
    -P X, --projection X: spatial projection as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
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
    -f X, --fill-value X: Invalid value for spatial fields
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
    pandas: Python Data Analysis Library
        https://pandas.pydata.org/
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    spatial: utilities for reading, writing and operating on spatial data
    utilities.py: download and management utilities for syncing files
    arguments.py: load the nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    crs.py: Coordinate Reference System (CRS) routines
    io/model.py: retrieves tide model parameters for named tide models
    io/OTIS.py: extract tidal harmonic constants from OTIS tide models
    io/ATLAS.py: extract tidal harmonic constants from netcdf models
    io/FES.py: extract tidal harmonic constants from FES tide models
    interpolate.py: interpolation routines for spatial data
    predict.py: predict tidal values using harmonic constants

UPDATE HISTORY:
    Updated 12/2023: use new crs class to get projection information 
    Updated 10/2023: can write datetime as time column for csv files
    Updated 08/2023: changed ESR netCDF4 format to TMD3 format
    Updated 05/2023: use timescale class for time conversion operations
    Updated 04/2023: check if datetime before converting to seconds
        using pathlib to define and expand paths
    Updated 02/2023: added functionality for time series type
    Updated 01/2023: added default field mapping for reading from netCDF4/HDF5
        added data type keyword for netCDF4 output
    Updated 12/2022: single implicit import of pyTMD tools
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 10/2022: added delimiter option and datetime parsing for ascii files
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        updated keyword arguments to read tide model programs
        added command line option to apply flexure for applicable models
    Updated 04/2022: use argparse descriptions within documentation
    Updated 03/2022: using static decorators to define available models
    Updated 02/2022: added Arctic 2km model (Arc2kmTM) to list of models
    Updated 01/2022: added option for changing the time standard
    Updated 12/2021: added TPXO9-atlas-v5 to list of available tide models
    Updated 11/2021: add function for attempting to extract projection
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: refactor to use model class for files and attributes
    Updated 07/2021: added tide model reference to output attributes
        can use prefix files to define command line arguments
    Updated 06/2021: added new Gr1km-v2 1km Greenland model from ESR
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: added TPXO9-atlas-v4 in binary OTIS format
        simplified netcdf inputs to be similar to binary OTIS read program
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
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
import logging
import pathlib
import argparse
import numpy as np
import pyTMD

# attempt imports
try:
    import pandas as pd
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    logging.critical("pandas not available")
try:
    import pyproj
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    logging.critical("pyproj not available")

# PURPOSE: try to get the projection information for the input file
def get_projection(attributes, PROJECTION):
    # coordinate reference system string from file
    try:
        crs = pyTMD.crs.from_input(attributes['projection'])
    except (ValueError,KeyError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # coordinate reference system from input argument
    try:
        crs = pyTMD.crs.from_input(PROJECTION)
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # no projection can be made
    raise pyproj.exceptions.CRSError

# compute tides at points and times using tidal model driver algorithms
def compute_tidal_elevations(tide_dir, input_file, output_file,
    TIDE_MODEL=None,
    ATLAS_FORMAT='netcdf',
    GZIP=True,
    DEFINITION_FILE=None,
    FORMAT='csv',
    VARIABLES=[],
    HEADER=0,
    DELIMITER=',',
    TYPE='drift',
    TIME_UNITS='days since 1858-11-17T00:00:00',
    TIME_STANDARD='UTC',
    TIME=None,
    PROJECTION='4326',
    METHOD='spline',
    EXTRAPOLATE=False,
    CUTOFF=None,
    APPLY_FLEXURE=False,
    FILL_VALUE=-9999.0,
    VERBOSE=False,
    MODE=0o775):

    # create logger for verbosity level
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # get parameters for tide model
    if DEFINITION_FILE is not None:
        model = pyTMD.io.model(tide_dir).from_file(DEFINITION_FILE)
    else:
        model = pyTMD.io.model(tide_dir, format=ATLAS_FORMAT,
            compressed=GZIP).elevation(TIDE_MODEL)

    # read input file to extract time, spatial coordinates and data
    if (FORMAT == 'csv'):
        parse_dates = (TIME_STANDARD.lower() == 'datetime')
        dinput = pyTMD.spatial.from_ascii(input_file, columns=VARIABLES,
            delimiter=DELIMITER, header=HEADER, parse_dates=parse_dates)
        attributes = dinput['attributes']
    elif (FORMAT == 'netCDF4'):
        field_mapping = pyTMD.spatial.default_field_mapping(VARIABLES)
        dinput = pyTMD.spatial.from_netCDF4(input_file,
            field_mapping=field_mapping)
        attributes = dinput['attributes']
    elif (FORMAT == 'HDF5'):
        field_mapping = pyTMD.spatial.default_field_mapping(VARIABLES)
        dinput = pyTMD.spatial.from_HDF5(input_file,
            field_mapping=field_mapping)
        attributes = dinput['attributes']
    elif (FORMAT == 'geotiff'):
        dinput = pyTMD.spatial.from_geotiff(input_file)
        attributes = dinput['attributes']
    elif (FORMAT == 'parquet'):
        logging.info(str(input_file))
        field_mapping = pyTMD.spatial.default_field_mapping(VARIABLES)
        remap = pyTMD.spatial.inverse_mapping(field_mapping)
        dinput = pd.read_parquet(input_file, columns=VARIABLES)
        dinput.rename(columns=remap, inplace=True)
        attributes = {}
    # update time variable if entered as argument
    if TIME is not None:
        dinput['time'] = np.copy(TIME)

    # converting x,y from projection to latitude/longitude
    crs1 = get_projection(attributes, PROJECTION)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    if (TYPE == 'grid'):
        ny, nx = (len(dinput['y']), len(dinput['x']))
        gridx, gridy = np.meshgrid(dinput['x'], dinput['y'])
        lon, lat = transformer.transform(gridx, gridy)
    elif (TYPE == 'drift'):
        lon, lat = transformer.transform(dinput['x'], dinput['y'])
    elif (TYPE == 'time series'):
        nstation = len(dinput['y'])
        lon, lat = transformer.transform(dinput['x'], dinput['y'])

    # extract time units from netCDF4 and HDF5 attributes or from TIME_UNITS
    try:
        time_string = attributes['time']['units']
        epoch1, to_secs = pyTMD.time.parse_date_string(time_string)
    except (TypeError, KeyError, ValueError):
        epoch1, to_secs = pyTMD.time.parse_date_string(TIME_UNITS)

    # convert delta times or datetimes objects to timescale
    if (TIME_STANDARD.lower() == 'datetime'):
        timescale = pyTMD.time.timescale().from_datetime(
            np.ravel(dinput['time']))
    else:
        # convert time to seconds
        delta_time = to_secs*np.ravel(dinput['time'])
        timescale = pyTMD.time.timescale().from_deltatime(delta_time,
            epoch=epoch1, standard=TIME_STANDARD)
    # number of time points
    nt = len(timescale)

    # read tidal constants and interpolate to grid points
    if model.format in ('OTIS','ATLAS','TMD3'):
        amp,ph,D,c = pyTMD.io.OTIS.extract_constants(np.ravel(lon), np.ravel(lat),
            model.grid_file, model.model_file, model.projection,
            type=model.type, method=METHOD, extrapolate=EXTRAPOLATE,
            cutoff=CUTOFF, grid=model.format, apply_flexure=APPLY_FLEXURE)
        deltat = np.zeros((nt))
    elif (model.format == 'netcdf'):
        amp,ph,D,c = pyTMD.io.ATLAS.extract_constants(np.ravel(lon), np.ravel(lat),
            model.grid_file, model.model_file, type=model.type,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            scale=model.scale, compressed=model.compressed)
        deltat = np.zeros((nt))
    elif (model.format == 'GOT'):
        amp,ph,c = pyTMD.io.GOT.extract_constants(np.ravel(lon), np.ravel(lat),
            model.model_file, method=METHOD, extrapolate=EXTRAPOLATE,
            cutoff=CUTOFF, scale=model.scale, compressed=model.compressed)
        # delta time (TT - UT1)
        deltat = timescale.tt_ut1
    elif (model.format == 'FES'):
        amp,ph = pyTMD.io.FES.extract_constants(np.ravel(lon), np.ravel(lat),
            model.model_file, type=model.type, version=model.version,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            scale=model.scale, compressed=model.compressed)
        # available model constituents
        c = model.constituents
        # delta time (TT - UT1)
        deltat = timescale.tt_ut1

    # calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    # calculate constituent oscillation
    hc = amp*np.exp(cph)

    # predict tidal elevations at time and infer minor corrections
    if (TYPE == 'grid'):
        tide = np.ma.zeros((ny,nx,nt), fill_value=FILL_VALUE)
        tide.mask = np.zeros((ny,nx,nt),dtype=bool)
        for i in range(nt):
            TIDE = pyTMD.predict.map(timescale.tide[i], hc, c,
                deltat=deltat[i], corrections=model.format)
            MINOR = pyTMD.predict.infer_minor(timescale.tide[i], hc, c,
                deltat=deltat[i], corrections=model.format)
            # add major and minor components and reform grid
            tide[:,:,i] = np.reshape((TIDE+MINOR), (ny,nx))
            tide.mask[:,:,i] = np.reshape((TIDE.mask | MINOR.mask), (ny,nx))
    elif (TYPE == 'drift'):
        tide = np.ma.zeros((nt), fill_value=FILL_VALUE)
        tide.mask = np.any(hc.mask,axis=1)
        tide.data[:] = pyTMD.predict.drift(timescale.tide, hc, c,
            deltat=deltat, corrections=model.format)
        minor = pyTMD.predict.infer_minor(timescale.tide, hc, c,
            deltat=deltat, corrections=model.format)
        tide.data[:] += minor.data[:]
    elif (TYPE == 'time series'):
        tide = np.ma.zeros((nstation,nt), fill_value=FILL_VALUE)
        tide.mask = np.zeros((nstation,nt),dtype=bool)
        for s in range(nstation):
            # calculate constituent oscillation for station
            TIDE = pyTMD.predict.time_series(timescale.tide, hc[s,None,:], c,
                deltat=deltat, corrections=model.format)
            MINOR = pyTMD.predict.infer_minor(timescale.tide, hc[s,None,:], c,
                deltat=deltat, corrections=model.format)
            tide.data[s,:] = TIDE.data[:] + MINOR.data[:]
            tide.mask[s,:] = (TIDE.mask | MINOR.mask)
    # replace invalid values with fill value
    tide.data[tide.mask] = tide.fill_value

    # output netCDF4 and HDF5 file attributes
    # will be added to YAML header in csv files
    attrib = {}
    # latitude
    attrib['lat'] = {}
    attrib['lat']['long_name'] = 'Latitude'
    attrib['lat']['units'] = 'Degrees_North'
    # longitude
    attrib['lon'] = {}
    attrib['lon']['long_name'] = 'Longitude'
    attrib['lon']['units'] = 'Degrees_East'
    # tides
    output_variable = model.variable
    attrib[output_variable] = {}
    attrib[output_variable]['description'] = model.description
    attrib[output_variable]['reference'] = model.reference
    attrib[output_variable]['model'] = model.name
    attrib[output_variable]['units'] = 'meters'
    attrib[output_variable]['long_name'] = model.long_name
    attrib[output_variable]['_FillValue'] = FILL_VALUE
    # time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['calendar'] = 'standard'

    # output data dictionary
    output = {'lon':lon, 'lat':lat, output_variable:tide}
    if (FORMAT == 'csv') and (TIME_STANDARD.lower() == 'datetime'):
        output['time'] = timescale.to_string()
    else:
        attrib['time']['units'] = 'days since 1992-01-01T00:00:00'
        output['time'] = timescale.tide

    # output to file
    if (FORMAT == 'csv'):
        pyTMD.spatial.to_ascii(output, attrib, output_file,
            delimiter=DELIMITER, header=False,
            columns=['time','lat','lon',output_variable])
    elif (FORMAT == 'netCDF4'):
        pyTMD.spatial.to_netCDF4(output, attrib, output_file, data_type=TYPE)
    elif (FORMAT == 'HDF5'):
        pyTMD.spatial.to_HDF5(output, attrib, output_file)
    elif (FORMAT == 'geotiff'):
        # copy global geotiff attributes for projection and grid parameters
        for att_name in ['projection','wkt','spacing','extent']:
            attrib[att_name] = attributes[att_name]
        pyTMD.spatial.to_geotiff(output, attrib, output_file,
            varname=output_variable)
    elif (FORMAT == 'parquet'):
        # write to parquet file
        logging.info(str(output_file))
        pd.DataFrame(output).to_parquet(output_file)
    # change the permissions level to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates tidal elevations for an input file
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    group = parser.add_mutually_exclusive_group(required=True)
    # command line options
    # input and output file
    parser.add_argument('infile',
        type=pathlib.Path, nargs='?',
        help='Input file to run')
    parser.add_argument('outfile',
        type=pathlib.Path, nargs='?',
        help='Computed output file')
    # set data directory containing the tidal data
    parser.add_argument('--directory','-D',
        type=pathlib.Path,
        help='Working data directory')
    # tide model to use
    choices = sorted(pyTMD.io.model.ocean_elevation() +
                     pyTMD.io.model.load_elevation())
    group.add_argument('--tide','-T',
        type=str, choices=choices,
        help='Tide model to use in correction')
    parser.add_argument('--atlas-format',
        type=str, choices=('OTIS','netcdf'), default='netcdf',
        help='ATLAS tide model format')
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='Tide model files are gzip compressed')
    # tide model definition file to set an undefined model
    group.add_argument('--definition-file',
        type=pathlib.Path,
        help='Tide model definition file')
    # input and output data format
    parser.add_argument('--format','-F',
        type=str, default='csv',
        choices=('csv','netCDF4','HDF5','geotiff','parquet'),
        help='Input and output data format')
    # variable names (for csv names of columns)
    parser.add_argument('--variables','-v',
        type=str, nargs='+', default=['time','lat','lon','data'],
        help='Variable names of data in input file')
    # number of header lines for csv files
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of header lines for csv files')
    # delimiter for csv or ascii files
    parser.add_argument('--delimiter',
        type=str, default=',',
        help='Delimiter for csv or ascii files')
    # input data type
    # drift: drift buoys or satellite/airborne altimetry (time per data point)
    # grid: spatial grids or images (single time for all data points)
    # time series: station locations with multiple time values
    parser.add_argument('--type','-t',
        type=str, default='drift',
        choices=('drift','grid','time series'),
        help='Input data type')
    # time epoch (default Modified Julian Days)
    # in form "time-units since yyyy-mm-dd hh:mm:ss"
    parser.add_argument('--epoch','-e',
        type=str, default='days since 1858-11-17T00:00:00',
        help='Reference epoch of input time')
    # input delta time for files without date information
    parser.add_argument('--deltatime','-d',
        type=float, nargs='+',
        help='Input delta time for files without date variables')
    # input time standard definition
    parser.add_argument('--standard','-s',
        type=str, choices=('UTC','GPS','TAI','LORAN','datetime'), default='UTC',
        help='Input time standard for delta times')
    # spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection','-P',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
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
    # fill value for output spatial fields
    parser.add_argument('--fill-value','-f',
        type=float, default=-9999.0,
        help='Invalid value for spatial fields')
    # verbose output of processing run
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    # permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of output file')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # set output file from input filename if not entered
    if not args.outfile:
        flexure_flag = '_FLEXURE' if args.apply_flexure else ''
        vars = (args.infile.stem,args.tide,flexure_flag,args.infile.suffix)
        args.outfile = args.infile.with_name('{0}_{1}{2}{3}'.format(*vars))

    # run tidal elevation program for input file
    compute_tidal_elevations(args.directory, args.infile, args.outfile,
        TIDE_MODEL=args.tide,
        ATLAS_FORMAT=args.atlas_format,
        GZIP=args.gzip,
        DEFINITION_FILE=args.definition_file,
        FORMAT=args.format,
        VARIABLES=args.variables,
        HEADER=args.header,
        DELIMITER=args.delimiter,
        TYPE=args.type,
        TIME_UNITS=args.epoch,
        TIME=args.deltatime,
        TIME_STANDARD=args.standard,
        PROJECTION=args.projection,
        METHOD=args.interpolate,
        EXTRAPOLATE=args.extrapolate,
        CUTOFF=args.cutoff,
        APPLY_FLEXURE=args.apply_flexure,
        FILL_VALUE=args.fill_value,
        VERBOSE=args.verbose,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
