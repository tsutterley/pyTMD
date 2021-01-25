#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (01/2021)

Utilities for reading, writing and operating on spatial data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format
        https://www.h5py.org/
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL
    PyYAML: YAML parser and emitter for Python
        https://github.com/yaml/pyyaml

UPDATE HISTORY:
    Updated 01/2020: add streaming from bytes for ascii, netCDF4, HDF5, geotiff
        set default time for geotiff files to 0
    Updated 12/2020: added module for converting ellipsoids
    Updated 11/2020: output data as masked arrays if containing fill values
        add functions to read from and write to geotiff image formats
    Written 09/2020
"""
import os
import re
import io
import gzip
import uuid
import h5py
import yaml
import netCDF4
import datetime
import numpy as np
import osgeo.gdal, osgeo.osr

def case_insensitive_filename(filename):
    """
    Searches a directory for a filename without case dependence
    """
    #-- check if file presently exists with input case
    if not os.access(os.path.expanduser(filename),os.F_OK):
        #-- search for filename without case dependence
        basename = os.path.basename(filename)
        directory = os.path.dirname(os.path.expanduser(filename))
        f = [f for f in os.listdir(directory) if re.match(basename,f,re.I)]
        if not f:
            raise IOError('{0} not found in file system'.format(filename))
        filename = os.path.join(directory,f.pop())
    return os.path.expanduser(filename)

def from_ascii(filename, compression=None, verbose=False,
    columns=['time','y','x','data'], header=0):
    """
    Read data from an ascii file
    Inputs: full path of input ascii file
    Options:
        ascii file is compressed or streamed from memory
        verbose output of file information
        column names of ascii file
        header lines to skip from start of file
    """
    #-- set filename
    print(filename) if verbose else None
    #-- open the ascii file and extract contents
    if (compression == 'gzip'):
        #-- read input ascii data from gzip compressed file and split lines
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            file_contents = f.read().decode('ISO-8859-1').splitlines()
    elif (compression == 'bytes'):
        #-- read input file object and split lines
        file_contents = filename.read().splitlines()
    else:
        #-- read input ascii file (.txt, .asc) and split lines
        with open(case_insensitive_filename(filename),'r') as f:
            file_contents = f.read().splitlines()
    #-- number of lines in the file
    file_lines = len(file_contents)
    #-- compile regular expression operator for extracting numerical values
    #-- from input ascii files of spatial data
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)
    #-- check if header has a known format
    if (str(header).upper() == 'YAML'):
        #-- counts the number of lines in the header
        YAML = False
        count = 0
        #-- Reading over header text
        while (YAML is False) & (count < file_lines):
            #-- file line at count
            line = file_contents[count]
            #-- if End of YAML Header is found: set YAML flag
            YAML = bool(re.search("\# End of YAML header",line))
            #-- add 1 to counter
            count += 1
        #-- parse the YAML header (specifying yaml loader)
        YAML_HEADER = yaml.load('\n'.join(file_contents[:count]),
           Loader=yaml.BaseLoader)
        #-- output spatial data and attributes
        dinput = {}
        #-- copy global attributes
        dinput['attributes'] = YAML_HEADER['header']['global_attributes']
        #-- allocate for each variable and copy variable attributes
        for c in columns:
            dinput[c] = np.zeros((file_lines-count))
            dinput['attributes'][c] = YAML_HEADER['header']['variables'][c]
        #-- update number of file lines to skip for reading data
        header = np.int(count)
    else:
        #-- output spatial data and attributes
        dinput = {c:np.zeros((file_lines-header)) for c in columns}
        dinput['attributes'] = {c:dict() for c in columns}
    #-- extract spatial data array
    #-- for each line in the file
    for i,line in enumerate(file_contents[header:]):
        #-- extract columns of interest and assign to dict
        #-- convert fortran exponentials if applicable
        column = {c:r.replace('D','E') for c,r in zip(columns,rx.findall(line))}
        #-- copy variables from column dict to output dictionary
        for c in columns:
            dinput[c][i] = np.float(column[c])
    #-- convert to masked array if fill values
    if '_FillValue' in dinput['attributes']['data'].keys():
        dinput['data'] = np.ma.asarray(dinput['data'])
        dinput['data'].fill_value = dinput['attributes']['data']['_FillValue']
        dinput['data'].mask = (dinput['data'].data == dinput['data'].fill_value)
    #-- return the spatial variables
    return dinput

def from_netCDF4(filename, compression=None, verbose=False,
    timename='time', xname='lon', yname='lat', varname='data'):
    """
    Read data from a netCDF4 file
    Inputs: full path of input netCDF4 file
    Options:
        netCDF4 file is compressed or streamed from memory
        verbose output of file information
        netCDF4 variable names of time, longitude, latitude, and data
    """
    #-- read data from netCDF4 file
    #-- Open the NetCDF4 file for reading
    if (compression == 'gzip'):
        #-- read as in-memory (diskless) netCDF4 dataset
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex,memory=f.read())
    elif (compression == 'bytes'):
        #-- read as in-memory (diskless) netCDF4 dataset
        fileID = netCDF4.Dataset(uuid.uuid4().hex,memory=filename.read())
    else:
        #-- read netCDF4 dataset
        fileID = netCDF4.Dataset(case_insensitive_filename(filename), 'r')
    #-- Output NetCDF file information
    if verbose:
        print(fileID.filepath())
        print(list(fileID.variables.keys()))
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {}
    #-- get attributes for the file
    for attr in ['title','description','projection']:
        #-- try getting the attribute
        try:
            ncattr, = [s for s in dir(fileID) if re.match(attr,s,re.I)]
            dinput['attributes'][attr] = getattr(fileID,ncattr)
        except (ValueError,AttributeError):
            pass
    #-- list of attributes to attempt to retrieve from included variables
    attributes_list = ['description','units','long_name','calendar',
        'standard_name','_FillValue']
    #-- mapping between netCDF4 variable names and output names
    variable_mapping = dict(x=xname,y=yname,data=varname,time=timename)
    #-- for each variable
    for key,nc in variable_mapping.items():
        #-- Getting the data from each NetCDF variable
        dinput[key] = fileID.variables[nc][:]
        #-- get attributes for the included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            #-- try getting the attribute
            try:
                ncattr, = [s for s in dir(fileID) if re.match(attr,s,re.I)]
                dinput['attributes'][key][attr] = \
                    getattr(fileID.variables[nc],ncattr)
            except (ValueError,AttributeError):
                pass
    #-- convert to masked array if fill values
    if '_FillValue' in dinput['attributes']['data'].keys():
        dinput['data'] = np.ma.asarray(dinput['data'])
        dinput['data'].fill_value = dinput['attributes']['data']['_FillValue']
        dinput['data'].mask = (dinput['data'].data == dinput['data'].fill_value)
    #-- Closing the NetCDF file
    fileID.close()
    #-- return the spatial variables
    return dinput

def from_HDF5(filename, compression=None, verbose=False,
    timename='time', xname='lon', yname='lat', varname='data'):
    """
    Read data from a HDF5 file
    Inputs: full path of input HDF5 file
    Options:
        HDF5 file is compressed or streamed from memory
        verbose output of file information
        HDF5 variable names of time, longitude, latitude, and data
    """
    #-- read data from HDF5 file
    #-- Open the HDF5 file for reading
    if (compression == 'gzip'):
        #-- read gzip compressed file and extract into in-memory file object
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            fid = io.BytesIO(f.read())
        #-- set filename of BytesIO object
        fid.filename = os.path.basename(filename)
        #-- rewind to start of file
        fid.seek(0)
        #-- read as in-memory (diskless) HDF5 dataset from BytesIO object
        fileID = h5py.File(fid, 'r')
    elif (compression == 'bytes'):
        #-- read as in-memory (diskless) HDF5 dataset
        fileID = h5py.File(filename, 'r')
    else:
        #-- read HDF5 dataset
        fileID = h5py.File(case_insensitive_filename(filename), 'r')
    #-- Output HDF5 file information
    if verbose:
        print(fileID.filename)
        print(list(fileID.keys()))
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {}
    #-- get attributes for the file
    for attr in ['title','description','projection']:
        #-- try getting the attribute
        try:
            dinput['attributes'][attr] = fileID.attrs[attr]
        except (KeyError,AttributeError):
            pass
    #-- list of attributes to attempt to retrieve from included variables
    attributes_list = ['description','units','long_name','calendar',
        'standard_name','_FillValue']
    #-- mapping between HDF5 variable names and output names
    variable_mapping = dict(x=xname,y=yname,data=varname,time=timename)
    #-- for each variable
    for key,h5 in variable_mapping.items():
        #-- Getting the data from each HDF5 variable
        dinput[key] = np.copy(fileID[h5][:])
        #-- get attributes for the included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            #-- try getting the attribute
            try:
                dinput['attributes'][key][attr] = fileID[h5].attrs[attr]
            except (KeyError,AttributeError):
                pass
    #-- convert to masked array if fill values
    if '_FillValue' in dinput['attributes']['data'].keys():
        dinput['data'] = np.ma.asarray(dinput['data'])
        dinput['data'].fill_value = dinput['attributes']['data']['_FillValue']
        dinput['data'].mask = (dinput['data'].data == dinput['data'].fill_value)
    #-- Closing the HDF5 file
    fileID.close()
    #-- return the spatial variables
    return dinput

def from_geotiff(filename, compression=None, verbose=False):
    """
    Read data from a geotiff file
    Inputs: full path of input geotiff file
    Options:
        geotiff file is compressed or streamed from memory
        verbose output of file information
    """
    #-- Open the geotiff file for reading
    if (compression == 'gzip'):
        #-- read gzip compressed file and extract into memory-mapped object
        mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            osgeo.gdal.FileFromMemBuffer(mmap_name, f.read())
        #-- read as GDAL memory-mapped (diskless) geotiff dataset
        ds = osgeo.gdal.Open(mmap_name)
    elif (compression == 'bytes'):
        #-- read as GDAL memory-mapped (diskless) geotiff dataset
        mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
        osgeo.gdal.FileFromMemBuffer(mmap_name, filename.read())
        ds = osgeo.gdal.Open(mmap_name)
    else:
        #-- read geotiff dataset
        ds = osgeo.gdal.Open(case_insensitive_filename(filename))
    #-- print geotiff file if verbose
    print(filename) if verbose else None
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {c:dict() for c in ['x','y','data']}
    #-- get the spatial projection reference information
    srs = ds.GetSpatialRef()
    dinput['attributes']['projection'] = srs.ExportToProj4()
    dinput['attributes']['wkt'] = srs.ExportToWkt()
    #-- get dimensions
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    bsize = ds.RasterCount
    #-- get geotiff info
    info_geotiff = ds.GetGeoTransform()
    dinput['attributes']['spacing'] = (info_geotiff[1],info_geotiff[5])
    #-- calculate image extents
    xmin = info_geotiff[0]
    ymax = info_geotiff[3]
    xmax = xmin + (xsize-1)*info_geotiff[1]
    ymin = ymax + (ysize-1)*info_geotiff[5]
    dinput['attributes']['extent'] = (xmin,xmax,ymin,ymax)
    #-- x and y pixel center coordinates (converted from upper left)
    dinput['x'] = xmin + info_geotiff[1]/2.0 + np.arange(xsize)*info_geotiff[1]
    dinput['y'] = ymax + info_geotiff[5]/2.0 + np.arange(ysize)*info_geotiff[5]
    #-- read full image with GDAL
    dinput['data'] = ds.ReadAsArray()
    #-- set default time to zero for each band
    dinput.setdefault('time', np.zeros((bsize)))
    #-- check if image has fill values
    if ds.GetRasterBand(1).GetNoDataValue():
        #-- convert to masked array if fill values
        dinput['data'] = np.ma.asarray(dinput['data'])
        #-- mask invalid values
        dinput['data'].fill_value = ds.GetRasterBand(1).GetNoDataValue()
        #-- create mask array for bad values
        dinput['data'].mask = (dinput['data'].data == dinput['data'].fill_value)
        #-- set attribute for fill value
        dinput['attributes']['data']['_FillValue'] = dinput['data'].fill_value
    #-- close the dataset
    ds = None
    #-- return the spatial variables
    return dinput

def to_ascii(output, attributes, filename, delimiter=',',
    columns=['time','lat','lon','tide'], header=False, verbose=False):
    """
    Write data to an ascii file
    Inputs:
        python dictionary of output data
        python dictionary of output attributes
        full path of output ascii file
    Options:
        delimiter for output spatial file
        order of columns for output spatial file
        create a YAML header with data attributes
        verbose output
    """
    filename = os.path.expanduser(filename)
    print(filename) if verbose else None
    #-- open the output file
    fid = open(filename, 'w')
    #-- create a column stack arranging data in column order
    data_stack = np.c_[[output[col] for col in columns]]
    ncol,nrow = np.shape(data_stack)
    #-- print YAML header to top of file
    if header:
        fid.write('{0}:\n'.format('header'))
        #-- data dimensions
        fid.write('\n  {0}:\n'.format('dimensions'))
        fid.write('    {0:22}: {1:d}\n'.format('time',nrow))
        #-- non-standard attributes
        fid.write('  {0}:\n'.format('non-standard_attributes'))
        #-- data format
        fid.write('    {0:22}: ({1:d}f0.8)\n'.format('formatting_string',ncol))
        fid.write('\n')
        #-- global attributes
        fid.write('\n  {0}:\n'.format('global_attributes'))
        today = datetime.datetime.now().isoformat()
        fid.write('    {0:22}: {1}\n'.format('date_created', today))
        # print variable descriptions to YAML header
        fid.write('\n  {0}:\n'.format('variables'))
        #-- print YAML header with variable attributes
        for i,v in enumerate(columns):
            fid.write('    {0:22}:\n'.format(v))
            for atn,atv in attributes[v].items():
                fid.write('      {0:20}: {1}\n'.format(atn,atv))
            #-- add precision and column attributes for ascii yaml header
            fid.write('      {0:20}: double_precision\n'.format('precision'))
            fid.write('      {0:20}: column {1:d}\n'.format('comments',i+1))
        #-- end of header
        fid.write('\n\n# End of YAML header\n')
    #-- write to file for each data point
    for line in range(nrow):
        line_contents = ['{0:0.8f}'.format(d) for d in data_stack[:,line]]
        print(delimiter.join(line_contents), file=fid)
    #-- close the output file
    fid.close()

def to_netCDF4(output, attributes, filename, verbose=False):
    """
    Write data to a netCDF4 file
    Inputs:
        python dictionary of output data
        python dictionary of output attributes
        full path of output netCDF4 file
    Options: verbose output
    """
    #-- opening NetCDF file for writing
    fileID = netCDF4.Dataset(os.path.expanduser(filename),'w',format="NETCDF4")
    #-- Defining the NetCDF dimensions
    fileID.createDimension('time', len(np.atleast_1d(output['time'])))
    #-- defining the NetCDF variables
    nc = {}
    for key,val in output.items():
        if '_FillValue' in attributes[key].keys():
            nc[key] = fileID.createVariable(key, val.dtype, ('time',),
                fill_value=attributes[key]['_FillValue'], zlib=True)
        else:
            nc[key] = fileID.createVariable(key, val.dtype, ('time',))
        #-- filling NetCDF variables
        nc[key][:] = val
        #-- Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            setattr(nc[key],att_name,att_val)
    #-- add attribute for date created
    fileID.date_created = datetime.datetime.now().isoformat()
    #-- Output NetCDF structure information
    if verbose:
        print(filename)
        print(list(fileID.variables.keys()))
    #-- Closing the NetCDF file
    fileID.close()

def to_HDF5(output, attributes, filename, verbose=False):
    """
    Write data to a HDF5 file
    Inputs:
        python dictionary of output data
        python dictionary of output attributes
        full path of output HDF5 file
    Options: verbose output
    """
    #-- opening HDF5 file for writing
    fileID = h5py.File(filename, 'w')
    #-- Defining the HDF5 dataset variables
    h5 = {}
    for key,val in output.items():
        if '_FillValue' in attributes[key].keys():
            h5[key] = fileID.create_dataset(key, val.shape, data=val,
                dtype=val.dtype, fillvalue=attributes[key]['_FillValue'],
                compression='gzip')
        else:
            h5[key] = fileID.create_dataset(key, val.shape, data=val,
                dtype=val.dtype, compression='gzip')
        #-- Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            h5[key].attrs[att_name] = att_val
    #-- add attribute for date created
    fileID.attrs['date_created'] = datetime.datetime.now().isoformat()
    #-- Output HDF5 structure information
    if verbose:
        print(filename)
        print(list(fileID.keys()))
    #-- Closing the HDF5 file
    fileID.close()

def to_geotiff(output, attributes, filename, verbose=False,
    varname='data', dtype=osgeo.gdal.GDT_Float64):
    """
    Write data to a geotiff file
    Inputs:
        python dictionary of output data
        python dictionary of output attributes
        full path of output geotiff file
    Options:
        verbose output
        output variable name
        GDAL data type
    """
    #-- verify grid dimensions to be iterable
    output = expand_dims(output, varname=varname)
    #-- grid shape
    ny,nx,nband = np.shape(output[varname])
    #-- output as geotiff
    driver = osgeo.gdal.GetDriverByName("GTiff")
    #-- set up the dataset with compression options
    ds = driver.Create(filename,nx,ny,nband,dtype,['COMPRESS=LZW'])
    #-- top left x, w-e pixel resolution, rotation
    #-- top left y, rotation, n-s pixel resolution
    xmin,xmax,ymin,ymax = attributes['extent']
    dx,dy = attributes['spacing']
    ds.SetGeoTransform([xmin,dx,0,ymax,0,dy])
    #-- set the spatial projection reference information
    srs = osgeo.osr.SpatialReference()
    srs.ImportFromWkt(attributes['wkt'])
    #-- export
    ds.SetProjection( srs.ExportToWkt() )
    #-- for each band
    for band in range(nband):
        #-- set fill value for band
        if '_FillValue' in attributes[varname].keys():
            fill_value = attributes[varname]['_FillValue']
            ds.GetRasterBand(band+1).SetNoDataValue(fill_value)
        #-- write band to geotiff array
        ds.GetRasterBand(band+1).WriteArray(output[varname][:,:,band])
    #-- print filename if verbose
    print(filename) if verbose else None
    #-- close dataset
    ds.FlushCache()

def expand_dims(obj, varname='data'):
    """
    Add a singleton dimension to a spatial dictionary if non-existent
    Options:
        variable name to modify
    """
    #-- change time dimensions to be iterableinformation
    try:
        obj['time'] = np.atleast_1d(obj['time'])
    except:
        pass
    #-- output spatial with a third dimension
    if isinstance(varname,list):
        for v in varname:
            obj[v] = np.atleast_3d(obj[v])
    elif isinstance(varname,str):
        obj[varname] = np.atleast_3d(obj[varname])
    #-- return reformed spatial dictionary
    return obj

def convert_ellipsoid(phi1, h1, a1, f1, a2, f2, eps=1e-12, itmax=10):
    """
    Convert latitudes and heights to a different ellipsoid using Newton-Raphson

    Inputs:
        phi1: latitude of input ellipsoid in degrees
        h1: height above input ellipsoid in meters
        a1: semi-major axis of input ellipsoid
        f1: flattening of input ellipsoid
        a2: semi-major axis of output ellipsoid
        f2: flattening of output ellipsoid

    Options:
        eps: tolerance to prevent division by small numbers
            and to determine convergence
        itmax: maximum number of iterations to use in Newton-Raphson

    Returns:
        phi2: latitude of output ellipsoid in degrees
        h2: height above output ellipsoid in meters

    References:
        Astronomical Algorithms, Jean Meeus, 1991, Willmann-Bell, Inc.
            pp. 77-82
    """
    if (len(phi1) != len(h1)):
        raise ValueError('phi and h have incompatable dimensions')
    #-- semiminor axis of input and output ellipsoid
    b1 = (1.0 - f1)*a1
    b2 = (1.0 - f2)*a2
    #-- initialize output arrays
    npts = len(phi1)
    phi2 = np.zeros((npts))
    h2 = np.zeros((npts))
    #-- for each point
    for N in range(npts):
        #-- force phi1 into range -90 <= phi1 <= 90
        if (np.abs(phi1[N]) > 90.0):
            phi1[N] = np.sign(phi1[N])*90.0
        #-- handle special case near the equator
        #-- phi2 = phi1 (latitudes congruent)
        #-- h2 = h1 + a1 - a2
        if (np.abs(phi1[N]) < eps):
            phi2[N] = np.copy(phi1[N])
            h2[N] = h1[N] + a1 - a2
        #-- handle special case near the poles
        #-- phi2 = phi1 (latitudes congruent)
        #-- h2 = h1 + b1 - b2
        elif ((90.0 - np.abs(phi1[N])) < eps):
            phi2[N] = np.copy(phi1[N])
            h2[N] = h1[N] + b1 - b2
        #-- handle case if latitude is within 45 degrees of equator
        elif (np.abs(phi1[N]) <= 45):
            #-- convert phi1 to radians
            phi1r = phi1[N] * np.pi/180.0
            sinphi1 = np.sin(phi1r)
            cosphi1 = np.cos(phi1r)
            #-- prevent division by very small numbers
            cosphi1 = np.copy(eps) if (cosphi1 < eps) else cosphi1
            #-- calculate tangent
            tanphi1 = sinphi1 / cosphi1
            u1 = np.arctan(b1 / a1 * tanphi1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinphi1
            hpr1cos = a1 * np.cos(u1) + h1[N] * cosphi1
            #-- set initial value for u2
            u2 = np.copy(u1)
            #-- setup constants
            k0 = b2 * b2 - a2 * a2
            k1 = a2 * hpr1cos
            k2 = b2 * hpr1sin
            #-- perform newton-raphson iteration to solve for u2
            #-- cos(u2) will not be close to zero since abs(phi1) <= 45
            for i in range(0, itmax+1):
                cosu2 = np.cos(u2)
                fu2 = k0 * np.sin(u2) + k1 * np.tan(u2) - k2
                fu2p = k0 * cosu2 + k1 / (cosu2 * cosu2)
                if (np.abs(fu2p) < eps):
                    i = np.copy(itmax)
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if (np.abs(delta) < eps):
                        i = np.copy(itmax)
            #-- convert latitude to degrees and verify values between +/- 90
            phi2r = np.arctan(a2 / b2 * np.tan(u2))
            phi2[N] = phi2r*180.0/np.pi
            if (np.abs(phi2[N]) > 90.0):
                phi2[N] = np.sign(phi2[N])*90.0
            #-- calculate height
            h2[N] = (hpr1cos - a2 * np.cos(u2)) / np.cos(phi2r)
        #-- handle final case where latitudes are between 45 degrees and pole
        else:
            #-- convert phi1 to radians
            phi1r = phi1[N] * np.pi/180.0
            sinphi1 = np.sin(phi1r)
            cosphi1 = np.cos(phi1r)
            #-- prevent division by very small numbers
            cosphi1 = np.copy(eps) if (cosphi1 < eps) else cosphi1
            #-- calculate tangent
            tanphi1 = sinphi1 / cosphi1
            u1 = np.arctan(b1 / a1 * tanphi1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinphi1
            hpr1cos = a1 * np.cos(u1) + h1[N] * cosphi1
            #-- set initial value for u2
            u2 = np.copy(u1)
            #-- setup constants
            k0 = a2 * a2 - b2 * b2
            k1 = b2 * hpr1sin
            k2 = a2 * hpr1cos
            #-- perform newton-raphson iteration to solve for u2
            #-- sin(u2) will not be close to zero since abs(phi1) > 45
            for i in range(0, itmax+1):
                sinu2 = np.sin(u2)
                fu2 = k0 * np.cos(u2) + k1 / np.tan(u2) - k2
                fu2p =  -1 * (k0 * sinu2 + k1 / (sinu2 * sinu2))
                if (np.abs(fu2p) < eps):
                    i = np.copy(itmax)
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if (np.abs(delta) < eps):
                        i = np.copy(itmax)
            #-- convert latitude to degrees and verify values between +/- 90
            phi2r = np.arctan(a2 / b2 * np.tan(u2))
            phi2[N] = phi2r*180.0/np.pi
            if (np.abs(phi2[N]) > 90.0):
                phi2[N] = np.sign(phi2[N])*90.0
            #-- calculate height
            h2[N] = (hpr1sin - b2 * np.sin(u2)) / np.sin(phi2r)

    #-- return the latitude and height
    return (phi2, h2)
