#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (09/2020)

Utilities for reading and writing spatial data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
        (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    h5py: Pythonic interface to the HDF5 binary data format.
        (https://www.h5py.org/)Data class
    PyYAML: YAML parser and emitter for Python
        https://github.com/yaml/pyyaml

UPDATE HISTORY:
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
        ascii file is compressed using gzip
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
        dinput['attributes'] = None
    #-- extract spatial data array
    #-- for each line in the file
    for i,line in enumerate(file_contents[header:]):
        #-- extract columns of interest and assign to dict
        #-- convert fortran exponentials if applicable
        column = {c:r.replace('D','E') for c,r in zip(columns,rx.findall(line))}
        #-- copy variables from column dict to output dictionary
        for c in columns:
            dinput[c][i] = np.float(column[c])
    #-- return the spatial variables
    return dinput

def from_netCDF4(filename, compression=None, verbose=False,
    timename='time', xname='lon', yname='lat', varname='data'):
    """
    Read data from a netCDF4 file
    Inputs: full path of input netCDF4 file
    Options:
        netCDF4 file is compressed using gzip
        verbose output of file information
        netCDF4 variable names of time, longitude, latitude, and data
    """
    #-- read data from netCDF4 file
    #-- Open the NetCDF4 file for reading
    if (compression == 'gzip'):
        #-- read as in-memory (diskless) netCDF4 dataset
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex,memory=f.read())
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
        HDF5 file is compressed using gzip
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
        #-- Getting the data from each NetCDF variable
        dinput[key] = np.copy(fileID[h5][:])
        #-- get attributes for the included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            #-- try getting the attribute
            try:
                dinput['attributes'][key][attr] = fileID[h5].attrs[attr]
            except (KeyError,AttributeError):
                pass
    #-- Closing the NetCDF file
    fileID.close()
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
