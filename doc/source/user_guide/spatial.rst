==========
spatial.py
==========

Utilities for reading and writing spatial data

 - Can read ascii, netCDF4, HDF5 files
 - Can output to ascii, netCDF4 or HDF5 files

Calling Sequence
================

Reading a netCDF4 file

.. code-block:: python

    import pyTMD.spatial
    dinput = pyTMD.spatial.from_netCDF4(path_to_netCDF4_file)

Reading a HDF5 file

.. code-block:: python

    import pyTMD.spatial
    dinput = pyTMD.spatial.from_HDF5(path_to_netCDF4_file)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/spatial.py

General Methods
===============


    .. method:: pyTMD.spatial.case_insensitive_filename(filename)

        Searches a directory for a filename without case dependence


    .. method:: pyTMD.spatial.from_ascii(filename, compression=None, verbose=False, columns=['time','y','x','data'], header=0)

        Read data from an ascii file

        Inputs: full path of input ascii file

        Options:
            `compression` ascii file is compressed using gzip

            `verbose` print ascii filename

            `columns` variable names for each column

            `header` lines to skip from start of file


    .. method:: pyTMD.spatial.from_netCDF4(filename, compression=None, verbose=False, timename='time', xname='lon', yname='lat', varname='data')

        Read data from a netCDF4 file

        Inputs: full path of input netCDF4 file

        Options:
            `compression` netCDF4 file is compressed using gzip

            `verbose` print netCDF4 file information

            `time` input time variable name in netCDF4 file

            `xname` input x variable name in netCDF4 file

            `yname` input y variable units in netCDF4 file

            `varname` input data variable units in netCDF4 file


    .. method:: pyTMD.spatial.from_HDF5(filename, compression=None, verbose=False, timename='time', xname='lon', yname='lat', varname='data')

        Read data from a HDF5 file

        Inputs: full path of input HDF5 file

        Options:
            `compression` HDF5 file is compressed using gzip

            `verbose` print HDF5 file information

            `time` input time variable name in HDF5 file

            `xname` input x variable name in HDF5 file

            `yname` input y variable units in HDF5 file

            `varname` input data variable units in HDF5 file


    .. method:: pyTMD.spatial.to_ascii(output, attributes, filename, delimiter=',', columns=['time','lat','lon','tide'], header=False, verbose=False)

        Write data to an ascii file

        Inputs:

            `output` python dictionary of output data

            `attributes` python dictionary of output attributes

            `filename` full path of output ascii file

        Options:

            `delimiter` for output spatial file

            `columns` order of columns for output spatial file

            `header` create a YAML header with data attributes

            `verbose` print ascii file name


    .. method:: pyTMD.spatial.to_netCDF4(output, attributes, filename, verbose=False)

        Write data to a netCDF4 file

        Inputs:

            `output` python dictionary of output data

            `attributes` python dictionary of output attributes

            `filename` full path of output netCDF4 file

        Options:

            `verbose` print netCDF4 file information


    .. method:: pyTMD.spatial.to_HDF5(output, attributes, filename, verbose=False)

        Write data to a HDF5 file

        Inputs:

            `output` python dictionary of output data

            `attributes` python dictionary of output attributes

            `filename` full path of output HDF5 file

        Options:

            `verbose` print HDF5 file information
