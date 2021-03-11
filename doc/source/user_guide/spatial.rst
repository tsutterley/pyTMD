==========
spatial.py
==========

Utilities for reading, writing and operating on spatial data

 - Can read ascii, netCDF4, HDF5 or geotiff files
 - Can output to ascii, netCDF4, HDF5 or geotiff files

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
        `compression` ascii file is compressed or streamed from memory

        `verbose` print ascii filename

        `columns` variable names for each column

        `header` lines to skip from start of file


.. method:: pyTMD.spatial.from_netCDF4(filename, compression=None, verbose=False, timename='time', xname='lon', yname='lat', varname='data')

    Read data from a netCDF4 file

    Inputs: full path of input netCDF4 file

    Options:
        `compression` netCDF4 file is compressed or streamed from memory

        `verbose` print netCDF4 file information

        `time` input time variable name in netCDF4 file

        `xname` input x variable name in netCDF4 file

        `yname` input y variable units in netCDF4 file

        `varname` input data variable units in netCDF4 file


.. method:: pyTMD.spatial.from_HDF5(filename, compression=None, verbose=False, timename='time', xname='lon', yname='lat', varname='data')

    Read data from a HDF5 file

    Inputs: full path of input HDF5 file

    Options:
        `compression` HDF5 file is compressed or streamed from memory

        `verbose` print HDF5 file information

        `time` input time variable name in HDF5 file

        `xname` input x variable name in HDF5 file

        `yname` input y variable units in HDF5 file

        `varname` input data variable units in HDF5 file


.. method:: pyTMD.spatial.from_geotiff(filename, compression=None, verbose=False)

    Read data from a geotiff file

    Inputs: full path of input geotiff file

    Options:
        `compression` geotiff file is compressed using gzip

        `verbose` print geotiff filename


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


.. method:: pyTMD.spatial.to_geotiff(output, attributes, filename, verbose=False, varname='data', dtype=osgeo.gdal.GDT_Float64)

    Write data to a HDF5 file

    Inputs:

        `output` python dictionary of output data

        `attributes` python dictionary of output attributes

        `filename` full path of output HDF5 file


    Options:

        `verbose` print geotiff filename

        `varname` output variable name

        `dtype` GDAL data type


.. method:: pyTMD.spatial.expand_dims(obj, varname='data')

    Add a singleton dimension to a spatial dictionary if non-existent

    Options:

        variable name to modify


.. method:: pyTMD.spatial.convert_ellipsoid(phi1, h1, a1, f1, a2, f2, eps=1e-12, itmax=10)

    Convert latitudes and heights to a different ellipsoid using Newton-Raphson

    Inputs:

        `phi1`: latitude of input ellipsoid in degrees

        `h1`: height above input ellipsoid in meters

        `a1`: semi-major axis of input ellipsoid

        `f1`: flattening of input ellipsoid

        `a2`: semi-major axis of output ellipsoid

        `f2`: flattening of output ellipsoid


    Options:

        `eps`: tolerance to prevent division by small numbers and to determine convergence

        `itmax`: maximum number of iterations to use in Newton-Raphson


    Returns:

        `phi2`: latitude of output ellipsoid in degrees

        `h2`: height above output ellipsoid in meters


.. method:: pyTMD.spatial.scale_areas(lat, flat=1.0/298.257223563, ref=70.0)

    Calculates area scaling factors for a polar stereographic projection

    Inputs:

        `lat`: latitude

    Options:

        `flat`: ellipsoidal flattening

        `ref`: reference latitude (true scale latitude)

    Returns:

        `scale`: area scaling factors at input latitudes
