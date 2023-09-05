=======
spatial
=======

Utilities for reading, writing and operating on spatial data

 - Can read ascii, netCDF4, HDF5 or geotiff files
 - Can output to ascii, netCDF4, HDF5 or geotiff files

Calling Sequence
----------------

Reading a netCDF4 file

.. code-block:: python

    import pyTMD.spatial
    dinput = pyTMD.spatial.from_netCDF4(path_to_netCDF4_file)

Reading a HDF5 file

.. code-block:: python

    import pyTMD.spatial
    dinput = pyTMD.spatial.from_HDF5(path_to_HDF5_file)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/spatial.py

General Methods
===============


.. autofunction:: pyTMD.spatial.case_insensitive_filename

.. autofunction:: pyTMD.spatial.data_type

.. autofunction:: pyTMD.spatial.from_file

.. autofunction:: pyTMD.spatial.from_ascii

.. autofunction:: pyTMD.spatial.from_netCDF4

.. autofunction:: pyTMD.spatial.from_HDF5

.. autofunction:: pyTMD.spatial.from_geotiff

.. autofunction:: pyTMD.spatial.to_ascii

.. autofunction:: pyTMD.spatial.to_netCDF4

.. autofunction:: pyTMD.spatial._drift_netCDF4

.. autofunction:: pyTMD.spatial._grid_netCDF4

.. autofunction:: pyTMD.spatial._time_series_netCDF4

.. autofunction:: pyTMD.spatial.to_HDF5

.. autofunction:: pyTMD.spatial.to_geotiff

.. autofunction:: pyTMD.spatial.expand_dims

.. autofunction:: pyTMD.spatial.default_field_mapping

.. autofunction:: pyTMD.spatial.inverse_mapping

.. autofunction:: pyTMD.spatial.convert_ellipsoid

.. autofunction:: pyTMD.spatial.compute_delta_h

.. autofunction:: pyTMD.spatial.wrap_longitudes

.. autofunction:: pyTMD.spatial.to_cartesian

.. autofunction:: pyTMD.spatial.to_sphere

.. autofunction:: pyTMD.spatial.to_geodetic

.. autofunction:: pyTMD.spatial._moritz_iterative

.. autofunction:: pyTMD.spatial._bowring_iterative

.. autofunction:: pyTMD.spatial._zhu_closed_form

.. autofunction:: pyTMD.spatial.scale_areas
