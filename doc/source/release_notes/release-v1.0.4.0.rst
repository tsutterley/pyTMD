####################
`Release v1.0.4.0`__
####################

- ``refactor``: use model class for files and attributes
- ``feat``: can use numpy datetime arrays as input time variable
- ``feat``: calculator for height differences between ellipsoids
- ``feat``: can use argparse prefix files to define command line arguments
- ``feat``: added generic list from Apache http server
- ``fix``: netCDF4 cases where there is no mask on constituent files
- ``fix``: USAP now requires capchas
- ``fix``: pole case in stereographic area scale calculation
- ``fix``: use ``getncattr`` to get attributes from netCDF4 files to prevent deprecation errors
- ``ci``: install proj from source for cartopy dependency
- ``test``: NSIDC no longer requires authentication headers
- ``test``: added test for model definition files

.. __: https://github.com/tsutterley/pyTMD/releases/tag/1.0.4.0
