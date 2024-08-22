============================
compute_OPT_displacements.py
============================

- Calculates radial ocean pole load tide displacements for an input file following IERS Convention (2010) guidelines

  * `https://iers-conventions.obspm.fr/chapter7.php <https://iers-conventions.obspm.fr/chapter7.php>`_
  * `ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz <ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz>`_
- Can read and write ascii, netCDF4, HDF5, (cloud optimized) geotiff and (geo)parquet formats

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_OPT_displacements.py

Calling Sequence
################

.. argparse::
    :filename: compute_OPT_displacements.py
    :func: arguments
    :prog: compute_OPT_displacements.py
    :nodescription:
    :nodefault:

    --format -F : @after
        * ``'csv'``: ascii with comma-separated values
        * ``'netCDF4'``: Network Common Data Format version 4
        * ``'HDF5'``: Hierarchical Data Format version 5
        * ``'GTiff'``: GeoTiff georeferenced raster imagery
        * ``'cog'``: Cloud Optimized GeoTIFF
        * ``'parquet'``: Apache (geo)parquet

    --variables : @after
        * for csv files: the order of the columns within the file
        * for HDF5, netCDF4 and parquet files: time, y, x and data variable names

    --type -t : @after
        * ``'drift'``: drift buoys or satellite/airborne altimetry (time per data point)
        * ``'grid'``: spatial grids or images (single time for all data points)
        * ``'time series'``: station locations with multiple time values

    --epoch -e : @after
        * ``'days since 1858-11-17T00:00:00'`` (default Modified Julian Days)

    --deltatime -d : @after
        * can be set to ``0`` to use exact calendar date from epoch

    --standard -s : @after
        * ``'UTC'``: Coordinate Universal Time
        * ``'GPS'``: GPS Time
        * ``'LORAN'``: Long Range Navigator Time
        * ``'TAI'``: International Atomic Time
        * ``'datetime'``: formatted datetime string in UTC

    --projection : @after
        * ``4326``: latitude and longitude coordinates on WGS84 reference ellipsoid
