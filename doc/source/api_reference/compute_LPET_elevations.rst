=========================
compute_LPET_elevation.py
=========================

- Calculates long-period equilibrium tides for an input file
- Uses the summation of fifteen tidal spectral lines from [Cartwright1971]_ [Cartwright1973]_
- Can read and write ascii, netCDF4, HDF5, (cloud optimized) geotiff and (geo)parquet formats

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPET_elevations.py

Calling Sequence
################

.. argparse::
    :filename: compute_LPET_elevations.py
    :func: arguments
    :prog: compute_LPET_elevations.py
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

References
##########

.. [Cartwright1971] D. E. Cartwright and R. J. Tayler,
    "New Computations of the Tide-generating Potential,"
    *Geophysical Journal of the Royal Astronomical Society*,
    23(1), 45--73. (1971). `doi: 10.1111/j.1365-246X.1971.tb01803.x
    <https://doi.org/10.1111/j.1365-246X.1971.tb01803.x>`_
.. [Cartwright1973] D. E. Cartwright and A. C. Edden,
    "Corrected Tables of Tidal Harmonics,"
    *Geophysical Journal of the Royal Astronomical Society*,
    33(3), 253--264, (1973). `doi: 10.1111/j.1365-246X.1973.tb03420.x
    <https://doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_
