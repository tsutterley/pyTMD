=========================
compute_LPET_elevation.py
=========================

- Calculates long-period equilibrium tides for an input file (ascii, netCDF4, HDF5, geotiff)
- For netCDF4 and HDF5 files, can extract the time units from attributes
- Uses the summation of fifteen tidal spectral lines from `Cartwright and Edden, (1973) <https://doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPET_elevations.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/compute_LPET_elevations.py
    :func: arguments
    :prog: compute_LPET_elevations.py
    :nodescription:
    :nodefault:

    --variables : @after
        * for csv files: the order of the columns within the file
        * for HDF5 and netCDF4 files: time, y, x and data variable names

    --type -t : @after
        * ``'drift'``: drift buoys or satellite/airborne altimetry (time per data point)
        * ``'grid'``: spatial grids or images (single time for all data points)

    --epoch -e : @after
        * ``'days since 1858-11-17T00:00:00'`` (default Modified Julian Days)

    --deltatime -d : @after
        * can be set to ``0`` to use exact calendar date from epoch

    --standard -s : @after
        * ``'UTC'``: Coordinate Universal Time
        * ``'GPS'``: GPS Time
        * ``'LORAN'``: Long Range Navigator Time
        * ``'TAI'``: International Atomic Time

    --projection : @after
        * ``4326``: latitude and longitude coordinates on WGS84 reference ellipsoid
