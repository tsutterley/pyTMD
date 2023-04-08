============================
compute_SET_displacements.py
============================

- Calculates radial solid Earth tide displacements for an input file following IERS Convention (2010) guidelines

  * `https://iers-conventions.obspm.fr/chapter7.php <https://iers-conventions.obspm.fr/chapter7.php>`_
- Can read and write ascii, netCDF4, HDF5 and geotiff formats

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_SET_displacements.py

.. argparse::
    :filename: compute_SET_displacements.py
    :func: arguments
    :prog: compute_SET_displacements.py
    :nodescription:
    :nodefault:

    --variables : @after
        * for csv files: the order of the columns within the file
        * for HDF5 and netCDF4 files: time, y, x and data variable names

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

    --tide-system -p : @replace
        Permanent tide system for output values

        * ``'tide_free'``: no permanent direct and indirect tidal potentials
        * ``'mean_tide'``: permanent tidal potentials (direct and indirect)
