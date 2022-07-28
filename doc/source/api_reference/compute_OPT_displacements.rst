============================
compute_OPT_displacements.py
============================

- Calculates radial ocean pole load tide displacements for an input file following IERS Convention (2010) guidelines
- Can read and write ascii, netCDF4, HDF5 and geotiff formats
- `http://maia.usno.navy.mil/conventions/2010officialinfo.php <http://maia.usno.navy.mil/conventions/2010officialinfo.php>`_
- `http://maia.usno.navy.mil/conventions/chapter7.php <http://maia.usno.navy.mil/conventions/chapter7.php>`_
- `ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz <ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz>`_

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_OPT_displacements.py

.. argparse::
    :filename: ../../scripts/compute_OPT_displacements.py
    :func: arguments
    :prog: compute_OPT_displacements.py
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