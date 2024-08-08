===========================
compute_tidal_elevations.py
===========================

- Calculates tidal elevations for an input file following [Egbert2002]_
- Can use OTIS format tidal solutions provided by Oregon State University and ESR
- Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
- Can use Finite Element Solution (FES) models provided by AVISO
- Can read and write ascii, netCDF4, HDF5, parquet and geotiff formats

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tidal_elevations.py

Calling Sequence
################

.. argparse::
    :filename: compute_tidal_elevations.py
    :func: arguments
    :prog: compute_tidal_elevations.py
    :nodescription:
    :nodefault:

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

    --cutoff -c : @after
        * set to ``'inf'`` to extrapolate for all points

    --apply-flexure : @after
        Only valid for models containing flexure fields

References
##########

.. [Egbert2002] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
    Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
    Technology*, 19(2), 183--204, (2002).
    `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

.. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2

