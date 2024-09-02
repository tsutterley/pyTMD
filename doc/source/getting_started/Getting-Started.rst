===============
Getting Started
===============

This documentation is intended to explain how to compute ocean, solid Earth, load and pole tide variations using the set of ``pyTMD`` programs.
See the `background material <./Background.html>`_ and `glossary <./Glossary.html>`_ for more information on the theory and methods used in ``pyTMD``.

Tide Model Formats
##################

Ocean and load tide constituent files are available from different modeling groups in different formats.
``pyTMD`` can access the harmonic constituents for the OTIS, GOT and FES families of ocean and load tide models.
OTIS and ATLAS formatted data use  binary files to store the constituent data for either heights (``z``) or zonal and meridional transports (``u``, ``v``).
They can be either a single file containing all the constituents (compact) or multiple files each containing a single constituent.
Arctic Ocean models can be downloaded from the NSF ArcticData server using the `arcticdata_tides.py <https://github.com/tsutterley/pyTMD/blob/main/scripts/arcticdata_tides.py>`_ program.
CATS2008 can be downloaded from the US Antarctic Program (USAP) using the `usap_cats_tides.py <https://github.com/tsutterley/pyTMD/blob/main/scripts/usap_cats_tides.py>`_ program.
ATLAS netCDF formatted data use netCDF4 files for each constituent and variable type (``z``, ``u``, ``v``).
GOT formatted data use ascii files for each height constituent (``z``).
FES formatted data use either ascii (1999, 2004) or netCDF4 (2012, 2014) files for each constituent and variable type (``z``, ``u``, ``v``).
The FES models can be downloaded using the `aviso_fes_tides.py <https://github.com/tsutterley/pyTMD/blob/main/scripts/aviso_fes_tides.py>`_ program for users registered with AVISO.

Directories
###########

``pyTMD`` uses a tree structure for storing the tidal constituent data.
This structure was chosen based on the different formats of each tide model.
Presently, the following models and their directories parameterized within ``pyTMD``.

- Circum-Antarctic Tidal Simulations [Padman2008]_

    * CATS0201: ``<path_to_tide_models>/cats0201_tmd/``
    * `CATS2008 <https://doi.org/10.15784/601235>`_: ``<path_to_tide_models>/CATS2008/``
    * CATS2008_load: ``<path_to_tide_models>/CATS2008a_SPOTL_Load/``
    * `CATS2008-v2023 <https://doi.org/10.15784/601772>`_: ``<path_to_tide_models>/CATS2008_v2023/``

- Arctic Ocean and Greenland Coast Tidal Simulations [Padman2004]_

    * `AODTM-5 <https://arcticdata.io/catalog/view/doi:10.18739/A2901ZG3N>`_: ``<path_to_tide_models>/aodtm5_tmd/``
    * `AOTIM-5 <https://arcticdata.io/catalog/view/doi:10.18739/A2S17SS80>`_: ``<path_to_tide_models>/aotim5_tmd/``
    * `AOTIM-5-2018 <https://arcticdata.io/catalog/view/doi:10.18739/A21R6N14K>`_: ``<path_to_tide_models>/Arc5km2018/``
    * `Arc2kmTM <https://arcticdata.io/catalog/view/doi:10.18739/A2D21RK6K>`_: ``<path_to_tide_models>/Arc2kmTM/``
    * Gr1km-v2: ``<path_to_tide_models>/greenlandTMD_v2/``

- TOPEX/POSEIDON global tide models [Egbert2002]_

    * TPXO7.2: ``<path_to_tide_models>/TPXO7.2_tmd/``
    * TPXO7.2_load: ``<path_to_tide_models>/TPXO7.2_load/``
    * `TPXO8-atlas <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/tpxo8_atlas/``
    * `TPXO9.1 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9.1/DATA/``
    * `TPXO9-atlas <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas/``
    * `TPXO9-atlas-v2 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas_v2/``
    * `TPXO9-atlas-v3 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas_v3/``
    * `TPXO9-atlas-v4 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas_v4/``
    * `TPXO9-atlas-v5 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas_v5/``
    * `TPXO10-atlas-v2 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO10_atlas_v2/``

- Global Ocean Tide models [Ray1999]_

    * GOT4.7: ``<path_to_tide_models>/GOT4.7/grids_oceantide/``
    * GOT4.7_load: ``<path_to_tide_models>/GOT4.7/grids_loadtide/``
    * `GOT4.8 <https://earth.gsfc.nasa.gov/sites/default/files/2022-07/got4.8.tar.gz>`_: ``<path_to_tide_models>/got4.8/grids_oceantide/``
    * `GOT4.8_load <https://earth.gsfc.nasa.gov/sites/default/files/2022-07/got4.8.tar.gz>`_: ``<path_to_tide_models>/got4.8/grids_loadtide/``
    * `GOT4.10 <https://earth.gsfc.nasa.gov/sites/default/files/2022-07/got4.10c.tar.gz>`_: ``<path_to_tide_models>/GOT4.10c/grids_oceantide/``
    * `GOT4.10_load <https://earth.gsfc.nasa.gov/sites/default/files/2022-07/got4.10c.tar.gz>`_: ``<path_to_tide_models>/GOT4.10c/grids_loadtide/``
    * `GOT5.5 <https://earth.gsfc.nasa.gov/sites/default/files/2024-07/GOT5.5.tar%201.gz>`_: ``<path_to_tide_models>/GOT5.5/ocean_tides/``
    * `GOT5.5_load <https://earth.gsfc.nasa.gov/sites/default/files/2024-07/GOT5.5.tar%201.gz>`_: ``<path_to_tide_models>/GOT5.5/load_tides/``

- Finite Element Solution tide models [Lyard2020]_

    * `FES2014 <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes/description-fes2014.html>`_: ``<path_to_tide_models>/fes2014/ocean_tide/``
    * `FES2014_load <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes/description-fes2014.html>`_: ``<path_to_tide_models>/fes2014/load_tide/``
    * `FES2022 <https://doi.org/10.24400/527896/A01-2024.004>`_: ``<path_to_tide_models>/fes2022b/ocean_tide/``
    * `FES2022_load <https://doi.org/10.24400/527896/A01-2024.004>`_: ``<path_to_tide_models>/fes2022b/load_tide/``

- Empirical Ocean Tide models [HartDavis2020]_

    * `EOT20 <https://doi.org/10.17882/79489>`_: ``<path_to_tide_models>/EOT20/ocean_tides/``
    * `EOT20_load <https://doi.org/10.17882/79489>`_: ``<path_to_tide_models>/EOT20/load_tides/``

- Hamburg direct data Assimilation Methods for Tides models [Taguchi2014]_

    * `HAMTIDE11 <https://www.cen.uni-hamburg.de/en/icdc/data/ocean/hamtide.html>`_: ``<path_to_tide_models>/hamtide/``

For other tide models, the model parameters can be set with a `model definition file <./Getting-Started.html#definition-files>`_.
Note that any alternatively defined model will have to fit the file standard of a currently supported model.

Programs
########

For users wanting to compute tide corrections for use with numpy arrays or pandas dataframes
`pyTMD.compute <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute.py>`_
is the place to start.
These are a series of functions that take ``x``, ``y``, and ``time`` coordinates and
compute the corresponding tidal elevation or currents.

.. code-block:: python

    import pyTMD
    tide_h = pyTMD.compute.tide_elevations(x, y, delta_time,
        DIRECTORY=path_to_tide_models,
        MODEL='CATS2008', EPSG=3031, EPOCH=(2000,1,1,0,0,0),
        TYPE='drift', TIME='GPS',
        METHOD='spline', FILL_VALUE=np.nan)
    tide_uv = pyTMD.compute.tide_currents(x, y, delta_time,
        DIRECTORY=path_to_tide_models,
        MODEL='CATS2008', EPSG=3031, EPOCH=(2000,1,1,0,0,0),
        TYPE='drift', TIME='GPS',
        METHOD='spline', FILL_VALUE=np.nan)


For users wanting to calculate tidal elevations or currents for a series of files, the
`compute_tidal_elevations.py <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tidal_elevations.py>`_ and
`compute_tidal_currents.py <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tidal_currents.py>`_ programs
cover most use cases.  They take an input file (in csv, netCDF4, HDF5, parquet or geotiff formats) and compute the tidal
elevations or currents (zonal and meridional) for each point.

.. code-block:: bash

    python compute_tidal_elevations.py --directory <path_to_tide_models> --tide CATS2008 \
        --format HDF5 --variables t_sec lat lon h_cor --projection 4326 \
        --epoch 'seconds since 1970-01-01T00:00:00' --verbose --mode 0o775 \
        input_file.H5 output_file.H5

    python compute_tidal_elevations.py --directory <path_to_tide_models> --tide CATS2008 \
        --format geotiff --projection 3031 --type grid --epoch '2000-01-01T12:00:00' \
        --verbose --mode 0o775 input_file.tif output_file.tif

    python compute_tidal_currents.py --directory <path_to_tide_models> --tide CATS2008 \
        --format HDF5 --variables t_sec lat lon h_cor --projection 4326 \
        --epoch 'seconds since 1970-01-01T00:00:00' --verbose --mode 0o775 \
        input_file.H5 output_file.H5

Definition Files
################

For models not parameterized within ``pyTMD``, the parameters can be set with a model definition file.
There are two types of definition files currently supported by ``pyTMD`` (``ascii`` and ``json``).

The tab-delimited text format (``ascii``) consists of two separate columns (parameter names and parameter values).
These definition files are read line by line to fill a python dictionary mapping the names with their values.
For parameters consisting of lists, the parameter values can be separated by spaces or commas.
For FES-type models of currents, the two lists of model files (``u`` and ``v``) should be separated in the list with a semi-colon.

The JSON format (``json``) directly maps the parameter names with their values stored in the appropriate data type (strings, lists, numbers, booleans, etc).
For FES-type models of currents, the two lists of model files (``u`` and ``v``) are stored in a name-value pair objects (similar to a python dictionary).
While still human readable, the JSON format is both interoperable and more easily machine readable.

Each definition file regardless of the format should have ``name``, ``format`` and ``type`` parameters.
Each model type also requires specific sets of parameters for the model reader.
For models with multiple constituent files, the files can be found using a ``glob`` string to search a directory.

- ``OTIS``, ``ATLAS-compact`` and ``TMD3``

    * ``format``: ``OTIS``, ``ATLAS-compact`` or ``TMD3``
    * ``grid_file``: path to model grid file
    * ``model_file``: path to model constituent file(s) or a ``glob`` string
    * ``name``: tide model name
    * ``projection``: `model spatial projection <./Getting-Started.html#spatial-coordinates>`_.
    * ``type``: ``z`` or ``u,v``

- ``ATLAS-netcdf``

    * ``compressed``: model files are gzip compressed
    * ``format``: ``ATLAS-netcdf``
    * ``grid_file``: path to model grid file
    * ``model_file``: path to model constituent files or a ``glob`` string
    * ``name``: tide model name
    * ``scale``: scaling factor for converting to output units
    * ``type``: ``z`` or ``u,v``

- ``GOT-ascii`` and ``GOT-netcdf``

    * ``compressed``: model files are gzip compressed
    * ``format``: ``GOT-ascii`` or ``GOT-netcdf``
    * ``model_file``: path to model constituent files or a ``glob`` string
    * ``name``: tide model name
    * ``scale``: scaling factor for converting to output units
    * ``type``: ``z``

- ``FES-ascii`` and ``FES-netcdf``

    * ``compressed``: model files are gzip compressed
    * ``format``: ``FES-ascii`` or ``FES-netcdf``
    * ``model_file``: path to model constituent files or a ``glob`` string
    * ``name``: tide model name
    * ``scale``: scaling factor for converting to output units
    * ``type``: ``z`` or ``u,v``
    * ``version``: tide model version

Time
####

The default time in ``pyTMD`` is days (UTC) since a given epoch.
For ocean, load and equilibrium tide programs, the epoch is 1992-01-01T00:00:00.
For pole tide programs, the epoch is 1858-11-17T00:00:00 (Modified Julian Days).
``pyTMD`` uses the ``timescale`` library to convert different time formats to the necessary time format of a given program.
``timescale`` can also parse date strings describing the units and epoch of relative times, or the calendar date of measurement for geotiff formats.
``timescale`` keeps updated `tables of leap seconds <https://github.com/tsutterley/timescale/blob/main/timescale/data/leap-seconds.list>`_ for converting from GPS, LORAN and TAI times.

- TAI time: International Atomic Time which is computed as the weighted average of several hundred atomic clocks.
- UTC time: Coordinated Universal Time which is `periodically adjusted <https://www.nist.gov/pml/time-and-frequency-division/leap-seconds-faqs>`_ to account for the difference between the definition of the second and the rotation of Earth.
- GPS time: Atomic timing system for the Global Positioning System constellation of satellites monitored by the United States Naval Observatory (USNO). GPS time and UTC time were equal on January 6, 1980. TAI time is ahead of GPS time by 19 seconds.
- LORAN time: Atomic timing system for the Loran-C chain transmitter sites used in terrestrial radionavigation. LORAN time and UTC time were equal on January 1, 1958. TAI time is ahead of LORAN time by 10 seconds.

``timescale`` also keeps updated `tables of delta times <https://github.com/tsutterley/timescale/blob/main/timescale/data/merged_deltat.data>`_ for converting between dynamic (TT) and universal (UT1) times.
Delta times (TT - UT1) are the differences between Dynamic Time (TT) and Universal Time (UT1) [Meeus1998]_.
Universal Time (UT1) is based on the rotation of the Earth,
which varies irregularly, and so UT1 is adjusted periodically.
Dynamic Time (TT) is a uniform, monotonically increasing time standard based on atomic clocks that is
used for the accurate calculation of celestial mechanics, orbits and ephemerides.
Delta times can be added to Universal Time (UT1) values to convert to Dynamic Time (TT) values.

Spatial Coordinates
###################

The default coordinate system in ``pyTMD`` is WGS84 geodetic coordinates in latitude and longitude.
``pyTMD`` uses `pyproj <https://pypi.org/project/pyproj/>`_ to convert from different coordinate systems and datums.
Some regional tide models are projected in a different coordinate system.
For these cases, ``pyTMD`` will `convert from latitude and longitude to the model coordinate system <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/convert_crs.py>`_.

OTIS models may be projected into a separate coordinate system. The available OTIS projections within ``pyTMD`` are

- ``'4326'`` (global latitude and longitude)
- ``'3031'`` (Antarctic Polar Stereographic in kilometers)
- ``'3413'`` (NSIDC Sea Ice Polar Stereographic North in kilometers)
- ``'CATS2008'`` (CATS2008 Polar Stereographic in kilometers)
- ``'3976'`` (NSIDC Sea Ice Polar Stereographic South in kilometers)
- ``'PSNorth'`` (idealized polar stereographic in kilometers)

For other model projections, a formatted coordinate reference system (CRS) descriptor (e.g. PROJ, WKT, or EPSG code) can be used.

Interpolation
#############

For converting from model coordinates, ``pyTMD`` uses spatial interpolation routines from `scipy <https://docs.scipy.org/doc/scipy/reference/interpolate.html>`_
along with a built-in `bilinear <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/interpolate.py>`_ interpolation routine.
The default interpolator uses a `biharmonic spline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html>`_
function to interpolate from the model coordinate system to the output coordinates.
There are options to use nearest and linear interpolators with the
`regular grid <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html>`_ function.
For coastal or near-grounded points, the model can be extrapolated using a
`nearest-neighbor <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/interpolate.py>`_ routine.
The default maximum extrapolation distance is 10 kilometers.
This default distance may not be a large enough extrapolation for some applications and models.
The extrapolation cutoff can be set to any distance in kilometers, but should be used with caution in cases such as narrow fjords or ice sheet grounding zones [Padman2018]_.

References
##########

.. [Egbert2002] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of Barotropic Ocean Tides", *Journal of Atmospheric and Oceanic Technology*, 19(2), 183--204, (2002). `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

.. [HartDavis2020] M. G. Hart-Davis, G. Piccioni, D. Dettmering, C. Schwatke, M. Passaro, and F. Seitz, "EOT20: a global ocean tide model from multi-mission satellite altimetry", *Earth System Science Data*, 13(8), 3869--3884, (2020). `doi: 10.5194/essd-13-3869-2021 <https://doi.org/10.5194/essd-13-3869-2021>`_

.. [Lyard2020] F. H. Lyard, D. J. Allain, M. Cancet, L. Carr\ |egrave|\ re, and N. Picot, "FES2014 global ocean tides atlas: design and performances", *Ocean Science Discussions*, in review, (2020). `doi: 10.5194/os-2020-96 <https://doi.org/10.5194/os-2020-96>`_

.. [Meeus1998] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).

.. [Padman2004] L. Padman and S. Y. Erofeeva, "A barotropic inverse tidal model for the Arctic Ocean", *Geophysical Research Letters*, 31(2), L02303. (2004). `doi: 10.1029/2003GL019003 <https://doi.org/10.1029/2003GL019003>`_

.. [Padman2008] L. Padman, S. Y. Erofeeva, and H. A. Fricker, "Improving Antarctic tide models by assimilation of ICESat laser altimetry over ice shelves", *Geophysical Research Letters*, 35, L22504, (2008). `doi: 10.1029/2008GL035592 <https://doi.org/10.1029/2008GL035592>`_

.. [Padman2018] L. Padman, M. R. Siegfried, and H. A. Fricker, "Ocean Tide Inﬂuences on the Antarctic and Greenland Ice Sheets", *Reviews of Geophysics*, 56, (2018). `doi: 10.1002/2016RG000546 <https://doi.org/10.1002/2016RG000546>`_

.. [Ray1999] R. D. Ray, "A Global Ocean Tide Model From TOPEX/POSEIDON Altimetry: GOT99.2", *NASA Technical Memorandum*, `NASA/TM--1999-209478 <https://ntrs.nasa.gov/citations/19990089548>`_.

.. [Taguchi2014] E. Taguchi, D. Stammer, and W. Zahel, "Inferring deep ocean tidal energy dissipation from the global high-resolution data-assimilative HAMTIDE model", *Journal of Geophysical Research: Oceans*, 119, 4573--4592, (2014). `doi: 10.1002/2013JC009766 <https://doi.org/10.1002/2013JC009766>`_

.. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2

.. |egrave|    unicode:: U+00E8 .. LATIN SMALL LETTER E WITH GRAVE
