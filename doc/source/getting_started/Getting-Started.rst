===============
Getting Started
===============

This documentation is intended to explain how to compute ocean, load and pole tide variations using the set of ``pyTMD`` programs.
The rise and fall of the oceanic tides are a major source of the vertical variability of the ocean surface.
Ocean tides are typically observed using float gauges, GPS stations, pressure recorders, and satellite altimetry.
Ocean and load tides are driven by gravitational undulations due to the relative positions of the Earth, moon and sun, and the centripetal acceleration due to the Earth's rotation.
The tidal oscillations can be decomposed into a series of tidal constituents (or partial tides) of particular frequencies.
Ocean and load tide constituent files are available from different modeling groups in different formats.
``pyTMD`` can access the harmonic constituents for the OTIS, GOT and FES families of ocean and load tide models.
These tide models will be one of following categories depending on the version: 1) an empirically adjusted model,
2) a barotropic hydrodynamic model constrained by data assimilation,
or 3) an unconstrained hydrodynamic model [Stammer2014]_.

Load and ocean pole tides are driven by variations in the Earth's figure axis.
These pole tides are due to Earth's ellipsoidal shape shifting as the rotation axis of the Earth
moves with respect to the mean pole location, and for the case of ocean pole tides the centripetal effects of polar motion on the ocean.
The Earth Orientation Parameters (EOPs) necessary for computing load pole and ocean pole tide variations are included within the ``pyTMD`` program.

Tide Model Formats
##################

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

- Circum-Antarctic Tidal Simulations [Padman2008]_
    * CATS0201: ``<path_to_tide_models>/cats0201_tmd/``
    * `CATS2008 <https://www.usap-dc.org/view/dataset/601235>`_: ``<path_to_tide_models>/CATS2008/``
    * CATS2008_load: ``<path_to_tide_models>/CATS2008a_SPOTL_Load/``

- Arctic Ocean Tidal Simulations [Padman2004]_
    * `AODTM-5 <https://arcticdata.io/catalog/view/doi:10.18739/A2901ZG3N>`_: ``<path_to_tide_models>/aodtm5_tmd/``
    * `AOTIM-5 <https://arcticdata.io/catalog/view/doi:10.18739/A2S17SS80>`_: ``<path_to_tide_models>/aotim5_tmd/``
    * `AOTIM-5-2018 <https://arcticdata.io/catalog/view/doi:10.18739/A21R6N14K>`_: ``<path_to_tide_models>/Arc5km2018/``

- TOPEX/POSEIDON global tide models [Egbert2002]_
    * `TPXO9-atlas <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas/``
    * `TPXO9-atlas-v2 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas_v2/``
    * `TPXO9-atlas-v3 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas_v3/``
    * `TPXO9-atlas-v4 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9_atlas_v4/``
    * `TPXO9.1 <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/TPXO9.1/DATA/``
    * `TPXO8-atlas <https://www.tpxo.net/tpxo-products-and-registration>`_: ``<path_to_tide_models>/tpxo8_atlas/``
    * TPXO7.2: ``<path_to_tide_models>/TPXO7.2_tmd/``
    * TPXO7.2_load: ``<path_to_tide_models>/TPXO7.2_load/``

- Global Ocean Tide models [Ray1999]_
    * GOT4.7: ``<path_to_tide_models>/GOT4.7/grids_oceantide/``
    * GOT4.7_load: ``<path_to_tide_models>/GOT4.7/grids_loadtide/``
    * GOT4.8: ``<path_to_tide_models>/got4.8/grids_oceantide/``
    * GOT4.8_load: ``<path_to_tide_models>/got4.8/grids_loadtide/``
    * GOT4.10: ``<path_to_tide_models>/GOT4.10c/grids_oceantide/``
    * GOT4.10_load: ``<path_to_tide_models>/GOT4.10c/grids_loadtide/``

- Finite Element Solution tide models [Lyard2020]_
    * `FES2014 <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes/description-fes2014.html>`_: ``<path_to_tide_models>/fes2014/ocean_tide/``
    * `FES2014_load <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes/description-fes2014.html>`_: ``<path_to_tide_models>/fes2014/load_tide/``

Programs
########

For users wanting to compute tide corrections for use with numpy arrays or pandas dataframes
`compute_tide_corrections.py <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute_tide_corrections.py>`_
is the place to start.  It is a function that takes ``x``, ``y``, and ``time`` coordinates and
computes the corresponding tidal elevation.

.. code-block:: python

    tide = compute_tide_corrections(x, y, delta_time, DIRECTORY=path_to_tide_models,
        MODEL='CATS2008', EPSG=3031, EPOCH=(2000,1,1,0,0,0), TYPE='drift', TIME='GPS',
        METHOD='spline', FILL_VALUE=np.nan)


For users wanting to calculate tidal elevations or currents for a series of files, the
`compute_tidal_elevations.py <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tidal_elevations.py>`_ and
`compute_tidal_currents.py <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tidal_currents.py>`_ programs
cover most use cases.  They take an input file (in csv, netCDF4, HDF5 or geotiff formats) and compute the tidal
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


There are specific programs for correcting some publicly available elevation datasets:

- `NASA Operation IceBridge data <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_icebridge_data.py>`_
- `ICESat GLA12 ice sheet altimetry data <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat_GLA12.py>`_
- `ICESat-2 ATL03 geolocated photon data <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL03.py>`_
- `ICESat-2 ATL06 land ice height data <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL06.py>`_
- `ICESat-2 ATL07 sea ice height data <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL07.py>`_
- `ICESat-2 ATL11 annual land ice height data <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL11.py>`_
- `ICESat-2 ATL12 ocean surface height data <https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL12.py>`_

Time
####

The default time in ``pyTMD`` is days (UTC) since a given epoch.
For ocean, load and equilibrium tide programs, the epoch is 1992-01-01T00:00:00.
For pole tide programs, the epoch is 1858-11-17T00:00:00 (Modified Julian Days).
The `time module <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/time.py>`_ within ``pyTMD`` can convert different time formats to the necessary time format of a given program.
The `time module <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/time.py>`_ can also parse date strings describing the units and epoch of relative times, or the calendar date of measurement for geotiff formats.
``pyTMD`` keeps updated `tables of leap seconds <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/data/leap-seconds.list>`_ for converting from GPS and TAI times.
``pyTMD`` keeps updated `tables of delta times <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/data/merged_deltat.data>`_ for converting between dynamic (TT) and universal (UT1) times.

Spatial Coordinates
###################

The default coordinate system in ``pyTMD`` is WGS84 geodetic coordinates in latitude and longitude.
``pyTMD`` uses `pyproj <https://pypi.org/project/pyproj/>`_ to convert from different coordinate systems and datums.
Some regional tide models are projected in a different coordinate system.
For these cases, ``pyTMD`` will `convert from latitude and longitude to the model coordinate system <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/convert_ll_xy.py>`_.

Interpolation
#############

For converting from model coordinates, ``pyTMD`` uses spatial interpolation routines from `scipy <https://docs.scipy.org/doc/scipy/reference/interpolate.html>`_
along with a built-in `bilinear <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/bilinear_interp.py>`_ interpolation routine.
The default interpolator uses a `biharmonic spline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html>`_
function to interpolate from the model coordinate system to the output coordinates.
There are options to use nearest and linear interpolators with the
`regular grid <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html>`_ function.
For coastal or near-grounded points, the model can be extrapolated using a
`nearest-neighbor <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/nearest_extrap.py>`_ routine.

References
##########

.. [Egbert2002] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of Barotropic Ocean Tides", *Journal of Atmospheric and Oceanic Technology*, 19(2), 183--204, (2002). `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

.. [Lyard2020] F. H. Lyard, D. J. Allain, M. Cancet, L. Carr\ |egrave|\ re, and N. Picot, "FES2014 global ocean tides atlas: design and performances", *Ocean Science Discussions*, in review, (2020). `doi: 10.5194/os-2020-96 <https://doi.org/10.5194/os-2020-96>`_

.. [Padman2004] L. Padman and S. Y. Erofeeva, "A barotropic inverse tidal model for the Arctic Ocean", *Geophysical Research Letters*, 31(2), L02303. (2004). `doi: 10.1029/2003GL019003 <https://doi.org/10.1029/2003GL019003>`_

.. [Padman2008] L. Padman, S. Y. Erofeeva, and H. A. Fricker, "Improving Antarctic tide models by assimilation of ICESat laser altimetry over ice shelves", *Geophysical Research Letters*, 35, L22504, (2008). `doi:10.1029/2008GL035592 <https://doi.org/10.1029/2008GL035592>`_

.. [Ray1999] R. D. Ray, "A Global Ocean Tide Model From TOPEX/POSEIDON Altimetry: GOT99.2", *NASA Technical Memorandum*, `NASA/TM--1999-209478 <https://ntrs.nasa.gov/search.jsp?R=19990089548>`_.

.. [Stammer2014] D. Stammer et al., "Accuracy assessment of global barotropic ocean tide models", *Reviews of Geophysics*, 52, 243--282, (2014). `doi:10.1002/2014RG000450 <https://doi.org/10.1002/2014RG000450>`_

.. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2

.. |egrave|    unicode:: U+00E8 .. LATIN SMALL LETTER E WITH GRAVE
