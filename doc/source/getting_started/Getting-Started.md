Getting Started
===============

#### Tide Model Formats
OTIS and ATLAS formatted data use a single binary file to store all the constituents for either heights (`z`) or transports (`u`, `v`).
Arctic Ocean models can be downloaded from the NSF ArcticData server using the [`arcticdata_tides.py`](https://github.com/tsutterley/pyTMD/blob/main/scripts/arcticdata_tides.py) program.
CATS2008 can be downloaded from the US Antarctic Program (USAP) using the [`usap_cats_tides.py`](https://github.com/tsutterley/pyTMD/blob/main/scripts/usap_cats_tides.py) program.
ATLAS netCDF formatted data use netCDF4 files for each constituent and variable type (`z`, `u`, `v`).
GOT formatted data use ascii files for each height constituent (`z`).
FES formatted data use either ascii (1999, 2004) or netCDF4 (2012, 2014) files for each constituent and variable type (`z`, `u`, `v`).
The FES models can be downloaded using the [`aviso_fes_tides.py`](https://github.com/tsutterley/pyTMD/blob/main/scripts/aviso_fes_tides.py) program for users registered with AVISO.

#### Directories
pyTMD uses a tree structure for storing the tidal constituent data.
This structure was chosen based on the different formats of each tide model.

- Circum-Antarctic Tidal Simulations
    * CATS0201: `<path_to_tide_models>/cats0201_tmd/`
    * [CATS2008](https://www.usap-dc.org/view/dataset/601235): `<path_to_tide_models>/CATS2008/`
    * CATS2008_load: `<path_to_tide_models>/CATS2008a_SPOTL_Load/`

- Arctic Ocean Tidal Simulations
    * [AODTM-5](https://arcticdata.io/catalog/view/doi:10.18739/A2901ZG3N): `<path_to_tide_models>/aodtm5_tmd/`
    * [AOTIM-5](https://arcticdata.io/catalog/view/doi:10.18739/A2S17SS80): `<path_to_tide_models>/aotim5_tmd/`
    * [AOTIM-5-2018](https://arcticdata.io/catalog/view/doi:10.18739/A21R6N14K): `<path_to_tide_models>/Arc5km2018/`

- TOPEX/POSEIDON global tide models
    * [TPXO9-atlas](https://www.tpxo.net/tpxo-products-and-registration): `<path_to_tide_models>/TPXO9_atlas/`
    * [TPXO9-atlas-v2](https://www.tpxo.net/tpxo-products-and-registration): `<path_to_tide_models>/TPXO9_atlas_v2/`
    * [TPXO9-atlas-v3](https://www.tpxo.net/tpxo-products-and-registration): `<path_to_tide_models>/TPXO9_atlas_v3/`
    * [TPXO9.1](https://www.tpxo.net/tpxo-products-and-registration): `<path_to_tide_models>/TPXO9.1/DATA/`
    * [TPXO8-atlas](https://www.tpxo.net/tpxo-products-and-registration): `<path_to_tide_models>/tpxo8_atlas/`
    * TPXO7.2: `<path_to_tide_models>/TPXO7.2_tmd/`
    * TPXO7.2_load: `<path_to_tide_models>/TPXO7.2_load/`

- Global Ocean Tide models
    * GOT4.7: `<path_to_tide_models>/GOT4.7/grids_oceantide/`
    * GOT4.7_load: `<path_to_tide_models>/GOT4.7/grids_loadtide/`
    * GOT4.8: `<path_to_tide_models>/got4.8/grids_oceantide/`
    * GOT4.8_load: `<path_to_tide_models>/got4.8/grids_loadtide/`
    * GOT4.10: `<path_to_tide_models>/GOT4.10c/grids_oceantide/`
    * GOT4.10_load: `<path_to_tide_models>/GOT4.10c/grids_loadtide/`

- Finite Element Solution tide models
    * [FES2014](https://www.aviso.altimetry.fr/data/products/auxiliary-products/global-tide-fes.html): `<path_to_tide_models>/fes2014/ocean_tide/`
    * [FES2014_load](https://www.aviso.altimetry.fr/data/products/auxiliary-products/global-tide-fes.html): `<path_to_tide_models>/fes2014/load_tide/`

#### Programs
For users wanting to compute tide corrections for use with numpy arrays or pandas data structures, [`compute_tide_corrections.py`](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute_tide_corrections.py) is the place to start.  It is a function that takes x, y, and time coordinates and computes the corresponding tidal elevation.
```python
tide = compute_tide_corrections(x, y, delta_time, DIRECTORY=path_to_tide_models,
    MODEL='CATS2008', EPSG=3031, EPOCH=(2000,1,1,0,0,0), TYPE='drift', TIME='GPS',
    METHOD='spline', FILL_VALUE=np.nan)
```

For users wanting to calculate tidal elevations or currents for a series of files, the [`compute_tidal_elevations.py`](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tidal_elevations.py) and [`compute_tidal_currents.py`](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tidal_currents.py) programs cover most use cases.  They take an input file (in csv, netCDF4 or HDF5) and compute the tidal elevations or currents (zonal and meridonal) for each point.
```bash
python compute_tidal_elevations.py --directory <path_to_tide_models> --tide CATS2008 \
    --format HDF5 --variables t_sec lat lon h_cor --projection 4326 \
    --epoch 'seconds since 1970-01-01T00:00:00' --verbose --mode 0o775 \
    input_file.H5 output_file.H5

python compute_tidal_currents.py --directory <path_to_tide_models> --tide CATS2008 \
    --format HDF5 --variables t_sec lat lon h_cor --projection 4326 \
    --epoch 'seconds since 1970-01-01T00:00:00' --verbose --mode 0o775 \
    input_file.H5 output_file.H5
```

There are specific programs for correcting NASA [Operation IceBridge](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_icebridge_data.py), [ICESat-2 ATL03 geolocated photon](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL03.py), [ICESat-2 ATL06 land ice](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL06.py), [ICESat-2 ATL07 sea ice](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL07.py) and [ICESat-2 ATL12 ocean surface](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL12.py) data.
