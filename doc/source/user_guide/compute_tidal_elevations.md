compute_tidal_elevations.py
===========================

- Calculates tidal elevations for an input file (ascii, netCDF4, HDF5, geotiff)
- For netCDF4 and HDF5 files, can extract the time units from attributes
- Can use OTIS format tidal solutions provided by Ohio State University and ESR
- Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
- Can use Finite Element Solution (FES) models provided by AVISO

#### Calling Sequence
```bash
python compute_tidal_elevations.py --directory <path_to_directory> --tide <model> input_file output_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tidal_elevations.py)

#### Inputs
 1. `input_file`: name of input file
 2. `output_file`: name of output file

#### Command Line Options
 - `-D X`, `--directory X`: Working data directory
 - `-T X`, `--tide X`: Tide model to use in correction
     * `'CATS0201'`
     * `'CATS2008'`
     * `'CATS2008_load'`
     * `'TPXO9-atlas'`
     * `'TPXO9-atlas-v2'`
     * `'TPXO9-atlas-v3'`
     * `'TPXO9-atlas-v4'`
     * `'TPXO9.1'`
     * `'TPXO8-atlas'`
     * `'TPXO7.2'`
     * `'TPXO7.2_load'`
     * `'AODTM-5'`
     * `'AOTIM-5'`
     * `'AOTIM-5-2018'`
     * `'GOT4.7'`
     * `'GOT4.7_load'`
     * `'GOT4.8'`
     * `'GOT4.8_load'`
     * `'GOT4.10'`
     * `'GOT4.10_load'`
     * `'FES2014'`
     * `'FES2014_load'`
 - `--format X`: input and output data format
     * `'csv'` (default)
     * `'netCDF4'`
     * `'HDF5'`
     * `'geotiff'`
 - `--variables X`: variable names of data in csv, HDF5 or netCDF4 file
     * for csv files: the order of the columns within the file
     * for HDF5 and netCDF4 files: time, y, x and data variable names
 - `-H X`, `--header X`: number of header lines for csv files
 - `-t X`, `--type X`: input data type
     * `'drift'`: drift buoys or satellite/airborne altimetry (time per data point)
     * `'grid'`: spatial grids or images (single time for all data points)
 - `-e X`, `--epoch X`: reference epoch of input time or calendar date of measurement
     * `'days since 1858-11-17T00:00:00'` (default Modified Julian Days)
 - `-d X`, `--deltatime X`: input delta time for files without date information
     * can be set to 0 to use exact calendar date from epoch
 - `--projection X`: spatial projection as EPSG code or PROJ4 string
     * `4326`: latitude and longitude coordinates on WGS84 reference ellipsoid
 - `-I X`, `--interpolate X`: Interpolation method
     * `'spline'`
     * `'linear'`
     * `'nearest'`
     * `'bilinear'`
 - `-E X`, `--extrapolate X`: Extrapolate with nearest-neighbors
 - `-V`, `--verbose`: Verbose output of processing run
 - `-M X`, `--mode X`: Permission mode of output file
