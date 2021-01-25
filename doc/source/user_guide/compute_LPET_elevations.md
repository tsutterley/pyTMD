compute_LPET_elevation.py
=========================

 - Calculates long-period equilibrium tides for an input file (ascii, netCDF4, HDF5, geotiff)
 - For netCDF4 and HDF5 files, can extract the time units from attributes
 - Uses the summation of fifteen tidal spectral lines from [Cartwright and Edden, (1973)](https://doi.org/10.1111/j.1365-246X.1973.tb03420.x)

#### Calling Sequence
```bash
python compute_LPET_elevation.py input_file output_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPET_elevation.py)

#### Inputs
 1. `input_file`: name of input file
 2. `output_file`: name of output file

#### Command Line Options
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
- `-V`, `--verbose`: Verbose output of processing run
 - `-M X`, `--mode X`: Permission mode of output file
