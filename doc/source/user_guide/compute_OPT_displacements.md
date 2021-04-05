compute_OPT_displacements.py
============================

 - Calculates radial ocean pole load tide displacements for an input file following IERS Convention (2010) guidelines
 - Can read and write ascii, netCDF4, HDF5 and geotiff formats
 - [http://maia.usno.navy.mil/conventions/2010officialinfo.php](http://maia.usno.navy.mil/conventions/2010officialinfo.php)
 - [http://maia.usno.navy.mil/conventions/chapter7.php](http://maia.usno.navy.mil/conventions/chapter7.php)
 - [ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz](ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz)

#### Calling Sequence
```bash
python compute_OPT_displacements.py input_file output_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_OPT_displacements.py)

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
 - `-I X`, `--interpolate X`: Interpolation method
     * `'spline'`
     * `'linear'`
     * `'nearest'`
 - `-V`, `--verbose`: Verbose output of processing run
 - `-M X`, `--mode X`: Permission mode of output file
