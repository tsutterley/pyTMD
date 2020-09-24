compute_OPT_displacements.py
============================

 - Calculates radial ocean pole load tide displacements for an input file following IERS Convention (2010) guidelines
 - http://maia.usno.navy.mil/conventions/2010officialinfo.php
 - http://maia.usno.navy.mil/conventions/chapter7.php

#### Calling Sequence
```bash
python compute_OPT_displacements.py --directory=<path_to_directory> input_file output_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_OPT_displacements.py)

#### Inputs
 1. `input_file`: name of input file
 2. `output_file`: name of output file

#### Command Line Options
 - `-D X`, `--directory=X`: Working data directory
 - `--format=X`: input and output data format
     * `'csv'` (default)
     * `'netCDF4'`
     * `'HDF5'`
 - `--variables=X`: variable names of data in csv, HDF5 or netCDF4 file
     * for csv files: the order of the columns within the file
     * for HDF5 and netCDF4 files: time, y, x and data variable names
 - `--epoch=X`: Reference epoch of input time
     * `'days since 1858-11-17T00:00:00'` (default Modified Julian Days)
 - `--projection=X`: spatial projection as EPSG code or PROJ4 string
     * `4326`: latitude and longitude coordinates on WGS84 reference ellipsoid
 - `-V`, `--verbose`: Verbose output of processing run
 - `-M X`, `--mode=X`: Permission mode of output file
