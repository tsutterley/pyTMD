merge_delta_time.py
==================

 - Connects to servers to download historic_deltat.data and deltat.data files
    * http://maia.usno.navy.mil/ser7/
    * ftp://cddis.nasa.gov/products/iers/
    * ftp://cddis.gsfc.nasa.gov/pub/products/iers/
 - Reads IERS Bulletin-A produced daily iers_deltat.data files
 - Creates a merged delta time file combining the historic, monthly and daily files

#### Calling Sequence
```bash
python merge_delta_time.py --verbose
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/scripts/merge_delta_time.py)

#### Command Line Options
 - `-D X`, `--directory=X`: Working data directory
 - `-V`, `--verbose`: Output information about each read file
 - `-M X`, `--mode=X`: Permission mode of output file

#### Outputs
 - `Y`: calendar year
 - `M`: calendar month
 - `D`: day of the month
 - `deltat`: delta time (UT1 - UTC)
