iers_delta_time.py
==================

 - Connects to the [IERS ftp server](ftp://ftp.iers.org/products/eop/rapid/bulletina/) to download Bulletin-A files
 - Reads the IERS Bulletin-A files and calculates the daily delta times
 - Delta times are the difference between universal time and dynamical time

#### Calling Sequence
```bash
python iers_delta_time.py --verbose
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/scripts/iers_delta_time.py)

#### Command Line Options
 - `-D X`, `--directory=X`: Working data directory
 - `-V`, `--verbose`: Output information about each read file
 - `-M X`, `--mode=X`: Permission mode of output file

#### Outputs
 - `Y`: calendar year
 - `M`: calendar month
 - `D`: day of the month  
 - `deltat`: delta time (UT1 - UTC)
