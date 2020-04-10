compute_OPT_displacements.py
============================

 - Calculates radial ocean pole load tide displacements for an input csv file following IERS Convention (2010) guidelines
 - http://maia.usno.navy.mil/conventions/2010officialinfo.php
 - http://maia.usno.navy.mil/conventions/chapter7.php

#### Calling Sequence
```bash
python compute_OPT_displacements.py --directory=<path_to_directory> input_file output_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/compute_OPT_displacements.py)

#### Inputs
 1. `input_file`: input csv file with columns:
    - Modified Julian Day (days since 1858-11-17 at 00:00:00)
    - latitude: degrees
    - longitude: degrees
    - elevation (height above or below WGS84 ellipsoid)
 2. `output_file`: name of output csv file that will have columns:
    - Modified Julian Day (days since 1858-11-17 at 00:00:00)
    - latitude: degrees
    - longitude: degrees
    - radial ocean pole tide displacements

#### Command Line Options
 - `-D X`, `--directory=X`: Working data directory
 - `-M X`, `--mode=X`: Permission mode of output file
