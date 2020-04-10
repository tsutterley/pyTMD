compute_tidal_elevations.py
===========================

- Calculates tidal elevations for an input csv file
- Can use OTIS format tidal solutions provided by Ohio State University and ESR
- Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

#### Calling Sequence
```bash
python compute_tidal_elevations.py --directory=<path_to_directory> --tide=<model> input_file output_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/compute_tidal_elevations.py)

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
   - tide heights

#### Command Line Options
 - `-D X`, `--directory=X`: Working data directory
 - `-T X`, `--tide=X`: Tide model to use in correction
     * CATS0201
     * CATS2008
     * CATS2008_load
     * TPXO9-atlas
     * TPXO9.1
     * TPXO8-atlas
     * TPXO7.2
     * TPXO7.2_load
     * AODTM-5
     * AOTIM-5
     * AOTIM-5-2018
     * GOT4.7
     * GOT4.7_load
     * GOT4.8
     * GOT4.8_load
     * GOT4.10
     * GOT4.10_load
 - `-M X`, `--mode=X`: Permission mode of output file
 - `-V`, `--verbose`: Output information about each created file
