compute_tides_ICESat_GLA12.py
=============================

- Calculates tidal elevations for correcting ICESat/GLAS L2 GLA12 Antarctic and Greenland Ice Sheet elevation data
- Can use OTIS format tidal solutions provided by Ohio State University and ESR
- Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
- Can use Finite Element Solution (FES) models provided by AVISO

#### Calling Sequence
```bash
python compute_tides_ICESat_GLA12.py --directory <path_to_directory> --tide <model> input_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat_GLA12.py)

#### Inputs
 1. `input_file`: input ICESat GLA12 file

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
 - `-I X`, `--interpolate X`: Interpolation method
     * `'spline'`
     * `'linear'`
     * `'nearest'`
     * `'bilinear'`
 - `-E X`, `--extrapolate X`: Extrapolate with nearest-neighbors
 - `-M X`, `--mode X`: Permission mode of output file
 - `-V`, `--verbose`: Output information about each created file
