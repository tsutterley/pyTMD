compute_tides_ICESat2_ATL06.py
==============================

- Calculates tidal elevations for correcting ICESat-2 land ice elevation data
- Can use OTIS format tidal solutions provided by Ohio State University and ESR
- Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

#### Calling Sequence
```bash
python compute_tides_ICESat2_ATL06.py --directory=<path_to_directory> --tide=<model> input_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/compute_tides_ICESat2_ATL06.py)

#### Inputs
 1. `input_file`: input ICESat-2 ATL06 file

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
