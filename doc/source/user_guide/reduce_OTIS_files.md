reduce_OTIS_files.py
====================

- Reads OTIS-format tidal files provided by Ohio State University and ESR and reduces to a regional subset

#### Calling Sequence
```bash
python reduce_OTIS_files.py --directory=<path_to_directory> --tide=<model> \
    --bounds=<xmin,xmax,ymin,ymax>
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/scripts/reduce_OTIS_files.py)

#### Command Line Options
 - `-D X`, `--directory=X`: Working data directory
 - `-B X`, `--bounds=X`: Grid Bounds (xmin,xmax,ymin,ymax)
 - `--date=X`: date to forecast in ISO format (YYYY-MM-DD)
 - `-T X`, `--tide=X`: Tide model to use
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
