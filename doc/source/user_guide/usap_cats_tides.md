usap_cats_tides.py
==================

 - Download Circum-Antarctic Tidal Simulations from the [US Antarctic Program](https://www.usap-dc.org)
 - [CATS2008](https://www.usap-dc.org/view/dataset/601235)

#### Calling Sequence
```bash
python usap_cats_tides.py --directory <path_to_tide_directory> --tide CATS2008
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/usap_cats_tides.py)

#### Command Line Options
 - `-D X`, `--directory X`: Working Data Directory
 - `--tide X`: Circum-Antarctic Tidal Simulation to download
    * `'CATS2008'`
 - `-M X`, `--mode X`: Permission mode of files downloaded
