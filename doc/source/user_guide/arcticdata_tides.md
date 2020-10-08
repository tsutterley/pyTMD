arcticdata_tides.py
===================

 - Download Arctic Ocean Tide Models from the [NSF ArcticData](https://arcticdata.io) archive
 - [AODTM-5](https://arcticdata.io/catalog/view/doi:10.18739/A2901ZG3N)
 - [AOTIM-5](https://arcticdata.io/catalog/view/doi:10.18739/A2S17SS80)
 - [AOTIM-5-2018](https://arcticdata.io/catalog/view/doi:10.18739/A21R6N14K)

#### Calling Sequence
```bash
python arcticdata_tides.py --directory <path_to_tide_directory> --tide AOTIM-5-2018
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/arcticdata_tides.py)

#### Command Line Options
 - `-D X`, `--directory X`: Working Data Directory
 - `-T X`, `--tide X`: Arctic Ocean tide model to download
    * `'AODTM-5'`
    * `'AOTIM-5'`
    * `'AOTIM-5-2018'`
 - `-M X`, `--mode X`: Permission mode of files downloaded
