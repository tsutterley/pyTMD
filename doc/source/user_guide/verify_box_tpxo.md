verify_box_tpxo.py
==================

 - Verifies TPXO9-atlas global tide models downloaded from the box file sharing service
 - Compares `sha1` hashes to verify the binary or netCDF4 files

#### Calling Sequence
```bash
python verify_box_tpxo.py --directory <path_to_tide_directory> \
   --tide TPXO9-atlas-v4 --token <box_api_token> --folder <box_folder_id>
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/verify_box_tpxo.py)

#### Command Line Options
 - `-D X`, `--directory X`: Working data directory
 - `-t X`, `--token X`: User access token for box API
 - `-F X`, `--folder X`: box folder id for model
 - `-T X`, `--tide X`: TPXO9-atlas model to verify
    * `'TPXO9-atlas'`
    * `'TPXO9-atlas-v2'`
    * `'TPXO9-atlas-v3'`
    * `'TPXO9-atlas-v4'`
 - `--currents`:  verify tide model current outputs
 - `-M X`, `--mode X`: Permission mode of files downloaded
