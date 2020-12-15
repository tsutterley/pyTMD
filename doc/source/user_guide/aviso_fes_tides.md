aviso_fes_tides.py
==================

 - Downloads the FES (Finite Element Solution) global tide model from [AVISO](https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html)
 - FES outputs are licensed for scientific purposes only
 - Decompresses the model tar files into the constituent files and auxiliary files
 - Must have [data access to tide models from AVISO](https://www.aviso.altimetry.fr/en/data/data-access.html)

#### Calling Sequence
```bash
python aviso_fes_tides.py --directory <path_to_tide_directory> --tide FES2014
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/aviso_fes_tides.py)

#### Command Line Options
 - `-D X`, `--directory X`: Working Data Directory
 - `--user X`: username for AVISO FTP servers
 - `-N X`, `--netrc X`: path to .netrc file for authentication
 - `-T X`, `--tide X`: FES tide model to download
    * `'FES1999'`
    * `'FES2004'`
    * `'FES2012'`
    * `'FES2014'`
 - `--load`: download load tide model outputs
 - `--currents`:  download tide model current outputs
 - `-M X`, `--mode X`: Permission mode of files downloaded
 - `-l`, `--log`: Output log of files downloaded
