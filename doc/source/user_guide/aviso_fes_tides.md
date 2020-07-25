aviso_fes_tides.py
==================

 - Downloads the FES (Finite Element Solution) global tide model from [AVISO](https://www.aviso.altimetry.fr/data/products/auxiliary-products/global-tide-fes.html)
 - Decompresses the model tar files into the constituent files and auxiliary files
 - Must have [data access to tide models from AVISO](https://www.aviso.altimetry.fr/en/data/data-access.html)

#### Calling Sequence
```bash
python aviso_fes_tides.py --directory=<path_to_tide_directory> --tide=fes2014
```
[Source code](https://github.com/tsutterley/pyTMD/blob/master/scripts/aviso_fes_tides.py)

#### Command Line Options
 - `-D X`, `--directory=X`: Working Data Directory
 - `--user=X`: username for AVISO FTP servers
 - `-N X, --netrc=X`: path to .netrc file for authentication
 - `--tide=X`: FES tide model to download
    * fes1999
    * fes2004
    * fes2012
    * fes2014
 - `--load`: download load tide model outputs
 - `--currents`:  download tide model current outputs
 - `-M X`, `--mode=X`: Permission mode of files downloaded
 - `-l`, `--log`: Output log of files downloaded
