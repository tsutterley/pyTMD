compute_LPT_icebridge_data.py
=============================

 - Calculates radial pole load tide displacements for correcting Operation IceBridge elevation data following IERS Convention (2010) guidelines
 - [http://maia.usno.navy.mil/conventions/2010officialinfo.php](http://maia.usno.navy.mil/conventions/2010officialinfo.php)
 - [http://maia.usno.navy.mil/conventions/chapter7.php](http://maia.usno.navy.mil/conventions/chapter7.php)

#### Calling Sequence
```bash
python compute_LPT_icebridge_data.py input_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPT_icebridge_data.py)

#### Inputs
 1. `input_file`: input ATM1B, ATM icessn or LVIS file from NSIDC

#### Command Line Options
 - `-M X`, `--mode X`: Permission mode of output file
 - `-V`, `--verbose`: Output information about each created file
