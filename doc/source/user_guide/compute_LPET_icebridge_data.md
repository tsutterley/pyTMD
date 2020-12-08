compute_LPET_icebridge_data.py
=============================

 - Calculates long-period equilibrium tides for correcting Operation IceBridge elevation data
 - Uses the summation of fifteen tidal spectral lines from [Cartwright and Edden, (1973)](https://doi.org/10.1111/j.1365-246X.1973.tb03420.x)

#### Calling Sequence
```bash
python compute_LPET_icebridge_data.py --verbose input_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPET_icebridge_data.py)

#### Inputs
 1. `input_file`: input ATM1B, ATM icessn or LVIS file from NSIDC

#### Command Line Options
 - `-M X`, `--mode X`: Permission mode of output file
 - `-V`, `--verbose`: Output information about each created file
