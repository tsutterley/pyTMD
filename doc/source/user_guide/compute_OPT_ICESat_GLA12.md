compute_OPT_ICESat_GLA12.py
===========================

- Calculates radial ocean pole tide displacements for correcting ICESat/GLAS L2 GLA12 Antarctic and Greenland Ice Sheet elevation data following IERS Convention (2010) guidelines
 - [http://maia.usno.navy.mil/conventions/2010officialinfo.php](http://maia.usno.navy.mil/conventions/2010officialinfo.php)
 - [http://maia.usno.navy.mil/conventions/chapter7.php](http://maia.usno.navy.mil/conventions/chapter7.php)
 - [ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz](ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz)

#### Calling Sequence
```bash
python compute_OPT_ICESat_GLA12.py input_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_OPT_ICESat_GLA12.py)

#### Inputs
 1. `input_file`: input ICESat GLA12 file

#### Command Line Options
 - `-I X`, `--interpolate X`: Interpolation method
     * `'spline'`
     * `'linear'`
     * `'nearest'`
 - `-M X`, `--mode X`: Permission mode of output file
 - `-V`, `--verbose`: Output information about each created file
