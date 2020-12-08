compute_LPET_ICESat_GLA12.py
============================

- Calculates long-period equilibrium tidal elevations for correcting ICESat/GLAS L2 GLA12 Antarctic and Greenland Ice Sheet elevation data
- Will calculate the long-period tides for all GLAS elevations and not just ocean elevations defined by the ocean tide mask

#### Calling Sequence
```bash
python compute_tides_ICESat_GLA12.py input_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPET_ICESat_GLA12.py)

#### Inputs
 1. `input_file`: input ICESat GLA12 file

#### Command Line Options
 - `-M X`, `--mode X`: Permission mode of output file
 - `-V`, `--verbose`: Output information about each created file
