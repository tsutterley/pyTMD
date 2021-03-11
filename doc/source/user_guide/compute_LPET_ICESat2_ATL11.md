compute_LPET_ICESat2_ATL11.py
=============================

- Calculates long-period equilibrium tidal elevations for correcting ICESat-2 annual land ice height data
- Calculates long-period equilibrium tidal elevations for both along-track and across-track locations
- Will calculate the long-period tides for all ATL11 segments and crossovers and not just ocean segments defined by the ocean tide mask

#### Calling Sequence
```bash
python compute_LPET_ICESat2_ATL11.py input_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPET_ICESat2_ATL11.py)

#### Inputs
 1. `input_file`: input ICESat-2 ATL11 file

#### Command Line Options
 - `-M X`, `--mode X`: Permission mode of output file
 - `-V`, `--verbose`: Output information about each created file
