compute_LPET_ICESat2_ATL06.py
=============================

- Calculates long-period equilibrium tidal elevations for correcting ICESat-2 land ice elevation data
- Will calculate the long-period tides for all ATL06 segments and not just ocean segments defined by the ocean tide mask

#### Calling Sequence
```bash
python compute_LPET_ICESat2_ATL06.py input_file
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPET_ICESat2_ATL06.py)

#### Inputs
 1. `input_file`: input ICESat-2 ATL06 file

#### Command Line Options
 - `-M X`, `--mode X`: Permission mode of output file
 - `-V`, `--verbose`: Output information about each created file
