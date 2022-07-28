==============================
compute_tides_ICESat2_ATL12.py
==============================

- Calculates tidal elevations for correcting ICESat-2 ocean surface height data
- Can use OTIS format tidal solutions provided by Ohio State University and ESR
- Can use Glgobal Tide Model (GOT) solutions provided by Richard Ray at GSFC
- Can use Finite Element Solution (FES) models provided by AVISO

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL12.py

.. argparse::
    :filename: ../../scripts/compute_tides_ICESat2_ATL12.py
    :func: arguments
    :prog: compute_tides_ICESat2_ATL12.py
    :nodescription:
    :nodefault:

    --cutoff -c : @after
        * set to ``'inf'`` to extrapolate for all points
