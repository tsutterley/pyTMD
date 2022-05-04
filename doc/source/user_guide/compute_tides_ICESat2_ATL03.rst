==============================
compute_tides_ICESat2_ATL03.py
==============================

- Calculates tidal elevations for correcting ICESat-2 photon height data
- Calculated at ATL03 segment level using reference photon geolocation and time
- Segment level corrections can be applied to the individual photon events (PEs)
- Can use OTIS format tidal solutions provided by Ohio State University and ESR
- Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
- Can use Finite Element Solution (FES) models provided by AVISO

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_tides_ICESat2_ATL03.py

.. argparse::
    :filename: ../../scripts/compute_tides_ICESat2_ATL03.py
    :func: arguments
    :prog: compute_tides_ICESat2_ATL03.py
    :nodescription:
    :nodefault:

    --cutoff -c : @after
        * set to ``'inf'`` to extrapolate for all points
