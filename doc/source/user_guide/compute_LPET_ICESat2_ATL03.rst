=============================
compute_LPET_ICESat2_ATL03.py
=============================

- Calculates long-period equilibrium tidal elevations for correcting ICESat-2 geolocated photon height data
- Calculated at ATL03 segment level using reference photon geolocation and time
- Segment level corrections can be applied to the individual photon events (PEs)
- Will calculate the long-period tides for all ATL03 segments and not just ocean segments defined by the ocean tide mask

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/compute_LPET_ICESat2_ATL03.py

.. argparse::
    :filename: ../../scripts/compute_LPET_ICESat2_ATL03.py
    :func: arguments
    :prog: compute_LPET_ICESat2_ATL03.py
    :nodescription:
    :nodefault:
