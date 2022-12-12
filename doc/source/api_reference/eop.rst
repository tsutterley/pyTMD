===
eop
===

Utilities for maintaining Earth Orientation Parameter (EOP) files

- Syncs mean pole files with IERS servers
- Can calculate update mean pole files using data from IERS servers
- Syncs finals orientation files with IERS servers

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/eop.py


General Methods
===============

.. autofunction:: pyTMD.eop.update_mean_pole

.. autofunction:: pyTMD.eop.calculate_mean_pole

.. autofunction:: pyTMD.eop.pull_pole_coordinates

.. autofunction:: pyTMD.eop.update_finals_file

.. autofunction:: pyTMD.eop.iers_mean_pole

.. autofunction:: pyTMD.eop.iers_daily_EOP
