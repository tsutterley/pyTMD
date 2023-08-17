##################
`Release v2.0.1`__
##################

* ``feat``: default geotiff output as cog
* ``feat``: added default field mapping for reading from netCDF4/HDF5 (`#152 <https://github.com/tsutterley/pyTMD/pull/152>`_)
* ``feat``: split netCDF4 output for ``grid`` and ``drift`` types to address `#154 <https://github.com/tsutterley/pyTMD/issues/154>`_ (`#159 <https://github.com/tsutterley/pyTMD/pull/159>`_)
* ``feat``: use debug level logging instead of import warnings in ``tools.py`` to address `#156 <https://github.com/tsutterley/pyTMD/issues/156>`_ (`#159 <https://github.com/tsutterley/pyTMD/pull/159>`_)
* ``feat``: add ``time series`` type to address `#153 <https://github.com/tsutterley/pyTMD/discussions/153>`_ (`#162 <https://github.com/tsutterley/pyTMD/pull/162>`_)
* ``fix``: verify warnings have type and only show once `#146 <https://github.com/tsutterley/pyTMD/issues/146>`_ (`#147 <https://github.com/tsutterley/pyTMD/pull/147>`_)
* ``fix``: pin ``scipy`` to 1.9.3 for `scipy/scipy#17716 <https://github.com/scipy/scipy/issues/17716>`_ (`#147 <https://github.com/tsutterley/pyTMD/pull/147>`_)
* ``fix``: use default context from ``utilities`` module (`#147 <https://github.com/tsutterley/pyTMD/pull/147>`_)
* ``fix``: include more possible dimension names for ``grid`` and ``time series`` outputs (`#163 <https://github.com/tsutterley/pyTMD/pull/163>`_)
* ``test``: add default field mapping test (`#152 <https://github.com/tsutterley/pyTMD/pull/152>`_)
* ``test``: validate gridded and time series netCDF/HDF5 io (`#163 <https://github.com/tsutterley/pyTMD/pull/163>`_)
* ``docs``: update documentation colors to match new logo
* ``docs``: don't have metavar for ``--tide`` to address `#155 <https://github.com/tsutterley/pyTMD/issues/155>`_ (`#159 <https://github.com/tsutterley/pyTMD/pull/159>`_)
* ``docs``: updated v2 link for Arc2km to address `#157 <https://github.com/tsutterley/pyTMD/issues/157>`_ (`#159 <https://github.com/tsutterley/pyTMD/pull/159>`_)

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.0.1
