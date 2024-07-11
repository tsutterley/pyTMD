##################
`Release v2.1.2`__
##################

* ``refactor``: use ``timescale`` for EOP and temporal operations (`#300 <https://github.com/tsutterley/pyTMD/pull/300>`_)
* ``docs``: fix ``netCDF4`` urls (`#300 <https://github.com/tsutterley/pyTMD/pull/300>`_)
* ``refactor``: remove older deprecated functions (`#300 <https://github.com/tsutterley/pyTMD/pull/300>`_)
* ``feat``: add debug mode printing input arguments and additional information (`#300 <https://github.com/tsutterley/pyTMD/pull/300>`_)
* ``feat``: wrapper to ``importlib`` for optional dependencies (`#300 <https://github.com/tsutterley/pyTMD/pull/300>`_)
* ``feat``: add io for (geo)parquet datasets with geometry columns (`#303 <https://github.com/tsutterley/pyTMD/pull/303>`_)
* ``feat``: make classes subscriptable and allow item assignment (`#303 <https://github.com/tsutterley/pyTMD/pull/303>`_)
* ``fix``: default module import as ``class`` (`#303 <https://github.com/tsutterley/pyTMD/pull/303>`_)
* ``feat``: add output ``to_file`` function (`#303 <https://github.com/tsutterley/pyTMD/pull/303>`_)
* ``fix``: deprecation update to replace ``np.longcomplex`` (`#303 <https://github.com/tsutterley/pyTMD/pull/303>`_)
* ``fix``: ``numpy`` 2.0 fix for time calculation (`#303 <https://github.com/tsutterley/pyTMD/pull/303>`_)
* ``test``: add parquet io tests (`#303 <https://github.com/tsutterley/pyTMD/pull/303>`_)
* ``fix``: prevent integer overflows with ``numpy`` 2.0 (`#304 <https://github.com/tsutterley/pyTMD/pull/304>`_)
* ``docs``: add references to new PERTH5 software from Richard Ray (`#305 <https://github.com/tsutterley/pyTMD/pull/305>`_)
* ``feat``: add wrapper function for normalizing angles (`#305 <https://github.com/tsutterley/pyTMD/pull/305>`_)
* ``feat``: add functions to convert to and from Degree-Minutes-Seconds (DMS) (`#305 <https://github.com/tsutterley/pyTMD/pull/305>`_)
* ``fix``: assert that data type is a known value (`#305 <https://github.com/tsutterley/pyTMD/pull/305>`_)
* ``feat``: added new ``FES2022`` and ``FES2022_load`` to list of models (`#307 <https://github.com/tsutterley/pyTMD/pull/307>`_)
* ``feat``: only download FES files if non-existent or updated (`#308 <https://github.com/tsutterley/pyTMD/pull/308>`_)
* ``refactor``: add ``_jd_j2000`` variable instead of hard coded (`#308 <https://github.com/tsutterley/pyTMD/pull/308>`_)

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.1.2
