##################
`Release v2.1.5`__
##################

* ``feat``: adding GOT5.5 and GOT download program (`#316 <https://github.com/tsutterley/pyTMD/pull/316>`_)
* ``refactor``: renamed format for ``ATLAS`` to ``ATLAS-compact`` (`#316 <https://github.com/tsutterley/pyTMD/pull/316>`_)
* ``refactor``: renamed format for ``netcdf`` to ``ATLAS-netcdf`` (`#316 <https://github.com/tsutterley/pyTMD/pull/316>`_)
* ``refactor``: renamed format for ``FES`` to ``FES-netcdf`` and added ``FES-ascii`` (`#316 <https://github.com/tsutterley/pyTMD/pull/316>`_)
* ``refactor``: renamed format for ``GOT`` to ``GOT-ascii`` and added ``GOT-netcdf`` (`#316 <https://github.com/tsutterley/pyTMD/pull/316>`_)
* ``feat``: add JSON definition files for GOT5.5D and GOT5.6 (`#316 <https://github.com/tsutterley/pyTMD/pull/316>`_)
* ``feat``: add support for constituents in PERTH5 tables (`#316 <https://github.com/tsutterley/pyTMD/pull/316>`_)
* ``ci``: use upstream matlab TMD for OTIS comparison (`#316 <https://github.com/tsutterley/pyTMD/pull/316>`_)
* ``feat``: include inference of ``eps2`` and ``eta2`` with GOT models (`#319 <https://github.com/tsutterley/pyTMD/pull/319>`_)
* ``feat``: add attribute for minor constituents to model object (`#319 <https://github.com/tsutterley/pyTMD/pull/319>`_)
* ``feat``: allow inferring only specific minor constituents (`#319 <https://github.com/tsutterley/pyTMD/pull/319>`_)
* ``feat``: allow searching over iterable ``glob`` strings in definition files for `#318 <https://github.com/tsutterley/pyTMD/issues/318>`_ (`#319 <https://github.com/tsutterley/pyTMD/pull/319>`_)
* ``feat``: add option to auto-detect definition file format for `#318 <https://github.com/tsutterley/pyTMD/issues/318>`_ (`#319 <https://github.com/tsutterley/pyTMD/pull/319>`_)
* ``feat``: add back nodal arguments from PERTH3 for backwards compatibility (`#319 <https://github.com/tsutterley/pyTMD/pull/319>`_)
* ``chore``: trim trim trailing whitespace (`#319 <https://github.com/tsutterley/pyTMD/pull/319>`_)
* ``docs``: add GOT5.5 to getting started (`#319 <https://github.com/tsutterley/pyTMD/pull/319>`_)
* ``refactor``: change ``'geotiff'`` to ``'GTiff'`` and ``'cog'`` for `#320 <https://github.com/tsutterley/pyTMD/issues/320>`_ (`#321 <https://github.com/tsutterley/pyTMD/pull/321>`_)
* ``chore``: create ``pyproject.toml`` (`#321 <https://github.com/tsutterley/pyTMD/pull/321>`_)
* ``refactor``: modernize build with ``pyproject.toml`` (`#322 <https://github.com/tsutterley/pyTMD/pull/322>`_)
* ``feat``: add functions to calculate pole tides in cartesian coordinates for `#323 <https://github.com/tsutterley/pyTMD/issues/323>`_ (`#324 <https://github.com/tsutterley/pyTMD/pull/324>`_)
* ``refactor``: renamed io for Desai ocean pole tide file to ``IERS`` (`#324 <https://github.com/tsutterley/pyTMD/pull/324>`_)
* ``docs``: update prediction functions (`#324 <https://github.com/tsutterley/pyTMD/pull/324>`_)
* ``fix``: don't overwrite ocean pole tide longitude in shift (`#325 <https://github.com/tsutterley/pyTMD/pull/325>`_)
* ``test``: add more ocean pole tide verifications (`#325 <https://github.com/tsutterley/pyTMD/pull/325>`_)
* ``feat``: add ECEF to ENU conversions (`#326 <https://github.com/tsutterley/pyTMD/pull/326>`_)
* ``refactor``: use rotation matrix to convert from cartesian to spherical (`#326 <https://github.com/tsutterley/pyTMD/pull/326>`_)

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.1.5
