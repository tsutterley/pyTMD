##################
`Release v2.1.1`__
##################

* ``refactor``: made the inferrence of minor constituents an option (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``refactor``: 1-liners in ``compute.py`` (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``refactor``: lunisolar ephemerides functions (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``feat``: added more constituent parameters for OTIS/ATLAS predictions (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``fix``: add option to return None and not raise error for Doodson numbers (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``docs``: add more definitions to the glossary (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``refactor``: moved constituent parameters function from ``predict`` to ``arguments`` (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``feat``: add functions for tide generating forces and potentials (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``fix``: variable typing for ``c`` in ``_constituent_parameters`` (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``test``: omit deprecated functions in coverage report (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``docs``: add ``toctree`` for ``io`` subdirectory (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``test``: add quick test for currents wrapper function (`#272 <https://github.com/tsutterley/pyTMD/pull/272>`_)
* ``fix``: construct OTIS currents masks differently if not global (`#273 <https://github.com/tsutterley/pyTMD/pull/273>`_)
* ``refactor``: renamed OTIS currents masks and bathymetry interpolation functions (`#273 <https://github.com/tsutterley/pyTMD/pull/273>`_)
* ``refactor``: renamed extend array and matrix functions (`#273 <https://github.com/tsutterley/pyTMD/pull/273>`_)
* ``docs``: add notebook showing tidal harmonic solver (`#275 <https://github.com/tsutterley/pyTMD/pull/275>`_)
* ``fix``: implicit import of ellipsoid constants class (`#275 <https://github.com/tsutterley/pyTMD/pull/275>`_)
* ``feat``: added inverse function to get currents from tide ellipse parameters (`#276 <https://github.com/tsutterley/pyTMD/pull/276>`_)
* ``refactor``: use complex algebra to calculate tidal ellipse parameters (`#276 <https://github.com/tsutterley/pyTMD/pull/276>`_)
* ``docs``: use ``importlib`` to prevent deprecation errors (`#276 <https://github.com/tsutterley/pyTMD/pull/276>`_)
* ``fix``: spelling mistake for solve notebook
* ``refactor``: changed class name for ellipsoid parameters to datum (`#287 <https://github.com/tsutterley/pyTMD/pull/287>`_)
* ``refactor``: move solve constants to subdirectory (`#287 <https://github.com/tsutterley/pyTMD/pull/287>`_)
* ``refactor``: move the immutable parameters in timescale class (`#287 <https://github.com/tsutterley/pyTMD/pull/287>`_)
* ``feat``: add capability to define a custom datum (`#287 <https://github.com/tsutterley/pyTMD/pull/287>`_)
* ``refactor``: changed variable for setting global grid flag to ``is_global`` (`#287 <https://github.com/tsutterley/pyTMD/pull/287>`_)
* ``fix``: doc strings for nodal arguments ``pu`` and ``pf`` (`#287 <https://github.com/tsutterley/pyTMD/pull/287>`_)
* ``refactor``: use numpy ``pad`` to interpolate data to u and v nodes (`#287 <https://github.com/tsutterley/pyTMD/pull/287>`_)
* ``feat``: can calculate polar stereographic distortion for distances (`#294 <https://github.com/tsutterley/pyTMD/pull/294>`_)
* ``docs``: update links to CATS2008-v2023 (`#294 <https://github.com/tsutterley/pyTMD/pull/294>`_)
* ``fix``: ``dtype`` suggestions (`#294 <https://github.com/tsutterley/pyTMD/pull/294>`_)
* ``fix``: append v currents in TPXO9 only if netcdf to address `#295 <https://github.com/tsutterley/pyTMD/issues/295>`_ (`#294 <https://github.com/tsutterley/pyTMD/pull/294>`_)

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.1.1
