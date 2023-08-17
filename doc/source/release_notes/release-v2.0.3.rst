##################
`Release v2.0.3`__
##################

* ``feat``: add basic variable typing to function inputs (`#174 <https://github.com/tsutterley/pyTMD/pull/174>`_)
* ``feat``: added "one-liners" for long-period equilibrium tide (LPET), load pole tide (LPT) and ocean pole tide (OPT) (`#176 <https://github.com/tsutterley/pyTMD/pull/176>`_)
* ``feat``: added option to change IERS mean or secular pole convention (`#176 <https://github.com/tsutterley/pyTMD/pull/176>`_)
* ``feat``: added 2018 IERS secular pole convention (`#176 <https://github.com/tsutterley/pyTMD/pull/176>`_)
* ``feat``: set ellipsoid name and output units as ``constants`` attributes (`#176 <https://github.com/tsutterley/pyTMD/pull/176>`_)
* ``feat``: add `HAMTIDE11` model to address `#179 <https://github.com/tsutterley/pyTMD/issues/179>`_ (`#180 <https://github.com/tsutterley/pyTMD/pull/180>`_)
* ``feat``: adding work for computing solid earth tides (`#186 <https://github.com/tsutterley/pyTMD/pull/186>`_)
* ``feat``: add solid Earth tide (SET) correction program for files (`#186 <https://github.com/tsutterley/pyTMD/pull/186>`_)
* ``feat``: add function for phase angles (`#186 <https://github.com/tsutterley/pyTMD/pull/186>`_)
* ``test``: add solid Earth tide (SET) checks vs IERS and ICESat-2 (`#186 <https://github.com/tsutterley/pyTMD/pull/186>`_)
* ``refactor``: renamed coordinate reference system conversion functions (`#174 <https://github.com/tsutterley/pyTMD/pull/174>`_)
* ``refactor``: mapping notebooks for matplotlib 3.5 (`#182 <https://github.com/tsutterley/pyTMD/pull/182>`_)
* ``fix``: setting directories for ``FES`` currents within ``model`` class (`#182 <https://github.com/tsutterley/pyTMD/pull/182>`_)
* ``fix``: check if ``datetime`` before converting to seconds (`#186 <https://github.com/tsutterley/pyTMD/pull/186>`_)
* ``fix``: copy inputs in cartesian to not modify original arrays (`#186 <https://github.com/tsutterley/pyTMD/pull/186>`_)
* ``docs``: remove deprecated ``.rst`` files (`#174 <https://github.com/tsutterley/pyTMD/pull/174>`_)
* ``docs``: update documentation to denote new solid Earth tide (SET) functionality (`#186 <https://github.com/tsutterley/pyTMD/pull/186>`_)

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.0.3
