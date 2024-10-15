##################
`Release v2.1.7`__
##################

* ``docs``: improve description of optional dependencies in examples
* ``fix``: use case insensitive assertions of string argument values (`#340 <https://github.com/tsutterley/pyTMD/pull/340>`_)
* ``feat``: added bounded options for least squares solvers (`#341 <https://github.com/tsutterley/pyTMD/pull/341>`_)
* ``feat``: add ``__models__`` with all model names in database (`#341 <https://github.com/tsutterley/pyTMD/pull/341>`_)
* ``feat``: add function lists as ``__all__`` (`#341 <https://github.com/tsutterley/pyTMD/pull/341>`_)
* ``feat``: add `Ray and Erofeeva (2014) <https://doi.org/10.1002/2013JB010830>`_ to the database for `#327 <https://github.com/tsutterley/pyTMD/issues/327>`_ (`#341 <https://github.com/tsutterley/pyTMD/pull/341>`_)
* ``feat``: add minor inference for long period tides to address `#327 <https://github.com/tsutterley/pyTMD/issues/327>`_ (`#342 <https://github.com/tsutterley/pyTMD/pull/342>`_)
* ``fix``: try inferring both long and short period tides for FES (`#343 <https://github.com/tsutterley/pyTMD/pull/343>`_)
* ``refactor``: using new JSON dictionary format for model projections for `#333 <https://github.com/tsutterley/pyTMD/issues/333>`_ (`#345 <https://github.com/tsutterley/pyTMD/pull/345>`_)
* ``docs``: add notebook with a cotidal chart for `#344 <https://github.com/tsutterley/pyTMD/issues/344>`_ and `#348 <https://github.com/tsutterley/pyTMD/discussions/348>`_ (`#345 <https://github.com/tsutterley/pyTMD/pull/345>`_)
* ``docs``: update descriptions of coordinate reference systems
* ``feat``: add new functions to infer semi-diurnal and diurnal tides (`#346 <https://github.com/tsutterley/pyTMD/pull/346>`_)
* ``feat``: use PREM as the default Earth model for Love numbers (`#347 <https://github.com/tsutterley/pyTMD/pull/347>`_)
* ``feat``: compute delta times based on corrections type (`#347 <https://github.com/tsutterley/pyTMD/pull/347>`_)
* ``feat``: add wrapper functions to read and interpolate constants (`#349 <https://github.com/tsutterley/pyTMD/pull/349>`_)
* ``tests``: use simplified wrapper functions (`#349 <https://github.com/tsutterley/pyTMD/pull/349>`_)
* ``feat``: updated computation of long-period equilibrium tides (`#349 <https://github.com/tsutterley/pyTMD/pull/349>`_)
* ``feat``: add functions to append node tide equilibrium values to amplitudes (`#349 <https://github.com/tsutterley/pyTMD/pull/349>`_)
* ``fix``: add messaging if there are no minor constituents to infer (`#349 <https://github.com/tsutterley/pyTMD/pull/349>`_)
* ``feat``: can convert Doodson numbers formatted as strings (`#349 <https://github.com/tsutterley/pyTMD/pull/349>`_)
* ``docs``: expand the glossary (`#349 <https://github.com/tsutterley/pyTMD/pull/349>`_)

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.1.7
