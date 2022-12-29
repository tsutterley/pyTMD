##################
`Release v2.0.0`__
##################

- ``refactor``: single implicit import of pyTMD tools (`#130 <https://github.com/tsutterley/pyTMD/pull/130>`_)
- ``refactor``: reorganization of tide model readers under ``io`` (`#132 <https://github.com/tsutterley/pyTMD/pull/132>`_)
- ``refactor``: placed interpolation routines into new module (`#141 <https://github.com/tsutterley/pyTMD/pull/141>`_)
- ``refactor``: move ``model`` class to ``io``
- ``feat``: add constants class for ellipsoidal and gravitational parameters (`#135 <https://github.com/tsutterley/pyTMD/pull/135>`_)
- ``feat``: functions to read and interpolate from all constituents (`#137 <https://github.com/tsutterley/pyTMD/pull/137>`_)
- ``feat``: new functions to output ATLAS, FES and GOT netCDF4 files (`#139 <https://github.com/tsutterley/pyTMD/pull/139>`_)
- ``feat``: output variables for some standard epochs used within tide programs
- ``feat``: update forecast notebook with dynamic plotting
- ``fix``: copy input coordinates within read functions so they are not transformed
- ``docs``: update documentation for new structure
- ``docs``: standardized citation format throughout docstrings
- ``docs``: add release notes for all prior public releases
- ``docs``: add new pyTMD logo
- ``test``: read header from OPT test file and compare more variables
- ``test``: add tests for ``io`` methods

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.0.0
