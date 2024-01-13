##################
`Release v2.1.0`__
##################

* ``fix``: revert TPXO9-atlas currents changes to separate dicts for `#258 <https://github.com/tsutterley/pyTMD/issues/258>`_ (`#259 <https://github.com/tsutterley/pyTMD/pull/259>`_)
* ``test``: fix ``u`` and ``v`` for TPXO9-atlas netCDF (`#259 <https://github.com/tsutterley/pyTMD/pull/259>`_)
* ``fix``: escape sequences in docstrings to raw (`#259 <https://github.com/tsutterley/pyTMD/pull/259>`_)
* ``fix``: updated ssl context to fix deprecation error (`#259 <https://github.com/tsutterley/pyTMD/pull/259>`_)
* ``docs``: update docstrings for ssl context (`#259 <https://github.com/tsutterley/pyTMD/pull/259>`_)
* ``refactor``: use doodson arguments tables to calculate values (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``refactor``: rename ``phase_angles`` to ``doodson_arguments`` (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``refactor``: coordinate reference system class ``crs.py`` (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``refactor``: pass through ``VBox`` and ``HBox`` in ``tools.py`` (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``docs``: add link to TMD3 in ``Resources.rst`` (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: made keyword argument for selecting M1 coefficients (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: add initial solver for harmonic constants (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``chore``: include whitespace after commas (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: add function to create arguments coefficients table (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``test``: add check that arguments match prior version (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``docs``: started creating a glossary (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``docs``: add citation to Simon et al. (1994) (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: create arguments coefficients table for minor constituents (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: add function to calculate Doodson numbers for `#263 <https://github.com/tsutterley/pyTMD/discussions/263>`_ (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``refactor``: use mean lunar time as independent variable (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``refactor``: moved minor arguments calculation into new function (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``refactor``: implicit import of arguments (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``test``: add check for constants solve (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: add functions to convert to and from Doodson numbers (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: add option to output Cartwright numbers (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: add properties for Doodson and Cartwright numbers to ``constituents`` class (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)
* ``feat``: try to get the constituents of FES files (`#270 <https://github.com/tsutterley/pyTMD/pull/270>`_)

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.1.0
