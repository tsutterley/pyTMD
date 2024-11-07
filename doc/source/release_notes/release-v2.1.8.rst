##################
`Release v2.1.8`__
##################

* ``docs``: fix repository url fetch from ``Project-URL``
* ``docs``: update ``CITATION.cff`` to add version information
* ``docs``: add more definitions to glossary
* ``refactor``: convert Doodson coefficients table to JSON `(#353) <https://github.com/tsutterley/pyTMD/pull/353>`_
* ``feat``: added option to use Munk-Cartwright admittance interpolation for minor `(#353) <https://github.com/tsutterley/pyTMD/pull/353>`_
* ``feat``: add `Cartwright and Edden (1973) <http://dx.doi.org/10.1111/j.1365-246X.1971.tb01803.x>`_ table 1 `(#353) <https://github.com/tsutterley/pyTMD/pull/353>`_
* ``feat``: add `Cartwright and Tayler (1971) <http://dx.doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_ table 5 `(#353) <https://github.com/tsutterley/pyTMD/pull/353>`_
* ``feat``: add function to parse Cartwright/Tayler/Edden tables `(#353) <https://github.com/tsutterley/pyTMD/pull/353>`_
* ``feat``: add functions to calculate UKHO Extended Doodson numbers for constituents `(#353) <https://github.com/tsutterley/pyTMD/pull/353>`_
* ``test``: add test for extended doodson `(#353) <https://github.com/tsutterley/pyTMD/pull/353>`_
* ``docs``: add citations to included data `(#353) <https://github.com/tsutterley/pyTMD/pull/353>`_
* ``fix``: remove default bounds being ``None`` for `#356 <https://github.com/tsutterley/pyTMD/issues/356>`_ `(#357) <https://github.com/tsutterley/pyTMD/pull/357>`_
* ``docs``: move notebooks to docs and use myst to render `(#359) <https://github.com/tsutterley/pyTMD/pull/359>`_
* ``fix``: correct error when using default bounds in `extract_constants` for `#356 <https://github.com/tsutterley/pyTMD/issues/356>`_ `(#359) <https://github.com/tsutterley/pyTMD/pull/359>`_
* ``fix``: correct ``TPXO10-atlas-v2`` binary grid filename for `#358 <https://github.com/tsutterley/pyTMD/issues/358>`_ `(#359) <https://github.com/tsutterley/pyTMD/pull/359>`_
* ``fix``: some `Cartwright and Edden (1973) <http://dx.doi.org/10.1111/j.1365-246X.1971.tb01803.x>`_ table entries `(#359) <https://github.com/tsutterley/pyTMD/pull/359>`_
* ``docs``: use cards for notebook examples page `(#360) <https://github.com/tsutterley/pyTMD/pull/360>`_
* ``fix``: GOT5.6 names in database and add ``'n2'`` `(#360) <https://github.com/tsutterley/pyTMD/pull/360>`_
* ``docs``: set card grid to be either 1, 2 or 4
* ``fix``: update ``pyproject.toml`` for ``doc`` build
* ``fix``: allow variable case for Doodson number formalisms `(#361) <https://github.com/tsutterley/pyTMD/pull/361>`_
* ``feat``: added property for Extended Doodson numbers `(#361) <https://github.com/tsutterley/pyTMD/pull/361>`_
* ``fix``: use Love numbers for long-period tides when inferring (won't affect tilt factors) `(#361) <https://github.com/tsutterley/pyTMD/pull/361>`_
* ``docs``: add form factor notebook for classifying regional tides `(#361) <https://github.com/tsutterley/pyTMD/pull/361>`_

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.1.8
