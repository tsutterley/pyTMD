##################
`Release v2.1.6`__
##################

* ``refactor``: use JSON database for known model parameters for `#328 <https://github.com/tsutterley/pyTMD/issues/328>`_ (`#329 <https://github.com/tsutterley/pyTMD/pull/329>`_)
* ``feat``: added new ``TPXO10-atlas-v2`` to list of models (`#329 <https://github.com/tsutterley/pyTMD/pull/329>`_)
* ``feat``: generalize hash function to use any available algorithm (`#329 <https://github.com/tsutterley/pyTMD/pull/329>`_)
* ``fix``: use model name in default output filename for definition file case (`#329 <https://github.com/tsutterley/pyTMD/pull/329>`_)
* ``refactor``: create database from providers (`#329 <https://github.com/tsutterley/pyTMD/pull/329>`_)
* ``refactor``: drop support for the ascii definition file format (`#329 <https://github.com/tsutterley/pyTMD/pull/329>`_)
* ``docs``: add model providers section to contributions (`#329 <https://github.com/tsutterley/pyTMD/pull/329>`_)
* ``ci``: add GitHub Action to update ``database.json`` (`#330 <https://github.com/tsutterley/pyTMD/pull/330>`_)
* ``feat``: Add additional models to provider JSON (`#331 <https://github.com/tsutterley/pyTMD/pull/331>`_)
* ``fix``: add parse constituents back to model load (`#331 <https://github.com/tsutterley/pyTMD/pull/331>`_)
* ``fix``: drop constituents from database (`#331 <https://github.com/tsutterley/pyTMD/pull/331>`_)
* ``ci``: only run ``pytest`` action if secrets are accessible (`#331 <https://github.com/tsutterley/pyTMD/pull/331>`_)
* ``refactor``: use idealized Azimuthal equidistant for Arctic models (`#332 <https://github.com/tsutterley/pyTMD/pull/332>`_)
* ``fix``: deprecation in case where an array is output to scalars (`#332 <https://github.com/tsutterley/pyTMD/pull/332>`_)
* ``fix``: j1 and theta for FES type models for `#335 <https://github.com/tsutterley/pyTMD/issues/335>`_ (`#336 <https://github.com/tsutterley/pyTMD/pull/336>`_)
* ``fix``: nodal corrections for ``eps2`` and ``eta2`` when inferring for GOT5.5 (`#336 <https://github.com/tsutterley/pyTMD/pull/336>`_)
* ``feat``: use model class attributes for file format and corrections (`#336 <https://github.com/tsutterley/pyTMD/pull/336>`_)
* ``feat``: add option to select nodal corrections type (`#336 <https://github.com/tsutterley/pyTMD/pull/336>`_)

.. __: https://github.com/tsutterley/pyTMD/releases/tag/2.1.6
