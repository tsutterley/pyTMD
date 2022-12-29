#####################
`Release v1.0.2.18`__
#####################

- Simplified inputs to netcdf, GOT and FES tide read programs to be similar to binary OTIS program
- Added `TPXO9-atlas-v4` in `single constituent OTIS file format <https://www.tpxo.net/global/tpxo9-atlas>`_
- Updated documentation for consistency and to add contribution guidelines
- Add checks to extrapolation program to prevent runtime exceptions
- Updated ``spatial`` module adding polar stereographic area scale calculation
- Updated ``spatial`` module adding routines for cartesian coordinate conversion
- Prevent ``ComplexWarning`` for fill values when calculating amplitudes
- Prevent ``numpy`` DeprecationWarning for ``np.bool``, ``np.int`` and ``np.float``

.. __: https://github.com/tsutterley/pyTMD/releases/tag/1.0.2.18
