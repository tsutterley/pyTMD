#####################
`Release v1.0.2.11`__
#####################

- replaced ``griddata`` with `scipy regular grid interpolators <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html>`_
- allow small extrapolations with bilinear interpolation
- better incorporation of masked arrays in bilinear interpolation
- update date conversions in jupyter notebooks
- use ``pyproj2`` transformations in jupyter notebooks
- add ``time`` conversion test program
- add ocean load tide test program
- add comparison test for FES models
- add wrapper function test to perth3 workflow

.. __: https://github.com/tsutterley/pyTMD/releases/tag/1.0.2.11
