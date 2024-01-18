===============
ocean_pole_tide
===============

- Reads ocean pole load tide coefficients provided by IERS as computed by `Desai et al. (2002) <https://doi.org/10.1029/2001JC001224>`_ and `Desai et al. (2015) <https://doi.org/10.1007/s00190-015-0848-7>`_
- See `materials from Chapter 7 of the IERS Conventions <https://webtai.bipm.org/iers/convupdt/convupdt_c7.html>`_
- `Ocean Pole Tide File <ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/additional_info/opoleloadcoefcmcor.txt.gz>`_

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.io
    import pyTMD.utilities
    ocean_pole_tide_file = pyTMD.utilities.get_data_path(['data','opoleloadcoefcmcor.txt.gz'])
    ur,un,ue,glon,glat = pyTMD.io.ocean_pole_tide(ocean_pole_tide_file)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/io/ocean_pole_tide.py

.. autofunction:: pyTMD.io.ocean_pole_tide
