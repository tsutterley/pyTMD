======
eop.py
======

Utilities for maintaining Earth Orientation Parameter (EOP) files

- Syncs mean pole files with IERS servers
- Can calculate update mean pole files using data from IERS servers
- Syncs finals orientation files with IERS servers

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/eop.py


General Methods
===============


.. method:: pyTMD.eop.update_mean_pole(verbose=False, mode=0o775)

    Connects to servers to download mean-pole.tab files from `IERS servers`__

.. __: ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab

    Keyword arguments:

        ``verbose``: print file information about output file

        ``mode``: permissions mode of output file


.. method:: pyTMD.eop.calculate_mean_pole(verbose=False, mode=0o775)

    Calculates the mean pole coordinates x and y are obtained by a Gaussian-weighted average of
    the `IERS pole coordinates <ftp://ftp.iers.org/products/eop/long-term/c01/eopc01.1900-now.dat>`_
    calculated following `IERS mean pole guidelines <ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.readme>`_

    Keyword arguments:

        ``verbose``: print file information about output file

        ``mode``: permissions mode of output file

    Returns:

        ``T``: date [decimal-years]

        ``xm``: mean pole coordinate x [arcsec]

        ``ym``: mean pole coordinate y [arcsec]


.. method:: pyTMD.eop.update_finals_file(username=None, password=None, verbose=False, mode=0o775)

    Connects to `servers`__ and downloads finals EOP files

.. __: ftp://cddis.nasa.gov/products/iers/

    Keyword arguments:

        ``username``: NASA Earthdata username

        ``password``: NASA Earthdata password

        ``verbose``: print file information about output file

        ``mode``: permissions mode of output file
