=======
compute
=======

- Calculates tidal elevations at points and times

  * Can use OTIS format tidal solutions provided by Oregon State University and ESR
  * Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
  * Can use Finite Element Solution (FES) models provided by AVISO
- Calculates tidal currents at points and times

  * Can use OTIS format tidal solutions provided by Oregon State University and ESR
  * Can use Finite Element Solution (FES) models provided by AVISO
- Calculates long-period equilibrium tides (LPET) at points and times

  * Uses the summation of fifteen tidal spectral lines from `Cartwright and Edden, (1973) <https://doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_
- Calculates radial pole load tides (LPT) at points and times

  * Following `IERS Convention (2010) guidelines <https://iers-conventions.obspm.fr/chapter7.php>`_
- Calculates radial ocean pole load tides (OPT) at points and times

  * Following `IERS Convention (2010) guidelines <https://iers-conventions.obspm.fr/chapter7.php>`_
- Calculates radial solid earth tides (SET) at points and times

  * Following `IERS Convention (2010) guidelines <https://iers-conventions.obspm.fr/chapter7.php>`_

Calling Sequence
----------------

.. code-block:: python

    import pyTMD
    tide_h = pyTMD.compute.tide_elevations(x, y, delta_time,
        DIRECTORY=DIRECTORY, MODEL=MODEL, EPOCH=(2000,1,1,0,0,0),
        EPSG=3031, TYPE='drift')
    tide_uv = pyTMD.compute.tide_currents(x, y, delta_time,
        DIRECTORY=DIRECTORY, MODEL=MODEL, EPOCH=(2000,1,1,0,0,0),
        EPSG=3031, TYPE='drift')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute.py

.. autofunction:: pyTMD.compute.corrections

.. autofunction:: pyTMD.compute.tide_elevations

.. autofunction:: pyTMD.compute.tide_currents

.. autofunction:: pyTMD.compute.LPET_elevations

.. autofunction:: pyTMD.compute.LPT_displacements

.. autofunction:: pyTMD.compute.OPT_displacements

.. autofunction:: pyTMD.compute.SET_displacements
