========================
compute_tide_corrections
========================

- Calculates tidal elevations for correcting elevation or imagery data

  * Can use OTIS format tidal solutions provided by Ohio State University and ESR
  * Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
  * Can use Finite Element Solution (FES) models provided by AVISO
- Calculates long-period equilibrium tides (LPET) for correcting elevation or imagery data

  * Uses the summation of fifteen tidal spectral lines from `Cartwright and Edden, (1973) <https://doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_
- Calculates radial pole load tides (LPT) for correcting elevation or imagery data following `IERS Convention (2010) guidelines <https://iers-conventions.obspm.fr/chapter7.php>`_
- Calculates radial ocean pole load tides (OPT) for correcting elevation or imagery data following `IERS Convention (2010) guidelines <https://iers-conventions.obspm.fr/chapter7.php>`_
- Calculates radial solid earth tides (SET) for correcting elevation or imagery data following `IERS Convention (2010) guidelines <https://iers-conventions.obspm.fr/chapter7.php>`_

Calling Sequence
----------------

.. code-block:: python

    from pyTMD.compute_tide_corrections import compute_tide_corrections
    tide = compute_tide_corrections(x, y, delta_time, DIRECTORY=DIRECTORY,
        MODEL=MODEL, EPOCH=(2000,1,1,0,0,0), EPSG=3031, TYPE='drift')

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/compute_tide_corrections.py

.. autofunction:: pyTMD.compute_corrections

.. autofunction:: pyTMD.compute_tide_corrections

.. autofunction:: pyTMD.compute_LPET_corrections

.. autofunction:: pyTMD.compute_LPT_corrections

.. autofunction:: pyTMD.compute_OPT_corrections

.. autofunction:: pyTMD.compute_SET_corrections
