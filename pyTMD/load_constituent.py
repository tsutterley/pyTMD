#!/usr/bin/env python
u"""
load_constituent.py (01/2024)
Loads parameters for a given tidal constituent

CALLING SEQUENCE:
    amplitude,phase,omega,alpha,species = load_constituent(c)

INPUTS:
    c: tidal constituent ID

OUTPUT:
    amplitude: amplitude of equilibrium tide in m for tidal constituent
    phase: phase of tidal constituent
    omega: angular frequency of constituent in radians
    alpha: load love number of tidal constituent
    species: spherical harmonic dependence of quadrupole potential

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 01/2024: moved constituent parameters function to arguments
    Updated 09/2023: deprecated in favor of pyTMD.predict function
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 07/2020: add more constituents from OTPSnc and function docstrings
    Updated 09/2017: Rewritten in Python
"""
import warnings
import pyTMD.arguments

def load_constituent(c):
    """
    Loads parameters for a given tidal constituent

    Parameters
    ----------
    c: list
        tidal constituent ID

    Returns
    -------
    amplitude: float
        amplitude of equilibrium tide in m for tidal constituent
    phase: float
        phase of tidal constituent
    omega: float
        angular frequency of constituent in radians
    alpha: float
        load love number of tidal constituent
    species: float
        spherical harmonic dependence of quadrupole potential

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    # raise warning for deprecated function call
    warnings.warn("Deprecated. Please use pyTMD.arguments instead",
        DeprecationWarning)
    # call updated function to not break current workflows
    return pyTMD.arguments._constituent_parameters(c)
