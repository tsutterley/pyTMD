#!/usr/bin/env python
u"""
test_love_numbers.py (11/2024)
Verify Love numbers correspond to those from Wahr et al. (1981)

UPDATE HISTORY:
    Updated 11/2024: moved love number calculator to arguments
    Written 09/2024
"""
import pytest
import numpy as np
import pyTMD.arguments

def test_love_numbers():
    """
    Tests the calculation of body tide Love numbers compared
    with the 1066A values from Wahr et al. (1981)
    """
    # expected values
    exp = {}
    # diurnal species
    exp['2q1'] = (0.604, 0.298, 0.0841)
    exp['sigma1'] = (0.604, 0.298, 0.0841)
    exp['q1'] = (0.603, 0.298, 0.0841)
    exp['rho1'] = (0.603, 0.298, 0.0841)
    exp['o1'] = (0.603, 0.298, 0.0841)
    exp['tau1'] = (0.603, 0.298, 0.0842)
    exp['m1'] = (0.600, 0.297, 0.0842)
    exp['chi1'] = (0.600, 0.296, 0.0843)
    exp['pi1'] = (0.587, 0.290, 0.0847)
    exp['p1'] = (0.581, 0.287, 0.0849)
    exp['s1'] = (0.568, 0.280, 0.0853)
    exp['k1'] = (0.520, 0.256, 0.0868)
    exp['psi1'] = (0.937, 0.466, 0.0736)
    exp['phi1'] = (0.662, 0.328, 0.0823)
    exp['theta1'] = (0.612, 0.302, 0.0839)
    exp['j1'] = (0.611, 0.302, 0.0839)
    exp['so1'] = (0.608, 0.301, 0.0840)
    exp['oo1'] = (0.608, 0.301, 0.0840)
    exp['ups1'] = (0.607, 0.300, 0.0840)
    # semi-diurnal species
    exp['m2'] = (0.609, 0.302, 0.0852)
    # long-period species
    exp['mm'] = (0.606, 0.299, 0.0840)
    # for each tidal constituent
    for c, v in exp.items():
        # calculate Love numbers
        omega, = pyTMD.arguments.frequency(c)
        h2, k2, l2 = pyTMD.arguments._love_numbers(
            omega, model='1066A')
        # check Love numbers
        assert np.isclose(h2, v[0], atol=15e-4)
        assert np.isclose(k2, v[1], atol=15e-4)
        assert np.isclose(l2, v[2], atol=15e-4)
