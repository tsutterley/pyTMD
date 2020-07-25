#!/usr/bin/env python
u"""
load_constituent.py (07/2020)
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

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 07/2020: add more constituents from OTPSnc and function docstrings
    Updated 09/2017: Rewritten in Python
"""
import numpy as np

def load_constituent(c):
    """
    Loads parameters for a given tidal constituent

    Arguments
    ---------
    c: tidal constituent ID

    Returns
    -------
    amplitude: amplitude of equilibrium tide in m for tidal constituent
    phase: phase of tidal constituent
    omega: angular frequency of constituent in radians
    alpha: load love number of tidal constituent
    species: spherical harmonic dependence of quadrupole potential
    """
    #-- constituents array that are included in tidal program
    cindex = ['m2','s2','k1','o1','n2','p1','k2','q1','2n2','mu2','nu2','l2',
        't2','j1','m1','oo1','rho1','mf','mm','ssa','m4','ms4','mn4','m6','m8',
        'mk3','s6','2sm2','2mk3']
    #-- species type (spherical harmonic dependence of quadrupole potential)
    species_all = np.array([2,2,1,1,2,1,2,1,2,2,2,2,2,1,1,1,1,0,0,0,0,0,0,
        0,0,0,0,0,0])
    #-- loading love number
    #-- alpha = correction factor for first order load tides
    alpha_all = np.array([0.693,0.693,0.736,0.695,0.693,0.706,0.693,0.695,0.693,
        0.693,0.693,0.693,0.693,0.695,0.695,0.695,0.695,0.693,0.693,0.693,
        0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693])
    #-- omega: angular frequency of constituent, in radians
    omega_all = np.array([1.405189e-04,1.454441e-04,7.292117e-05,6.759774e-05,
        1.378797e-04,7.252295e-05,1.458423e-04,6.495854e-05,1.352405e-04,
        1.355937e-04,1.382329e-04,1.431581e-04,1.452450e-04,7.556036e-05,
        7.028195e-05,7.824458e-05,6.531174e-05,0.053234e-04,0.026392e-04,
        0.003982e-04,2.810377e-04,2.859630e-04,2.783984e-04,4.215566e-04,
        5.620755e-04,2.134402e-04,4.363323e-04,1.503693e-04,2.081166e-04])
    #-- Astronomical arguments (relative to t0 = 1 Jan 0:00 1992)
    #-- phases for each constituent are referred to the time when the phase of
    #-- the forcing for that constituent is zero on the Greenich meridian
    phase_all = np.array([1.731557546,0.000000000,0.173003674,1.558553872,
        6.050721243,6.110181633,3.487600001,5.877717569,4.086699633,
        3.463115091,5.427136701,0.553986502,0.052841931,2.137025284,
        2.436575100,1.929046130,5.254133027,1.756042456,1.964021610,
        3.487600001,3.463115091,1.731557546,1.499093481,5.194672637,
        6.926230184,1.904561220,0.000000000,4.551627762,3.809122439])
    #-- amplitudes of equilibrium tide in m
    # amplitude_all = np.array([0.242334,0.112743,0.141565,0.100661,0.046397,
    amplitude_all = np.array([0.2441,0.112743,0.141565,0.100661,0.046397,
        0.046848,0.030684,0.019273,0.006141,0.007408,0.008811,0.006931,0.006608,
        0.007915,0.007915,0.004338,0.003661,0.042041,0.022191,0.019567,0.,0.,0.,
        0.,0.,0.,0.,0.,0.])

    #-- map between input constituent and cindex
    j = [j for j,val in enumerate(cindex) if (val == c.lower())]
    #-- set the values for the constituent
    if j:
        amplitude, = amplitude_all[j]
        phase, = phase_all[j]
        omega, = omega_all[j]
        alpha, = alpha_all[j]
        species, = species_all[j]
    else:
        amplitude = 0.0; phase = 0.0; omega = 0.0; alpha = 0.0; species = 0
    #-- return the values for the constituent
    return (amplitude,phase,omega,alpha,species)
