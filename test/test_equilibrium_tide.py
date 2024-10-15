"""
test_solid_earth.py (10/2024)
Tests the calculation of long-period equilibrium tides with respect
to the LPEQMT subroutine

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    timescale: Python tools for time and astronomical calculations
        https://pypi.org/project/timescale/

UPDATE HISTORY:
    Written 10/2024
"""
import pytest
import numpy as np
import pyTMD.predict
import timescale.time

# PURPOSE: test the estimation of long-period equilibrium tides
@pytest.mark.parametrize("TYPE", ['grid','drift'])
def test_equilibrium_tide(TYPE):
    """
    Test the computation of the long-period equilibrium tides
    from the summation of fifteen tidal spectral lines from
    Cartwright-Tayler-Edden tables [1]_ [2]_

    References
    ----------
    .. [1] D. E. Cartwright and R. J. Tayler,
        "New Computations of the Tide-generating Potential,"
        *Geophysical Journal of the Royal Astronomical Society*,
        23(1), 45--73. (1971). `doi: 10.1111/j.1365-246X.1971.tb01803.x
        <https://doi.org/10.1111/j.1365-246X.1971.tb01803.x>`_
    .. [2] D. E. Cartwright and A. C. Edden,
        "Corrected Tables of Tidal Harmonics,"
        *Geophysical Journal of the Royal Astronomical Society*,
        33(3), 253--264, (1973). `doi: 10.1111/j.1365-246X.1973.tb03420.x
        <https://doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_
    """
    # create a test dataset for data type
    if (TYPE == 'drift'):
        # number of data points
        n_time = 3000
        lat = np.random.randint(-90,90,size=n_time)
        delta_time = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'grid'):
        # number of data points
        n_lat,n_time = (181,100)
        lat = np.linspace(-90,90,n_lat)
        delta_time = np.random.randint(0,31557600,size=n_time)

    # convert from seconds since 2018 to tide time
    EPOCH = (2018, 1, 1, 0, 0, 0)
    t = timescale.from_deltatime(delta_time, epoch=EPOCH, standard='GPS')
    # calculate long-period equilibrium tides
    lpet = pyTMD.predict.equilibrium_tide(t.tide, lat)

    # longitude of moon
    # longitude of sun
    # longitude of lunar perigee
    # longitude of ascending lunar node
    PHC = np.array([290.21,280.12,274.35,343.51])
    DPD = np.array([13.1763965,0.9856473,0.1114041,0.0529539])

    # number of input points
    nt = len(np.atleast_1d(t.tide))
    nlat = len(np.atleast_1d(lat))
    # compute 4 principal mean longitudes in radians at delta time (SHPN)
    SHPN = np.zeros((4,nt))
    for N in range(4):
        # convert time from days relative to 1992-01-01 to 1987-01-01
        ANGLE = PHC[N] + (t.tide + 1826.0)*DPD[N]
        SHPN[N,:] = np.pi*pyTMD.astro.normalize_angle(ANGLE)/180.0

    # assemble long-period tide potential from 15 CTE terms greater than 1 mm
    # nodal term is included but not the constant term.
    PH = np.zeros((nt))
    Z = np.zeros((nt))
    Z += 2.79*np.cos(SHPN[3,:]) - 0.49*np.cos(SHPN[1,:] - \
        283.0*np.pi/180.0) - 3.10*np.cos(2.0*SHPN[1,:])
    PH += SHPN[0,:]
    Z += -0.67*np.cos(PH - 2.0*SHPN[1,:] + SHPN[2,:]) - \
        (3.52 - 0.46*np.cos(SHPN[3,:]))*np.cos(PH - SHPN[2,:])
    PH += SHPN[0,:]
    Z += - 6.66*np.cos(PH) - 2.76*np.cos(PH + SHPN[3,:]) - \
        0.26 * np.cos(PH + 2.*SHPN[3,:]) - 0.58 * np.cos(PH - 2.*SHPN[1,:]) - \
        0.29 * np.cos(PH - 2.*SHPN[2,:])
    PH += SHPN[0,:]
    Z += - 1.27*np.cos(PH - SHPN[2,:]) - \
        0.53*np.cos(PH - SHPN[2,:] + SHPN[3,:]) - \
        0.24*np.cos(PH - 2.0*SHPN[1,:] + SHPN[2,:])

    # Multiply by gamma_2 * normalization * P20(lat)
    k2 = 0.302
    h2 = 0.609
    gamma_2 = (1.0 + k2 - h2)
    P20 = 0.5*(3.0*np.sin(lat*np.pi/180.0)**2 - 1.0)
    # calculate long-period equilibrium tide and convert to meters
    if (nlat != nt):
        exp = gamma_2*np.sqrt((4.0+1.0)/(4.0*np.pi))*np.outer(P20,Z/100.0)
    else:
        exp = gamma_2*np.sqrt((4.0+1.0)/(4.0*np.pi))*P20*(Z/100.0)
    # compare with functional values
    eps = np.finfo(np.float16).eps
    assert np.all(np.abs(lpet - exp) < eps)

# PURPOSE: test the estimation of long-period equilibrium tides
def test_node_tide():
    """
    Test the computation of the equilibrium node tides

    References
    ----------
    .. [1] D. E. Cartwright and R. J. Tayler,
        "New Computations of the Tide-generating Potential,"
        *Geophysical Journal of the Royal Astronomical Society*,
        23(1), 45--73. (1971). `doi: 10.1111/j.1365-246X.1971.tb01803.x
        <https://doi.org/10.1111/j.1365-246X.1971.tb01803.x>`_
    .. [2] D. E. Cartwright and A. C. Edden,
        "Corrected Tables of Tidal Harmonics,"
        *Geophysical Journal of the Royal Astronomical Society*,
        33(3), 253--264, (1973). `doi: 10.1111/j.1365-246X.1973.tb03420.x
        <https://doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_
    """
    # create a test dataset for data type
    # number of data points
    n_time = 3000
    lat = np.random.randint(-90,90,size=n_time)
    delta_time = np.random.randint(0,31557600,size=n_time)
    # convert from seconds since 2018 to tide time
    EPOCH = (2018, 1, 1, 0, 0, 0)
    t = timescale.from_deltatime(delta_time,
        epoch=EPOCH, standard='GPS')
    # calculate long-period equilibrium tides
    lpet = pyTMD.predict.equilibrium_tide(t.tide, lat,
        constituents='node', corrections='GOT')
    # calculate node tide amplitude
    amp, ph = pyTMD.io.model.node_equilibrium(lat)
    # calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    hc = amp*np.exp(cph)
    tide = pyTMD.predict.drift(t.tide, hc[:,None],
        constituents=['node'], corrections='GOT')
    # compare with functional values
    eps = np.finfo(np.float16).eps
    assert np.all(np.abs(lpet - tide) < eps)
