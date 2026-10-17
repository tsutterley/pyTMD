#!/usr/bin/env python
"""
math.py
Written by Tyler Sutterley (10/2026)
Special functions of mathematical physics

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

UPDATE HISTORY:
    Updated 10/2026: added option to calculate exact factorial values
        set to False for cases with large arguments that would cause overflow
    Updated 09/2026: added recurrence function for Legendre polynomials
        added option to calculate fully-normalized Legendre polynomials
    Updated 06/2026: standardize use of lambda (lmda) to denote longitudes
    Updated 05/2026: added kronecker delta function and updated docstrings
    Updated 03/2026: add radius and scalar product functions
        calculate Legendre polynomials using Hofmann-Wellenhof (2006) eq. 1.67
        split out Condon-Shortley phase and Legendre normalization functions
    Updated 02/2026: add inverse functions for converting radians to arcseconds
    Updated 11/2025: legendre now returns both polynomials and derivatives
    Updated 09/2025: added degree 4 to legendre polynomials function
    Updated 08/2025: add asec2rad and masec2rad functions for arcseconds
    Updated 07/2025: add deriv and phase arguments to sph_harm function
        add Legendre polynomial derivatives with respect to theta
    Updated 04/2025: use numpy power function over using pow for consistency
    Updated 01/2025: added function for fully-normalized Legendre polynomials
    Updated 12/2024: added function to calculate an aliasing frequency
    Written 11/2024
"""

from __future__ import annotations

import numpy as np
from scipy.special import factorial, gammaln

__all__ = [
    "asec2rad",
    "masec2rad",
    "rad2asec",
    "rad2masec",
    "polynomial_sum",
    "normalize_angle",
    "radius",
    "rotate",
    "scalar_product",
    "aliasing",
    "legendre",
    "_assoc_legendre",
    "_condon_shortley",
    "_kronecker_delta",
    "_legendre_norm",
    "legendreP",
    "sph_harm",
]


def asec2rad(x: float | np.ndarray):
    """
    Convert angles from arcseconds to radians

    Parameters
    ----------
    x: float or np.ndarray
        Input angle (arcseconds)
    """
    return np.radians(x / 3600.0)


def masec2rad(x: float | np.ndarray):
    """
    Convert angles from microarcseconds to radians

    Parameters
    ----------
    x: float or np.ndarray
        Input angle (microarcseconds)
    """
    return np.radians(x / 3.6e9)


def rad2asec(x: float | np.ndarray):
    """
    Convert angles from radians to arcseconds

    Parameters
    ----------
    x: float or np.ndarray
        Input angle (radians)
    """
    return 3600.0 * np.degrees(x)


def rad2masec(x: float | np.ndarray):
    """
    Convert angles from radians to microarcseconds

    Parameters
    ----------
    x: float or np.ndarray
        Input angle (radians)
    """
    return 3.6e9 * np.degrees(x)


# PURPOSE: calculate the sum of a polynomial function of time
def polynomial_sum(
    coefficients: list | np.ndarray,
    t: np.ndarray,
):
    """
    Calculates the sum of a polynomial function using Horner's method
    :cite:p:`Horner:1819br`

    Parameters
    ----------
    coefficients: list or np.ndarray
        Leading coefficient of polynomials of increasing order
    t: np.ndarray
        Time for a given astronomical longitudes calculation
    """
    # convert time to array if importing a single value
    t = np.atleast_1d(t)
    return np.sum([c * (t**i) for i, c in enumerate(coefficients)], axis=0)


def normalize_angle(
    theta: float | np.ndarray,
    circle: float = 360.0,
):
    """
    Normalize an angle to a single rotation

    Parameters
    ----------
    theta: float or np.ndarray
        Angle to normalize
    circle: float, default 360.0
        Circle of the angle
    """
    return np.mod(theta, circle)


def radius(
    x: float | np.ndarray,
    y: float | np.ndarray,
    z: float | np.ndarray,
):
    """
    Calculate the radius from the origin to a point in 3-dimensional space

    Parameters
    ----------
    x: float or np.ndarray
        x-coordinate of the point
    y: float or np.ndarray
        y-coordinate of the point
    z: float or np.ndarray
        z-coordinate of the point
    """
    return np.sqrt(x**2 + y**2 + z**2)


def rotate(
    theta: float | np.ndarray,
    axis: str = "x",
):
    """
    Rotate a 3-dimensional matrix about a given axis

    Parameters
    ----------
    theta: float or np.ndarray
        Angle of rotation (radians)
    axis: str, default 'x'
        Axis of rotation (``'x'``, ``'y'``, or ``'z'``)
    """
    # allocate for output rotation matrix
    R = np.zeros((3, 3, len(np.atleast_1d(theta))))
    if axis.lower() == "x":
        # rotate about x-axis
        R[0, 0, :] = 1.0
        R[1, 1, :] = np.cos(theta)
        R[1, 2, :] = np.sin(theta)
        R[2, 1, :] = -np.sin(theta)
        R[2, 2, :] = np.cos(theta)
    elif axis.lower() == "y":
        # rotate about y-axis
        R[0, 0, :] = np.cos(theta)
        R[0, 2, :] = -np.sin(theta)
        R[1, 1, :] = 1.0
        R[2, 0, :] = np.sin(theta)
        R[2, 2, :] = np.cos(theta)
    elif axis.lower() == "z":
        # rotate about z-axis
        R[0, 0, :] = np.cos(theta)
        R[0, 1, :] = np.sin(theta)
        R[1, 0, :] = -np.sin(theta)
        R[1, 1, :] = np.cos(theta)
        R[2, 2, :] = 1.0
    else:
        raise ValueError(f"Invalid axis {axis}")
    # return the rotation matrix
    return R


def scalar_product(
    x: float | np.ndarray,
    y: float | np.ndarray,
    z: float | np.ndarray,
    u: float | np.ndarray,
    v: float | np.ndarray,
    w: float | np.ndarray,
):
    """
    Calculate the scalar product of two vectors in 3-dimensional space

    Parameters
    ----------
    x: float or np.ndarray
        x-coordinate of the first vector
    y: float or np.ndarray
        y-coordinate of the first vector
    z: float or np.ndarray
        z-coordinate of the first vector
    u: float or np.ndarray
        x-coordinate of the second vector
    v: float or np.ndarray
        y-coordinate of the second vector
    w: float or np.ndarray
        z-coordinate of the second vector
    """
    return x * u + y * v + z * w


def aliasing(
    f: float,
    fs: float,
) -> float:
    """
    Calculate the aliasing frequency of a signal

    Parameters
    ----------
    f: float
        Frequency of the signal
    fs: float
        Sampling frequency of the signal

    Returns
    -------
    fa: float
        Aliasing frequency of the signal
    """
    fa = np.abs(f - fs * np.round(f / fs))
    return fa


def legendre(
    l: int,
    x: np.ndarray,
    m: int = 0,
    norm: float = 1.0,
    exact: bool = True,
):
    r"""
    Computes associated Legendre functions and their first-derivatives
    for a particular degree and order
    :cite:p:`Munk:1966go,HofmannWellenhof:2006hy`

    Parameters
    ----------
    l: int
        Degree of the Legendre polynomials
    x: np.ndarray
        Elements ranging from -1 to 1

        Typically :math:`\cos(\theta)`, where :math:`\theta`
        is the colatitude
    m: int, default 0
        Order of the Legendre polynomials (:math:`0` to :math:`l`)
    norm: float, default 1.0
        Normalization to apply to outputs
    exact: bool, default True
        Calculate factorials using integer arithmetic

    Returns
    -------
    Plm: np.ndarray
        Legendre polynomials of degree :math:`l` and order :math:`m`
    dPlm: np.ndarray
        First derivative of Legendre polynomials with respect to
        :math:`\theta`
    """
    # verify values are integers
    l = np.int64(l)
    m = np.int64(m)
    # assert values for spherical harmonic orders
    if (m < 0) or (m > l):
        raise ValueError("Order must be between 0 and l")
    # verify x is array
    if isinstance(x, list):
        x = np.atleast_1d(x)
    # equation 1.67 from Hofmann-Wellenhof (2006)
    Plm = _assoc_legendre(l, m, x, exact=exact)
    # if x is the cos of colatitude, u is the sine
    u = np.sqrt(1.0 - x**2)
    # calculate first derivative
    # this will initially have a singularity at the poles
    Pm1 = _assoc_legendre(l - 1, m, x, exact=exact)
    # ignore divide by zero and invalid value warnings
    with np.errstate(divide="ignore", invalid="ignore"):
        dPlm = (l * x * Plm - (l + m) * Pm1) / u
    # handle singularity at the poles (x = +/-1) for order 1 terms
    pole = -0.5 * np.power(x, l) * l * (l + 1) if m == 1 else 0.0
    dPlm = np.where(np.isclose(u, 0.0), pole, dPlm)
    # return the associated legendre functions
    # and their first derivatives with respect to theta
    return norm * Plm, norm * dPlm


def _assoc_legendre(
    l: int,
    m: int,
    x: np.ndarray,
    exact: bool = True,
):
    r"""
    Computes associated Legendre polynomials following equation 1.67
    from :cite:t:`HofmannWellenhof:2006hy`

    Parameters
    ----------
    l: int
        Degree of the Legendre polynomials
    m: int
        Order of the Legendre polynomials (:math:`0` to :math:`l`)
    x: np.ndarray
        Elements ranging from -1 to 1

        Typically :math:`\cos(\theta)`, where :math:`\theta`
        is the colatitude
    exact: bool, default True
        Calculate factorials using integer arithmetic

    Returns
    -------
    Plm: np.ndarray
        Legendre polynomials of degree :math:`l` and order :math:`m`
    """
    # verify values are integers
    l = np.int64(l)
    m = np.int64(m)
    # return 0 for invalid values
    if m > l:
        return 0.0
    # verify x is array
    if isinstance(x, list):
        x = np.atleast_1d(x)
    # if x is the cos of colatitude, u is the sine
    u = np.sqrt(1.0 - x**2)
    # calculate un-normalized polynomials
    # equation 1.67 from Hofmann-Wellenhof (2006)
    P = 0.0
    r = np.floor_divide(l - m, 2)
    for k in range(r + 1):
        # calculate factorial terms
        if exact:
            # using integer arithmetic for exact calculations
            fact = (
                factorial(2 * l - 2 * k)
                / factorial(k)
                / factorial(l - k)
                / factorial(l - m - 2 * k)
            )
        else:
            # using gamma functions for cases with potential overflow
            fact = np.exp(
                gammaln(2.0 * l - 2.0 * k + 1.0)
                - gammaln(k + 1.0)
                - gammaln(l - k + 1.0)
                - gammaln(l - m - 2.0 * k + 1.0)
            )
        # add terms
        P += np.power(-1.0, k) * fact * np.power(x, l - m - 2.0 * k)
    # calculate for degree l and order m
    Plm = P * np.power(2.0, -l) * np.power(u, m)
    # apply Condon-Shortley phase
    Plm *= _condon_shortley(m)
    # return the associated legendre polynomials
    return Plm


def _condon_shortley(m: int | np.ndarray):
    r"""
    Computes the Condon-Shortley phase :math:`(-1)^m` for order :math:`m`

    Parameters
    ----------
    m: int or np.ndarray
        Order of the Legendre polynomials
    """
    return np.power(-1.0, m)


def _kronecker_delta(
    i: int | np.ndarray,
    j: int | np.ndarray,
):
    r"""
    Computes the Kronecker delta :math:`\delta_{ij}` function

    .. math::
        \delta_{ij} =
            \begin{cases}
                1 & \text{if } i = j \\
                0 & \text{if } i \neq j
            \end{cases}

    Parameters
    ----------
    i: int or np.ndarray
        First index
    j: int or np.ndarray
        Second index
    """
    return 1.0 * (i == j)


def _legendre_norm(
    l: int,
    m: int,
    normalization: str = "schmidt",
    exact: bool = True,
):
    r"""
    Calculates Legendre polynomials normalizations following
    :cite:t:`Munk:1966go`

    - Schmidt semi-normalized :cite:p:`Schmidt:1917vh`

        .. math::
            N_l^m = \sqrt{\frac{(l - m)!}{(l + m)!}}
    - Fully-normalized :cite:p:`HofmannWellenhof:2006hy`

        .. math::
            N_l^m = \sqrt{(2 - \delta_{m0}) (2 l + 1) \frac{(l - m)!}{(l + m)!}}

    Parameters
    ----------
    l: int
        Degree of the Legendre polynomials
    m: int
        Order of the Legendre polynomials (:math:`0` to :math:`l`)
    normalization: str, default schmidt
        Normalization scheme

        - ``'schmidt'`` from :cite:t:`Schmidt:1917vh`
        - ``'full'`` from :cite:t:`HofmannWellenhof:2006hy`
    exact: bool, default True
        Calculate factorials using integer arithmetic

    Returns
    -------
    norm: float
        Normalization for degree :math:`l` and order :math:`m`
    """
    # calculate factorial terms for the normalization
    if exact:
        # using integer arithmetic for exact calculations
        F = factorial(l - m) / factorial(l + m)
    else:
        # using gamma functions for cases with potential overflow
        F = np.exp(gammaln(l - m + 1) - gammaln(l + m + 1))
    # normalization from Munk and Cartwright (1966) equation A5
    if normalization.lower() == "schmidt":
        norm = np.sqrt(F)
    elif normalization.lower() == "full":
        # fully-normalized spherical harmonics (geodesy formulation)
        # from equation 1.97 of Hofmann-Wellenhof (2006)
        # also equation 6.2b of IERS conventions (2010)
        kron = _kronecker_delta(m, 0)
        norm = np.sqrt(2.0 - kron) * np.sqrt(2.0 * l + 1.0) * np.sqrt(F)
    else:
        raise ValueError(f"Unknown normalization {normalization}")
    # return the normalization
    return norm


def legendreP(lmax, x):
    r"""
    Computes fully-normalized associated Legendre Polynomials
    using the recursion relation from :cite:t:`Mohlenkamp:2016vv`

    Derived from the :cite:t:`Szego:1939tn` recurrence formula for
    Jacobi Polynomials

    Parameters
    ----------
    lmax: int
        Maximum degree of Legendre polynomials
    x: np.ndarray
        Elements ranging from -1 to 1

        Typically :math:`\cos(\theta)`, where :math:`\theta`
        is the colatitude

    Returns
    -------
    Plm: np.ndarray
        Fully-normalized Legendre polynomials
    dPlm: np.ndarray
        First derivative of Legendre polynomials with respect to
        :math:`\theta`
    """
    # verify values are integers
    lmax = np.int64(lmax)
    # verify length of the x array
    x = np.atleast_1d(x)
    n = len(x)
    # if x is the cos of colatitude, u is the sine
    u = np.sqrt(1.0 - x**2)
    # Jacobi polynomials for degree and order lmax
    Jlm = np.zeros((lmax + 1, lmax + 1, n))
    # allocate for output Legendre Polynomials
    # fully-normalized to geodesy convention
    Plm = np.zeros((lmax + 1, lmax + 1, n))
    dPlm = np.zeros((lmax + 1, lmax + 1, n))
    # iterate over spherical harmonic orders
    for m in range(0, lmax + 1):
        # Jacobi terms for degree 0
        Jlm[0, m, :] = 1.0 / np.sqrt(2.0)
        # apply factors for orders greater than 0
        for j in range(1, m + 1):
            Jlm[0, m, :] *= np.sqrt((2.0 * j + 1.0) / (2.0 * j))
        # Jk,m,m Terms
        for k in range(1, lmax + 1):
            # degree dependent factors
            f1 = np.sqrt((k + m - 0.5) / k)
            f2 = np.sqrt((k + m + 0.5) / (k + 2.0 * m))
            # calculate Jacobi general terms
            if k == 1:
                # for degree 1 terms
                Jlm[k, m, :] = 2.0 * x * f1 * f2 * Jlm[k - 1, m, :]
            else:
                # for all other spherical harmonic degrees
                f3 = np.sqrt(
                    (2.0 * k + 2.0 * m + 1.0) / (2.0 * k + 2.0 * m - 3.0)
                )
                f4 = np.sqrt((k + 2.0 * m - 1.0) / (k + 2.0 * m))
                Jm1 = 2.0 * x * f1 * f2 * Jlm[k - 1, m, :]
                Jm2 = np.sqrt((k - 1.0) / k) * f3 * f4 * Jlm[k - 2, m, :]
                Jlm[k, m, :] = Jm1 - Jm2
        # calculate and normalize Legendre polynomials
        for l in range(m, lmax + 1):
            # handle cases for normalization
            if (l == 0) & (m == 0):
                # degree and order 0
                Plm[l, m, :] = np.sqrt(2.0) * Jlm[l, m, :]
                dPlm[l, m, :] = 0.0 * x
                continue
            elif m == 0:
                # positive degree and order 0
                Plm[l, m, :] = np.sqrt(2.0) * Jlm[l, m, :]
                # calculate normalization for order 0
                # ignore divide by zero and invalid value warnings
                with np.errstate(divide="ignore", invalid="ignore"):
                    norm = np.sqrt(2.0) / u
                # handle singularity at the poles
                norm = np.where(np.isclose(u, 0.0), 0.0, norm)
            else:
                # all other degrees and orders
                Plm[l, m, :] = 2.0 * np.power(u, m) * Jlm[l - m, m, :]
                # calculate normalization for order m
                norm = 2.0 * np.power(u, m - 1)
            # calculate first derivatives w.r.t. colatitude
            f5 = np.sqrt(l**2.0 - m**2.0)
            f6 = np.sqrt((2.0 * l + 1.0) / (2.0 * l - 1.0))
            dPlm[l, m, :] = norm * (
                l * x * Jlm[l - m, m, :] - f5 * f6 * Jlm[l - m - 1, m, :]
            )
    # return the associated Legendre polynomials
    return Plm, dPlm


def sph_harm(
    l: int,
    theta: np.ndarray,
    lmda: np.ndarray,
    m: int = 0,
    phase: float = 0.0,
    exact: bool = True,
):
    r"""
    Computes the spherical harmonics for a particular degree
    and order :cite:p:`Munk:1966go,HofmannWellenhof:2006hy`

    .. math::
        Y_l^m = \sqrt{\frac{(2l + 1)(l - m)!}{4\pi(l + m)!}}
        P_l^m(\cos\theta)e^{i(m\lambda + \varphi)}

    Parameters
    ----------
    l: int
        Degree of the spherical harmonics
    theta: np.ndarray
        Colatitude :math:`\theta` (radians)
    lmda: np.ndarray
        Longitude :math:`\lambda` (radians)
    m: int, default 0
        Order of the spherical harmonics (:math:`0` to :math:`l`)
    phase: float, default 0.0
        Phase shift :math:`\varphi` (radians)
    exact: bool, default True
        Calculate factorials using integer arithmetic

    Returns
    -------
    Ylm: np.ndarray
        Complex spherical harmonics of degree :math:`l` and order :math:`m`
    dYlm: np.ndarray
        First derivative of spherical harmonics with respect to
        :math:`\theta`
    """
    # normalization from Munk and Cartwright (1966) equation A5
    norm = _legendre_norm(l, m, exact=exact)
    # calculate associated Legendre functions and derivatives
    Plm, dPlm = legendre(l, np.cos(theta), m=m, norm=norm, exact=exact)
    # normalized spherical harmonics of degree l and order m
    dfactor = np.sqrt((2.0 * l + 1.0) / (4.0 * np.pi))
    Ylm = dfactor * Plm * np.exp(1j * m * lmda + 1j * phase)
    dYlm = dfactor * dPlm * np.exp(1j * m * lmda + 1j * phase)
    # return values
    return Ylm, dYlm
