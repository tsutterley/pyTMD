.. _spherical-harmonics:

===================
Spherical Harmonics
===================

The tide potential at colatitude :math:`\theta` and longitude :math:`\lambda` is often expressed in terms of spherical harmonic functions :math:`Y_l^m(\theta,\lambda)` of degree :math:`l` and order :math:`m` :cite:p:`Doodson:1921kt`.
These harmonics are solutions of Laplace's equation in spherical coordinates using separation of variables :cite:p:`HofmannWellenhof:2006hy`.
The degree 2 spherical harmonic terms are the dominant source of tidal excitation, and induce the semi-diurnal (from :math:`Y_2^2`), diurnal (from :math:`Y_2^1`) and long-period (from :math:`Y_2^0`) tides :cite:p:`Ray:2020gn`.

.. plot:: ./background/sphharm.py
    :caption: Spherical harmonics for degree 2
    :align: center


Mathematical Definition
=======================

Spherical harmonics of degree :math:`l` and order :math:`m` are defined as:

.. math::
    :label: 9.1
    :name: eq:9.1

    Y_l^m(\theta, \lambda) = N_l^m\, P_l^m(\cos\theta)\, e^{im\lambda}

where :math:`P_l^m(x)` are the associated Legendre polynomials and :math:`N_l^m` is a normalization factor :cite:p:`HofmannWellenhof:2006hy,Munk:1960uk`.
:py:func:`pyTMD.math.sph_harm` calculates spherical harmonics using the normalization factor from Equation A5 in :cite:t:`Munk:1966go` following :cite:t:`Schmidt:1917vh`.
The associated Legendre polynomials of degree :math:`l` and order :math:`m` are calculated in :py:func:`pyTMD.math.legendre` using the explicit formula (1-67) from :cite:t:`HofmannWellenhof:2006hy`:

.. math::
    :label: 9.2
    :name: eq:9.2

    P_l^m(x) = 2^{-l}(1-x^2)^{m/2}\sum_{k=0}^{n}(-1)^k\frac{(2l-2k)!}{k!(l-k)!(l-m-2k)!}\, x^{l-m-2k}

where :math:`n` is the largest integer less than or equal to :math:`(l-m)/2`.
The first few (unnormalized) Legendre polynomials are:

.. list-table::
    :header-rows: 0
    :align: center
    :class: flat

    * - :math:`P_0^0(x) = 1`
      -
      -
    * - :math:`P_1^0(x) = x`
      - :math:`P_1^1(x) = \sqrt{1-x^2}`
      -
    * - :math:`P_2^0(x) = \tfrac{1}{2}(3x^2 - 1)`
      - :math:`P_2^1(x) = 3x\sqrt{1-x^2}`
      - :math:`P_2^2(x) = 3(1-x^2)`

The infinite sums of the Legendre polynomials as given in :cite:t:`Hobson:1955ur` and :cite:t:`Farrell:1972cm` are:

.. math::
    :label: 9.3
    :name: eq:9.3

    \sum_{l=0}^\infty P_l(\cos{\theta}) &= \frac{1}{\sin{(\theta/2)}}\\
    \sum_{l=0}^\infty l\,P_l(\cos{\theta}) &= -\frac{1}{4\sin{(\theta/2)}}\\
    \sum_{l=0}^\infty \frac{\partial P_l(\cos{\theta})}{\partial\theta} &= -\frac{\cos{(\theta/2)}}{4\sin^2{(\theta/2)}}

These infinite sums can be used to calculate the :ref:`greens-functions` of deformation using a truncated set of Legendre polynomials.

.. _tide-generating-potential:

Tide-Generating Potential
=========================

The gravitational potential :math:`W` at a location on Earth's surface :math:`(\varphi, \lambda)` from a planetary body (such as the Moon or Sun) can be expanded into spherical harmonics as the following:

.. math::
    :label: 9.4
    :name: eq:9.4

    W(\varphi, \lambda) = \frac{G M}{R} \sum_{l=2}^{\infty} \left(\frac{r}{R}\right)^l P_l(\cos\psi)

where :math:`G` is the gravitational constant, :math:`M` is the mass of the body, :math:`r` is the radius of the Earth, :math:`R` is the distance to the planetary body, and :math:`\psi` is the :term:`zenith angle <Zenith Angle>` between the planetary body and the position on the Earth's surface :cite:p:`Munk:1960uk,Ray:2020gn`.
The *standard gravitational parameter* (:math:`GM`) of the planetary body can be derived using that of the Earth (:math:`GM_E`) and their ratio of masses (:math:`M/M_E`).
The cosine of the :ref:`zenith angle <zenith-angle>` (:math:`\cos\psi`) can also be calculated using the dot product between the geocentric unit vectors of the :ref:`tide-generating body <ephemerides>` (:math:`\hat{\mathbf{R}}`) and the point on the Earth's surface (:math:`\hat{\mathbf{r}}`) :cite:p:`Mathews:1997js`:

.. math::
    :label: 9.5
    :name: eq:9.5

    \cos\psi = \hat{\mathbf{r}}\cdot\hat{\mathbf{R}}
             = \frac{\mathbf{r}\cdot\mathbf{R}}{|\mathbf{r}\mathbf{R}|}

Higher-Degree Terms
===================

Global asymmetry in the tide potential can lead to a dependence on higher degree harmonics, most notably the degree 3 and 4 terms.
These terms are small compared to those of degree 2, but have been detected at both local and global scales :cite:p:`Ray:2020gn`.
Both GNSS stations and superconducting gravimeters can be sensitive enough to detect the signals from these higher degree terms :cite:p:`Hartmann:1995jp,Petit:2010tp`.
The formalism for computing solid Earth tides within the IERS Conventions include the component of deformation induced by the degree 3 terms :cite:p:`Petit:2010tp,Mathews:1997js`.
Catalogs of tide potential, such as ``HW1995`` [see :ref:`tab-catalogs`], can include even higher degree terms, as well as the potentials related to planetary motion :cite:p:`Hartmann:1995jp`.

.. plot:: ./background/sphharm34.py
    :caption: Spherical harmonics for degrees 3 and 4
    :align: center


.. list-table:: Maximum Tide Potential from :cite:t:`Hartmann:1995jp`
    :header-rows: 1
    :align: center

    * - Body
      - Degree
      - V [m\ :sup:`2`/s\ :sup:`2`]
    * - Moon
      - 2
      - 4.41
    * - Moon
      - 3
      - 7.88\ |times|\ 10\ :sup:`-2`
    * - Moon
      - 4
      - 1.41\ |times|\ 10\ :sup:`-3`
    * - Moon
      - 5
      - 2.53\ |times|\ 10\ :sup:`-5`
    * - Moon
      - 6
      - 4.52\ |times|\ 10\ :sup:`-7`
    * - Moon
      - 7
      - 8.06\ |times|\ 10\ :sup:`-9`
    * - Sun
      - 2
      - 1.60
    * - Sun
      - 3
      - 6.80\ |times|\ 10\ :sup:`-5`
    * - Sun
      - 4
      - 2.89\ |times|\ 10\ :sup:`-9`

.. |times|      unicode:: U+00D7 .. MULTIPLICATION SIGN
