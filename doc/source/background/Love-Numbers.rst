.. _love-and-shida-numbers:

======================
Love and Shida Numbers
======================

When the Moon or Sun exerts a tidal force on the Earth, the solid body of the Earth deforms both elastically and inelastically [see :ref:`solid-earth-tides`].
This deformation has three distinct components:

1. a vertical (radial) displacement of the surface
2. a horizontal (tangential) displacement of the surface
3. a change in the gravitational potential

The magnitude of each individual component can be characterised by a dimensionless scaling factor :cite:p:`Love:1909eh,Shida:1912dj`.
These factors are collectively known as :term:`Love/Shida Numbers <Love and Shida Numbers>`, and are defined for each :ref:`spherical harmonic degree <spherical-harmonics>`.

.. list-table::
    :header-rows: 1
    :align: center

    * - Symbol
      - Name
    * - :math:`h_l`
      - Love number (vertical)
    * - :math:`k_l`
      - Love number (potential)
    * - :math:`l_l`
      - Shida number (horizontal)

.. _tidal-love-numbers:

Tidal Love Numbers
------------------

Tidal-effective Love numbers describe the deformation of the solid Earth in response to the *direct gravitational potential* of a remote body (such as the Sun or Moon).
The tidal forcing is distributed *throughout* the entirety of the Earth's volume :cite:p:`Munk:1960uk`.

Elastic Love numbers assume that the Earth responds instantaneously to the tidal forcing and without any dissipation.
For a spherical, non-rotating Earth, the Love/Shida numbers are largely independent of tidal frequency as the tidal periods are longer than the Earth's free oscillation periods :cite:p:`Baker:1984tq,Wahr:1979vx,Wahr:1981ea`.
However, for a rotating, ellipsoidal Earth, the Love/Shida numbers have some dependence on tidal frequency, with resonances particularly in the diurnal band :cite:p:`Wahr:1979vx,Wahr:1981ea,Lambeck:1980ic,Ray:2017jx`.

.. plot:: ./background/love-numbers.py
    :caption: Diurnal frequency dependence of :term:`Love/Shida numbers <Love and Shida Numbers>` from :cite:t:`Wahr:1979vx`
    :align: center

Additionally, the Earth's mantle is not perfectly elastic, and there is a small phase lag between the tidal forcing and the deformation response.
Complex Love numbers contain a *real part* describing the in-phase (elastic) response and an *imaginary part* describing the out-of-phase (dissipative) response :cite:p:`Wahr:1981ea,Petit:2010tp`.
``pyTMD`` computes both the frequency-dependent corrections and the dissipative responses following :cite:t:`Mathews:1997js` and :cite:t:`Wahr:1981ea`.

Combinations of Love/Shida numbers can derive additional quantities :cite:p:`Baker:1984tq,Farrell:1970tn,Cartwright:1999tj,Merriam:1973wi,Munk:1960uk`:

.. list-table:: :term:`Love/Shida numbers <Love and Shida Numbers>` combinations
    :header-rows: 1
    :align: center

    * - Name
      - Expression
    * - :term:`Gravimetric Factor`
      - :math:`\delta_l = 1 + \dfrac{2h_l}{l} - \dfrac{(l + 1)k_l}{l}`
    * - :term:`Tilt Factor`
      - :math:`\gamma_l = 1 + k_l - h_l`

``pyTMD`` uses these factors in a few different computations, including the calculation of :ref:`long-period equilibrium tides <equilibrium-theory>`, :ref:`ocean pole Tides <pole-tides>` and :ref:`gravity tides <gravity-tides>`. 
Similar to the original Love/Shida numbers, the tidal factors will have a dependence on tidal frequency :cite:p:`Melchior:1983wd,Wahr:1979vx` and contain imaginary components for the dissipative response :cite:p:`Wahr:1981ea`.

.. plot:: ./background/tidal-factors.py
    :caption: Diurnal frequency dependence of the :term:`gravimetric <Gravimetric Factor>` and :term:`tilt <Tilt Factor>` factors
    :align: center

.. _load-love-numbers:

Load Love Numbers
-----------------

Load Love numbers describe the deformation of the solid Earth in response to a change in *surface mass load* :cite:p:`Munk:1960uk`.
The loading change, such as the redistribution of ocean mass from tides or circulation, acts upon the *surface of the Earth* :cite:p:`Wahr:1998hy`.
These ":term:`load tides <Load Tide>`" are computed through a convolution of the tidal constituents and load Love numbers typically by either a :ref:`Green's function <greens-functions>` or :ref:`spherical harmonic <spherical-harmonics>` approach :cite:p:`Agnew:2013vx,Farrell:1972cm,Ray:1989wf`.
In either of those cases, the calculation uses a set of load Love numbers to high :ref:`spherical harmonic degree <spherical-harmonics>`.

.. role:: raw-html(raw)
   :format: html

.. plot:: ./background/load-love-numbers.py
    :align: center
    
    Elastic load :term:`Love/Shida numbers <Love and Shida Numbers>` from :cite:t:`Wang:2012gc` :raw-html:`<br />`
    :math:`k_l` and :math:`l_l` are both scaled by spherical harmonic degree :math:`l`

.. tip::

    Tidal Love numbers are for *external* potential perturbations and Load Love numbers are for *surface* mass changes :cite:p:`Michel:2021bn`

.. _greens-functions:

Green's Functions
-----------------

The point response of the Earth to a surface load change can be computed using Green's functions :cite:p:`Farrell:1972cm`. 
Convolving the Green's function for vertical displacement with the mass of an ocean tide height change provides the induced :term:`load tides <Load Tide>`.

.. math::
    :label: 10.1
    :name: eq:10.1

    u(\theta) &= \frac{a}{m_e}\sum_{l=0}^\infty h_l\,P_l(\cos{\theta})\\
              &= \frac{a h_\infty}{m_e}\sum_{l=0}^\infty P_l(\cos{\theta}) + \frac{a}{m_e}\sum_{l=0}^\infty (h_l - h_\infty)\,P_l(\cos{\theta})\\
              &= \frac{a}{m_e}\left(\frac{h_\infty}{2\sin{(\theta/2)}} + \sum_{l=0}^\infty (h_l - h_\infty)\,P_l(\cos{\theta})\right)

.. role:: raw-html(raw)
   :format: html

.. plot:: ./background/greens-functions.py
    :align: center
    
    Green's function for elastic vertical deformation using :term:`Love numbers <Love and Shida Numbers>` from :cite:t:`Wang:2012gc` :raw-html:`<br />`
    Function is normalized as described in Table A3 of :cite:t:`Farrell:1972cm`
