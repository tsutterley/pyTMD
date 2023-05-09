==========
Background
==========

Ocean and Load Tides
####################

The rise and fall of the oceanic tides are a major source of the vertical variability of the ocean surface.
Ocean tides are typically observed using float gauges, GPS stations, gravimeters, tiltmeters, pressure recorders, and satellite altimeters.
For each of these measurements, it is important to note the `vertical datum of the measurement technique <https://www.esr.org/data-products/antarctic_tg_database/ocean-tide-and-ocean-tide-loading/>`_.
Ocean tides are driven by gravitational undulations due to the relative positions of the Earth, moon and sun, and the centripetal acceleration due to the Earth's rotation [Meeus1998]_.
A secondary tidal effect, known as load tides, is due to the elastic response of the Earth's crust to ocean tidal loading, which produces deformation of both the sea floor and adjacent land areas.
Tidal oscillations for both ocean and load tides can be decomposed into a series of tidal constituents (or partial tides) of particular frequencies.

Ocean tide models are typically one of following categories:
1) an empirically adjusted model,
2) a barotropic hydrodynamic model constrained by data assimilation,
or 3) an unconstrained hydrodynamic model [Stammer2014]_.
``pyTMD`` is not an ocean or load tide model, but rather a tool for using constituents from ocean and load tide models to calculate the tide deflections or currents at particular locations and times [Egbert2002]_.

Solid Earth Tides
#################

Similar to ocean tides, solid Earth tides are tidal deformations due to gravitational undulations based on the relative positions of the Earth, moon and sun [Meeus1998]_ [Montenbruck1989]_.
Basic ephemerides for the sun and moon can be calculated within ``pyTMD``, and high-resolution ephemerides for the sun and moon can be downloaded from the Jet Propulsion Laboratory.

Within ``pyTMD``, the tidal deformation of the Earth is modeled using a Load Love/Shida number formalism taking into account the effects of mantle anelasticity [Mathews1997]_ [Wahr1981]_.
The formalism for estimating the solid Earth tides is based upon `IERS Conventions <https://iers-conventions.obspm.fr/>`_.
Love and Shida numbers describe the elastic response of the Earth in terms of vertical displacement (*h*), gravitational potential (*k*) and horizontal displacement (*l*).

In addition to the ups and downs of tides, there is a considerable portion of tidal potential and displacement that does not vary in time, a *permanent tide* that is due to the Earth being in the presence of the Sun and Moon (and other planetary bodies).
The `Earth is lower in polar areas and higher in equatorial areas <https://www.ngs.noaa.gov/PUBS_LIB/EGM96_GEOID_PAPER/egm96_geoid_paper.html>`_ than it would without those gravitational effects.
The `IERS formalism <https://iers-conventions.obspm.fr/>`_ for determining station locations is to remove all cyclical and permanent components of the tides, which is known as a "tide-free" system.
This is the default tide-system within ``pyTMD``.
Alternatively, the permanent tide components can be added back in order to calculate the station locations in a "mean-tide" state.
The radial difference in terms of latitude between the mean-tide and tide-free systems is:

.. math::
    :label: 1

    \delta r(\varphi) = -0.120582 \left(\frac{3}{2} sin^2 \varphi - \frac{1}{2} \right)

Pole Tides
##########

Load and ocean pole tides are driven by variations in the Earth's figure axis (e.g. Chandler wobble and annual variations) [Wahr1985]_ [Desai2002]_.
These pole tides are due to Earth's ellipsoidal shape shifting as the rotation axis of the Earth
moves with respect to the mean pole location, and for the case of ocean pole tides the centripetal effects of polar motion on the ocean [Desai2002]_ [Desai2015]_.
The formalism for estimating the pole tides is also based upon `IERS Conventions <https://iers-conventions.obspm.fr/>`_.
The Earth Orientation Parameters (EOPs) necessary for computing load pole and ocean pole tide variations are included within ``pyTMD``.
The currently accepted formalism for estimating the reference position of the Earth's figure axis at a given date is the `IERS 2018 secular pole model <https://iers-conventions.obspm.fr/chapter7.php>`_:

.. math::
    :label: 2

    \bar{x}_s(t) &= 0.055 + 0.001677(t - 2000.0)\\
    \bar{y}_s(t) &= 0.3205 + 0.00346(t - 2000.0)


The time-dependent offsets from the reference rotation pole position, are then calculated using instantaneous values of the Earth Orientation Parameters.


.. math::
    :label: 3

    m_1(t) &= x_p(t) - \bar{x}_s(t)\\
    m_2(t) &= -(y_p(t) - \bar{y}_s(t))

Terrestrial Reference Systems
#############################

Locations of planetary bodies and satellites can be determined in an Earth-centered Earth-Fixed (ECEF) coordinate system [Montenbruck1989]_.
ECEF is a Cartesian coordinate system representing *X*, *Y*, and *Z* measurements from the Earth's center of mass.
The *Z* axis is aligned with the Earth's rotation axis, the *X* axis is aligned with the intersection of the prime meridian and the equator, and the *Y* axis is aligned with 90 degrees east longitude and the equator.

As opposed to simple vertical offsets, changing the terrestial reference system can involve both `translation and rotation of the reference system <https://itrf.ign.fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt>`_.
This involves converting from a geographic coordinate system into a Cartesian coordinate system.
Within ``pyTMD``, solid Earth tides are calculated using ECEF coordinates, and pole tides are calculated using geocentric coordinates.

References
##########

.. [Desai2002] S Desai, "Observing the pole tide with satellite altimetry", *Journal of Geophysical Research: Oceans*, 107(C11), (2002). `doi: 10.1029/2001JC001224 <https://doi.org/10.1029/2001JC001224>`_

.. [Desai2015] S Desai, J Wahr and B Beckley "Revisiting the pole tide for and from satellite altimetry", *Journal of Geodesy*, 89(12), p1233-1243, (2015). `doi: 10.1007/s00190-015-0848-7 <https://doi.org/10.1007/s00190-015-0848-7>`_

.. [Egbert2002] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of Barotropic Ocean Tides", *Journal of Atmospheric and Oceanic Technology*, 19(2), 183--204, (2002). `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

.. [Mathews1997] P. M. Mathews, V. Dehant and J. M. Gipson, "Tidal station displacements", *Journal of Geophysical Research: Solid Earth*, 102(B9), 20469--20477, (1997). `doi: 10.1029/97JB01515 <https://doi.org/10.1029/97JB01515>`_

.. [Meeus1998] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).

.. [Montenbruck1989] O. Montenbruck, *Practical Ephemeris Calculations*, 146 pp., (1989).

.. [Stammer2014] D. Stammer et al., "Accuracy assessment of global barotropic ocean tide models", *Reviews of Geophysics*, 52, 243--282, (2014). `doi: 10.1002/2014RG000450 <https://doi.org/10.1002/2014RG000450>`_

.. [Wahr1981] J. M. Wahr, "Body tides on an elliptical, rotating, elastic and oceanless Earth", *Geophysical Journal of the Royal Astronomical Society*, 64(3), 677--703, (1981). `doi: 10.1111/j.1365-246X.1981.tb02690.x <https://doi.org/10.1111/j.1365-246X.1981.tb02690.x>`_

.. [Wahr1985] J. M. Wahr, "Deformation induced by polar motion", *Journal of Geophysical Research: Solid Earth*, 90(B11), 9363--9368, (1985). `doi: 10.1029/JB090iB11p09363 <https://doi.org/10.1029/JB090iB11p09363>`_

.. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
