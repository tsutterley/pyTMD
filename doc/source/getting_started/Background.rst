==========
Background
==========

Ocean and Load Tides
####################

The rise and fall of the oceanic tides are a major source of the vertical variability of the ocean surface.
Ocean tides are typically observed using float gauges, GPS stations, gravimeters, tiltmeters, pressure recorders, and satellite altimeters.
For each of these measurements, it is important to note the `vertical datum of the measurement technique <https://www.esr.org/data-products/antarctic_tg_database/ocean-tide-and-ocean-tide-loading/>`_.
Ocean tides are driven by gravitational undulations due to the relative positions of the Earth, moon and sun, and the centripetal acceleration due to the Earth's rotation :cite:p:`Doodson:1921kt` :cite:p:`Meeus:1991vh`.
A secondary tidal effect, known as load tides, is due to the elastic response of the Earth's crust to ocean tidal loading, which produces deformation of both the sea floor and adjacent land areas.
Tidal oscillations for both ocean and load tides can be decomposed into a series of tidal constituents (or partial tides) of particular frequencies that are associated with the relative positions of the sun, moon and Earth.
These tidal constituents are typically classified into different "species" based on their approximate period: short-period, semi-diurnal, diurnal, and long-period.

The amplitude and phase of major constituents are provided by ocean tide models, which can be used for tidal predictions.
Ocean tide models are typically one of following categories:
1) empirically adjusted models,
2) barotropic hydrodynamic models constrained by data assimilation, and
3) unconstrained hydrodynamic models :cite:p:`Stammer:2014ci`.
``pyTMD`` is not an ocean or load tide model, but rather a tool for using constituents from ocean and load tide models to calculate the tide deflections or currents at particular locations and times :cite:p:`Egbert:2002ge`.

``pyTMD.io`` contains routines for reading major constituent values from commonly available tide models, and interpolating those values to spatial locations.
For any given time, ``pyTMD.astro`` can calculate the longitudes of the sun (`S`), moon (`H`), lunar perigree (`P`), ascending lunar node (`N`) and solar perigree (`Ps`), which are used in combination with the lunar hour angle (\ |tau|\ ) in a six-dimensional Fourier series :cite:p:`Doodson:1921kt` :cite:p:`Dietrich:1980ua`.
Each constituent has a particular "Doodson number" describing the polynomial coefficients of each of these astronomical terms in the Fourier series :cite:p:`Doodson:1921kt`. 
``pyTMD`` stores these coefficients in an easily accessible `JSON database <https://github.com/tsutterley/pyTMD/blob/main/pyTMD/data/doodson.json>`_ supplied with the program.
Together these coefficients can be used to calculate the frequencies and 18.6-year modulations of the tidal constituents, and allow for the accurate determination of the tidal amplitudes :cite:p:`Schureman:1958ty` :cite:p:`Dietrich:1980ua`.

Solid Earth Tides
#################

Similar to ocean tides, solid Earth tides (or body tides) are tidal deformations due to gravitational undulations based on the relative positions of the Earth, moon and sun :cite:p:`Agnew:2015kw` :cite:p:`Doodson:1921kt` :cite:p:`Meeus:1991vh` :cite:p:`Montenbruck:1989uk`.
However, while ocean tides are apparent to observers on the coast, solid Earth tides are typically more difficult to observe due to the reference frame of the observer moving.
The tidal deformation of the Earth is to a very high degree instantaneous, with the Earth's response to the gravitational potential of the moon and sun being nearly immediate.
The total gravitational potential at a position on the Earth's surface due to a celestial object is directly related to the distance between the Earth and the object, and the mass of that object :cite:p:`Agnew:2015kw` :cite:p:`Wahr:1981ea`.
Analytical approximate positions for the sun and moon can be calculated within ``pyTMD``, and high-resolution numerical ephemerides for the sun and moon can be downloaded from the `Jet Propulsion Laboratory <https://ssd.jpl.nasa.gov/planets/orbits.html>`_.

Within ``pyTMD``, the tidal deformation of the Earth is modeled using the Load Love/Shida numbers formalism described in the `IERS Conventions <https://iers-conventions.obspm.fr/>`_, which are based on :cite:p:`Mathews:1997js`.
Love and Shida numbers describe the elastic response of the Earth in terms of vertical displacement (*h*), gravitational potential (*k*) and horizontal displacement (*l*).
For a spherical, non-rotating Earth, the Love and Shida numbers are largely independent of tidal frequency :cite:p:`Wahr:1981ea`.
However, for a rotating, ellipsoidal Earth, the Love and Shida numbers are dependent on tidal frequency, with resonances in the diurnal and semi-diurnal bands :cite:p:`Wahr:1981ea`.
``pyTMD`` computes these frequency-dependent corrections along with the dissipative mantle anelasticity corrections following :cite:p:`Mathews:1997js`.

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

Load and ocean pole tides are driven by variations in the Earth's figure axis (e.g. Chandler wobble and annual variations) :cite:p:`Wahr:1985gr` :cite:p:`Desai:2002ev` :cite:p:`Agnew:2015kw`.
These pole tides are due to Earth's ellipsoidal shape shifting as the rotation axis of the Earth
moves with respect to the mean pole location, and for the case of ocean pole tides the centripetal effects of polar motion on the ocean :cite:p:`Desai:2002ev` :cite:p:`Desai:2015jr`.
The formalism for estimating the pole tides is also based upon `IERS Conventions <https://iers-conventions.obspm.fr/>`_.
``pyTMD`` uses the ``timescale`` library for reading the Earth Orientation Parameters (EOPs) necessary for computing load pole and ocean pole tide variations.
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

Locations of planetary bodies and satellites can be determined in an Earth-centered Earth-Fixed (ECEF) coordinate system :cite:p:`Montenbruck:1989uk`.
ECEF is a Cartesian coordinate system representing *X*, *Y*, and *Z* measurements from the Earth's center of mass.
The *Z* axis is aligned with the Earth's rotation axis, the *X* axis is aligned with the intersection of the prime meridian and the equator, and the *Y* axis is aligned with 90 degrees east longitude and the equator.

As opposed to simple vertical offsets, changing the terrestial reference system can involve both `translation and rotation of the reference system <https://itrf.ign.fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt>`_.
This involves converting from a geographic coordinate system into a Cartesian coordinate system.
Within ``pyTMD``, solid Earth tides are calculated using ECEF coordinates, and pole tides are calculated using geocentric coordinates.

Nutation is the periodic oscillation of the Earth's rotation axis around its mean position.
Nutation is often split into two components, the nutation in longitude and the nutation in obliquity.
The angle between the equator and the orbital plane of Earth around the Sun (the ecliptic) defines the inclination of the Earth's rotation axis (obliquity of the ecliptic).

Time
####

The Julian Day (JD) is the continuous count of days starting at noon on January 1, 4713 B.C (-4712-01-01T12:00:00).
The Modified Julian Day (MJD) differs from the Julian Day by reducing the number of digits for modern periods, and by beginning at midnight.
The MJD is calculated from the Julian Day by

.. math::
    :label: 4

    MJD = JD - 2400000.5

The start of the Modified Julian Day calendar is 1858-11-17T00:00:00.
Time in Julian centuries (36525 days) are calculated relative to noon on January 1, 2000 (2000-01-01T12:00:00).

.. math::
    :label: 5

    T = \frac{JD - 2451545.0}{36525}

.. |tau|    unicode:: U+1D70F .. MATHEMATICAL ITALIC SMALL TAU
