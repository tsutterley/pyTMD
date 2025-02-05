Reference Systems
#################

Locations of planetary bodies and satellites can be determined in an Earth-centered Earth-Fixed (ECEF) coordinate system :cite:p:`Montenbruck:1989uk`.
ECEF is a Cartesian coordinate system representing *X*, *Y*, and *Z* measurements from the Earth's center of mass.
The *Z* axis is aligned with the Earth's rotation axis, the *X* axis is aligned with the intersection of the prime meridian and the equator, and the *Y* axis is aligned with 90 degrees east longitude and the equator.

As opposed to simple vertical offsets, changing the terrestial reference system can involve both `translation and rotation of the reference system <https://itrf.ign.fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt>`_.
This involves converting from a geographic coordinate system into a Cartesian coordinate system.
Within ``pyTMD``, solid Earth tides are calculated using ECEF coordinates, and pole tides are calculated using geocentric coordinates.

Nutation is the periodic oscillation of the Earth's rotation axis around its mean position.
Nutation is often split into two components, the nutation in longitude and the nutation in obliquity.
The angle between the equator and the orbital plane of Earth around the Sun (the ecliptic) defines the inclination of the Earth's rotation axis (obliquity of the ecliptic).
