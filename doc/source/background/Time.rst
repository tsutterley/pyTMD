Time
####

The Julian Day (JD) is the continuous count of days starting at noon on January 1, 4713 B.C (-4712-01-01T12:00:00).
The Modified Julian Day (MJD) differs from the Julian Day by reducing the number of digits for modern periods, and by beginning at midnight.
The MJD is calculated from the Julian Day by

.. math::
    :label: 5

    MJD = JD - 2400000.5

The start of the Modified Julian Day calendar is 1858-11-17T00:00:00.
Time in Julian centuries (36525 days) are calculated relative to noon on January 1, 2000 (2000-01-01T12:00:00).

.. math::
    :label: 6

    T = \frac{JD - 2451545.0}{36525}

Standards
---------

`Tables of leap seconds <https://github.com/tsutterley/timescale/blob/main/timescale/data/leap-seconds.list>`_ are used to convert between GPS, LORAN and TAI times.

- TAI time: International Atomic Time which is computed as the weighted average of several hundred atomic clocks.
- UTC time: Coordinated Universal Time which is `periodically adjusted <https://www.nist.gov/pml/time-and-frequency-division/leap-seconds-faqs>`_ to account for the difference between the definition of the second and the rotation of Earth. UTC is based off of atomic clocks and 1 day is exactly 86,400 seconds.
- GPS time: Atomic timing system for the Global Positioning System constellation of satellites monitored by the United States Naval Observatory (USNO). GPS time and UTC time were equal on January 6, 1980. TAI time is ahead of GPS time by 19 seconds.
- LORAN time: Atomic timing system for the Loran-C chain transmitter sites used in terrestrial radionavigation. LORAN time and UTC time were equal on January 1, 1958. TAI time is ahead of LORAN time by 10 seconds.

`Tables of delta times <https://github.com/tsutterley/timescale/blob/main/timescale/data/merged_deltat.data>`_ are used to convert between dynamic (TT) and universal (UT1) times :cite:p:`Meeus:1991vh`.
Universal Time (UT1) is based on the rotation of the Earth, which varies irregularly, and so UT1 is adjusted periodically.
Dynamic Time (TT) is a uniform, monotonically increasing time standard based on atomic clocks that is used for the accurate calculation of celestial mechanics, orbits and ephemerides.
Delta times (TT - UT1) can be added to Universal Time (UT1) values to convert to Dynamic Time (TT) values.
`Coordinated Universal Time (UTC) <https://crf.usno.navy.mil/ut1-utc>`_ is based on International Atomic Time (TAI) with leap seconds added to keep it within 0.9 seconds of UT1.
