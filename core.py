"""Planet classes for astronomical module.

Library Release 1.1

Copyright 2013 Paul Griffiths
Email: mail@paulgriffiths.net

Distributed under the terms of the GNU General Public License.
http://www.gnu.org/licenses/

The calculations for all planets, excluding the moon, are valid the 
time period 1800 through 2050, and are subject to the following
approximate errors:

  -- Right Ascension, +/- 50 arcsecs for all planets except Jupiter and
     Saturn, +/- 400 arcsecs for Jupiter, +/- 600 arcsecs for Saturn
  -- Declination, +/- 8 arcsecs for all planets except Jupiter and
     Saturn, +/- 10 arcsecs for Jupiter, +/- 25 arcsecs for Saturn
  -- Distance, +/- 25 km for the inner planets, and +/- between 200
     and 1,500 km for Jupiter through Pluto

Right ascension and declination calculations for the moon are accurate to
within approximately 4 arcminutes. 

Some of the moon astronomical functions here adapted from
Brett Hamilton's Astro-MoonPhase Perl module, available at the time
of writing at:
  http://cpansearch.perl.org/src/BRETT/Astro-MoonPhase-0.60/MoonPhase.pm

Other planetary data and mathematical formulae obtained from NASA's
Jet Propulsion Laboratory at:
  http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf

Much useful information and assistance also gained from
Keith Burnett's pages at:
  http://www.stargazing.net/kepler/ellipse.html
  http://www.stargazing.net/kepler/moon2.html

and from Paul Schlyter's page at:
  http://www.stjarnhimlen.se/comp/ppcomp.html

"""


from __future__ import division
from __future__ import print_function

from datetime import datetime
from math import radians, sin, cos, sqrt, atan2, hypot
from collections import namedtuple

from functions import J2000, RectCoords, _SECONDS_IN_A_DAY
from functions import _get_current_utc_datetime
from functions import julian_date, rec_to_sph, rasc_string, decl_string
from functions import kepler, zodiac_sign, UTC, rasc_to_zodiac


# Non-public named tuples

# Disable pylint message about constant names, we use named tuples
# like classes, so we'll follow the normal naming convention for those.
#
# pylint: disable=C0103

_OrbitalElements = namedtuple("_OrbitalElements", ["semi_major_axis",
                                                   "eccentricity",
                                                   "inclination",
                                                   "mean_longitude",
                                                   "longitude_perihelion",
                                                   "longitude_asc_node"])

# pylint: enable=C0103


# Non-public constants

# Keplerian element data for all planets excluding the moon, based
# on a J2000 epoch and with adjustment factors per Julian Century.
# Data obtained from:
#   http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf

_OEJ2000 = {
    "sun": _OrbitalElements(0, 0, 0, 0, 0, 0),
    "mercury": _OrbitalElements(0.387009927, 0.20563593, 7.00497902,
                                252.25032350, 77.45779628, 48.33076593),
    "venus": _OrbitalElements(0.72333566, 0.00677672, 3.39467605,
                              181.97909950, 131.60246718, 76.67984255),
    "embary": _OrbitalElements(1.00000261, 0.01671123, -0.00001531,
                               100.46457166, 102.93768193, 0.0),
    "mars": _OrbitalElements(1.52371034, 0.09339410, 1.84969142,
                             -4.55343205, -23.94362959, 49.55953891),
    "jupiter": _OrbitalElements(5.20288700, 0.04838624, 1.30439695,
                                34.39644051, 14.72847983, 100.47390909),
    "saturn": _OrbitalElements(9.53667594, 0.05386179, 2.48599187,
                               49.95424423, 92.59887831, 113.66242448),
    "uranus": _OrbitalElements(19.18916464, 0.04725744, 0.77263783,
                               313.23810451, 170.95427630, 74.01692503),
    "neptune": _OrbitalElements(30.06992276, 0.00859048, 1.77004347,
                                -55.12002969, 44.96476227, 131.78422574),
    "pluto": _OrbitalElements(39.48211675, 0.24882730, 17.14001206,
                              238.92903833, 224.06891629, 110.30393684)
}

_OE_CENT = {
    "sun": _OrbitalElements(0, 0, 0, 0, 0, 0),
    "mercury": _OrbitalElements(0.00000037, 0.00001906, -0.00594749,
                                149472.67411175, 0.16047689, -0.12534081),
    "venus": _OrbitalElements(0.00000390, -0.00004107, -0.00078890,
                              58517.81538729, 0.00268329, -0.27769418),
    "embary": _OrbitalElements(0.00000562, -0.00004392, -0.01294668,
                               35999.37244981, 0.32327364, 0.0),
    "mars": _OrbitalElements(0.00001847, 0.00007882, -0.00813131,
                             19140.30268499, 0.44441088, -0.29257343),
    "jupiter": _OrbitalElements(-0.00011607, -0.00013253, -0.00183714,
                                3034.74612775, 0.21252668, 0.20469106),
    "saturn": _OrbitalElements(-0.00125060, -0.00050991, 0.00193609,
                               1222.49362201, -0.41897216, -0.28867794),
    "uranus": _OrbitalElements(-0.00196176, -0.00004397, -0.00242939,
                               428.48202785, 0.40805281, 0.04240589),
    "neptune": _OrbitalElements(0.00026291, 0.00005105, 0.00035372,
                                218.45945325, -0.32241464, -0.00508664),
    "pluto": _OrbitalElements(-0.00031596, 0.00005170, 0.00004818,
                              145.20780515, -0.04062942, -0.01183482)
}

# Keplerian element data for the moon and the sun (when used for
# the moon calculations) based on a December 31, 1999 epoch with
# adjustment factors per day.
# Data obtained from:
#   http://www.stjarnhimlen.se/comp/ppcomp.html

_OE_Y2000 = {
    "moon": _OrbitalElements(60.2666, 0.0549, 5.1454, 198.5516,
                             83.1862, 125.1228),
    "sun": _OrbitalElements(1, 0.016709, 0, 278.9874, -77.0596, 0)
}

_OE_DAY = {
    "moon": _OrbitalElements(0, 0, 0, 13.1763964649, 0.111403514,
                             -0.0529538083),
    "sun": _OrbitalElements(0, -0.000000001151, 0, 0.98564735200,
                            0.00004709350, 0)
}

# Other non-public constants

_ZODIAC_SIGNS = [
    "Aries", "Taurus", "Gemini", "Cancer",
    "Leo", "Virgo", "Libra", "Scorpio",
    "Sagittarius", "Capricorn", "Aquarius", "Pisces"
]

_ZODIAC_SIGNS_SHORT = [
    "AR", "TA", "GE", "CN", "LE", "VI",
    "LI", "SC", "SG", "CP", "AQ", "PI"
]


# Classes

class _Planet(object):

    """Abstract planet class."""

    def __init__(self, dtime=None):

        """Initializer."""

        # Get current datetime if one not provided

        self._dtime = dtime
        if self._dtime is None:
            self._dtime = _get_current_utc_datetime()

        # Initialize _rhc for moon calculations

        self._rhc = 0

        # Get and store Keplerian elements

        self._oes = self._get_orbital_elements()

    def _get_orbital_elements(self):

        """Returns the Keplerian elements for the datetime
        specified at object construction.

        """

        # Calculate the number of Julian centuries since J2000

        jdc = (julian_date(self._dtime) - J2000) / 36525

        # Calculate and return Keplerian elements

        # Disable pylint warning for abstract _Planet class
        # having no _pname member.
        #
        # pylint: disable=E1101

        tmp = []
        for elem in range(6):
            tmp.append(_OEJ2000[self._pname][elem] +
                       _OE_CENT[self._pname][elem] * jdc)

        # pylint: enable=E1101

        return _OrbitalElements(*tmp)            # pylint: disable=W0142

    def helio_ecl_coords(self):

        """Returns a three-element tuple containing the
        heliocentric ecliptic rectangular coordinates of the planet.

        """

        # Step 1: get Keplerian elements

        oea = self._oes.semi_major_axis
        oee = self._oes.eccentricity
        oei = radians(self._oes.inclination)
        oel = radians(self._oes.mean_longitude)
        oew = radians(self._oes.longitude_perihelion)
        oeo = radians(self._oes.longitude_asc_node)

        # Step 2: calculate argument of perihelion and mean anomaly

        arp = oew - oeo
        m_anom = oel - oew

        # Step 3: obtain eccentric anomaly by solving Kepler

        e_anom = kepler(m_anom, oee)

        # Step 4: Compute heliocentric coordinates in the orbital plane

        xhc = oea * (cos(e_anom) - oee)
        yhc = oea * sqrt(1 - oee ** 2) * sin(e_anom)
        self._rhc = hypot(xhc, yhc)

        # Step 5: compute ecliptic coordinates in J2000 ecliptic plane

        xec = ((cos(arp) * cos(oeo) - sin(arp) * sin(oeo) * cos(oei)) * xhc +
               (-sin(arp) * cos(oeo) - cos(arp) * sin(oeo) * cos(oei)) * yhc)
        yec = ((cos(arp) * sin(oeo) + sin(arp) * cos(oeo) * cos(oei)) * xhc +
               (-sin(arp) * sin(oeo) + cos(arp) * cos(oeo) * cos(oei)) * yhc)
        zec = sin(arp) * sin(oei) * xhc + cos(arp) * sin(oei) * yhc

        return RectCoords(xec, yec, zec)

    def geo_ecl_coords(self):

        """Returns a three-element tuple containing the
        geocentric ecliptic rectangular coordinates of the planet.

        """

        phc = self.helio_ecl_coords()
        ehc = Earth(self._dtime).helio_ecl_coords()
        return RectCoords(phc.x - ehc.x, phc.y - ehc.y, phc.z - ehc.z)

    def geo_equ_coords(self):

        """Returns a three-element tuple containing the
        geocentric equatorial rectangular coordinates of the planet.

        """

        ecl = self.geo_ecl_coords()
        obliquity = radians(23.43928)

        xeq = ecl.x
        yeq = ecl.y * cos(obliquity) - ecl.z * sin(obliquity)
        zeq = ecl.y * sin(obliquity) + ecl.z * cos(obliquity)

        return RectCoords(xeq, yeq, zeq)

    def right_ascension(self, formatted=False, zodiac=False):

        """Returns the right ascension of the planet at the time
        specified at object construction.

        Arguments:
        formatted -- set to True to return the right ascension in
        a string of 'HHh MMm SSs' format, otherwise the number is
        returned in degrees.
        zodiac -- set to True to return the right ascension as a
        a string showing a right ascension in terms of the zodiac,
        e.g. 20CN15, meaning the 15th minute of the 20th degree in
        the sign of Cancer.

        Returns: either the right ascension in degrees, or in the
        specified string format.

        """

        sph_coords = rec_to_sph(self.geo_equ_coords())
        if formatted:
            return rasc_string(sph_coords.right_ascension)
        elif zodiac:
            return rasc_to_zodiac(sph_coords.right_ascension)
        else:
            return sph_coords.right_ascension

    def declination(self, formatted=False):

        """Returns the declination of the planet at the time specified
        at object construction.

        Arguments:
        formatted -- set to True to return the declination ascension in
        a string of '[+/-]DDd MMm SSs' format, otherwise the number is
        returned in degrees.

        Returns: either the declination in degrees, or the
        declination in a formatted DDMMSS string.

        """

        sph_coords = rec_to_sph(self.geo_equ_coords())
        if formatted:
            return decl_string(sph_coords.declination)
        else:
            return sph_coords.declination

    def distance(self):

        """Returns the distance of the planet at the time specified
        at object construction.

        Returns: the distance of the planet from Earth in astronomical
        units (AUs).

        """

        sph_coords = rec_to_sph(self.geo_equ_coords())
        return sph_coords.distance

    def zodiac_sign(self):

        """Returns a string containing the sign of the zodiac
        the planet is in at the time specified at object
        construction.

        """

        sph_coords = rec_to_sph(self.geo_equ_coords())
        return zodiac_sign(sph_coords.right_ascension)

    def name(self):

        """Returns the planet's name."""

        # Disable pylint warning about no _pname member
        # pylint: disable=E1101

        return self._pname

        # pylint: enable=E1101


class Mercury(_Planet):

    """Mercury class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "mercury"
        _Planet.__init__(self, dtime)


class Venus(_Planet):

    """Venus class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "venus"
        _Planet.__init__(self, dtime)


class Earth(_Planet):

    """Earth class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "embary"
        _Planet.__init__(self, dtime)

    def geo_ecl_coords(self):

        """Override this method to simply return zeros, since
        any geocentric coordinates for Earth will be zero.

        """

        return RectCoords(0, 0, 0)

    def geo_equ_coords(self):

        """Override this method to simply return zeros, since
        any geocentric coordinates for Earth will be zero.

        """

        return RectCoords(0, 0, 0)


class Mars(_Planet):

    """Mars class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "mars"
        _Planet.__init__(self, dtime)


class Jupiter(_Planet):

    """Jupiter class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "jupiter"
        _Planet.__init__(self, dtime)


class Saturn(_Planet):

    """Saturn class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "saturn"
        _Planet.__init__(self, dtime)


class Uranus(_Planet):

    """Uranus class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "uranus"
        _Planet.__init__(self, dtime)


class Neptune(_Planet):

    """Neptune class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "neptune"
        _Planet.__init__(self, dtime)


class Pluto(_Planet):

    """Pluto class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "pluto"
        _Planet.__init__(self, dtime)


class Sun(_Planet):

    """Sun class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "sun"
        _Planet.__init__(self, dtime)

    def helio_ecl_coords(self):

        """Override this method to simply return zeros, since
        any heliocentric coordinates for the Sun will be zero.

        """

        return RectCoords(0, 0, 0)


class Moon(_Planet):

    """Moon class."""

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "moon"
        _Planet.__init__(self, dtime)

    def _get_orbital_elements(self):

        """Returns the Keplerian elements for the datetime
        specified at object construction.

        Overriden for moon class, as the moon calculations use a
        different epoch and adjustment factor for calculating
        the Keplerian elements than is used for the other planets.

        """

        # Calculate the number of days since Y2000

        epoch_dt = datetime(1999, 12, 31, 0, 0, 0, 0, UTC())
        timediff = self._dtime - epoch_dt
        days = timediff.days + timediff.seconds / _SECONDS_IN_A_DAY

        # Calculate and return Keplerian elements

        # Disable pylint warning for abstract _Planet class
        # having no _pname member.
        #
        # pylint: disable=E1101

        tmp = []
        for elem in range(6):
            tmp.append(_OE_Y2000[self._pname][elem] +
                       _OE_DAY[self._pname][elem] * days)

        # pylint: enable=E1101

        return _OrbitalElements(*tmp)            # pylint: disable=W0142

    def geo_ecl_coords(self):

        """Returns a three-element tuple containing the
        geocentric ecliptic rectangular coordinates of the planet.

        Overriden for Moon class to adjust for perturbations. Note that
        for the Moon class, helio_ecl_coords() does not, in fact, return
        the heliocentric ecliptic rectangular coordinates, but the
        geocentric ecliptic rectangular coordinates unadjusted for
        perturbations.

        """

        hec = self.helio_ecl_coords()

        # Get latitude, longitude and rhc

        lon = atan2(hec.y, hec.x)
        lat = atan2(hec.z, hypot(hec.x, hec.y))
        rhc = self._rhc

        # Get selected Keplerian elements for the sun and the moon

        m_oel = radians(self._oes.mean_longitude)
        m_oew = radians(self._oes.longitude_perihelion)
        m_oeo = radians(self._oes.longitude_asc_node)

        # Disable pylint message about access to _get_orbital_elements()
        # pylint: disable=W0212

        soe = _SunForMoon(self._dtime)._get_orbital_elements()
        s_oel = radians(soe.mean_longitude)
        s_oew = radians(soe.longitude_perihelion)

        # pylint: enable=W0212

        # Calculate mean anomalies for the moon and the sun

        m_m_anom = m_oel - m_oew
        s_m_anom = s_oel - s_oew

        # Calculate mean elongation and argument of latitude of the moon

        m_mel = m_oel - s_oel
        m_arl = m_oel - m_oeo

        # Adjust for longitude perturbations

        dlon = -1.274 * sin(m_m_anom - 2 * m_mel)
        dlon += 0.658 * sin(2 * m_mel)
        dlon -= 0.186 * sin(s_m_anom)
        dlon -= 0.059 * sin(2 * m_m_anom - 2 * m_mel)
        dlon -= 0.057 * sin(m_m_anom - 2 * m_mel + s_m_anom)
        dlon += 0.053 * sin(m_m_anom + 2 * m_mel)
        dlon += 0.046 * sin(2 * m_mel - s_m_anom)
        dlon += 0.041 * sin(m_m_anom - s_m_anom)
        dlon -= 0.035 * sin(m_mel)
        dlon -= 0.031 * sin(m_m_anom + s_m_anom)
        dlon -= 0.015 * sin(2 * m_arl - 2 * m_mel)
        dlon += 0.011 * sin(m_m_anom - 4 * m_mel)
        lon = radians(dlon) + lon

        # Adjust for latitude perturbations

        dlat = -0.173 * sin(m_arl - 2 * m_mel)
        dlat -= 0.055 * sin(m_m_anom - m_arl - 2 * m_mel)
        dlat -= 0.046 * sin(m_m_anom + m_arl - 2 * m_mel)
        dlat += 0.033 * sin(m_arl + 2 * m_mel)
        dlat += 0.017 * sin(2 * m_m_anom + m_arl)
        lat = radians(dlat) + lat

        # Adjust for rhc pertubations

        rhc = rhc - 0.58 * cos(m_m_anom - 2 * m_mel)
        rhc -= 0.46 * cos(2 * m_mel)

        # Return cartersian coordinates of geocentric lunar position

        xgc = rhc * cos(lon) * cos(lat)
        ygc = rhc * sin(lon) * cos(lat)
        zgc = rhc * sin(lat)

        return RectCoords(xgc, ygc, zgc)


class _SunForMoon(Moon):

    """Sun class for moon.

    This class is provided purely to enable the calculations in
    the Moon class to calculate orbital elements for the Sun using
    the same methodology as used for its own orbital elements
    (i.e. using Dec 31, 1999 as an Epoch rather than J2000, and
    calculating the elements using a daily adjustment factor,
    rather than a Julian Century adjustment factor.

    """

    # Disable pylint messages about calling _Planet.__init__()
    # rather than Moon.__init__(), this is deliberate to avoid
    # overwriting self._pname which would happen if Moon.__init__()
    # was called. Moon.__init__() does nothing other than set its
    # own self._pname and call _Planet.__init__(), so no functionality
    # is lost.
    #
    # pylint: disable=W0233
    # pylint: disable=W0231

    def __init__(self, dtime=None):

        """Initializer."""

        self._pname = "sun"
        _Planet.__init__(self, dtime)

    # pylint: enable=W0233
    # pylint: enable=W0231


# Functions

def positions(dtime=None):

    """Prints information for all supported planets for a
    specified datetime.

    """

    if dtime is None:
        dtime = _get_current_utc_datetime()

    print("Current planetary data:\n")
    print("PLANET    R.ASCENSION   DECLINATION  DIST (AU)* ZODIAC ZODIAC SIGN")
    print("=======   ===========  ============= ========== ====== ===========")

    for planet in [Sun, Mercury, Venus, Mars, Jupiter,
                   Saturn, Uranus, Neptune, Pluto, Moon]:
        plnt = planet(dtime)
        print("{0:8}: {1}, {2}, {3:10.7f} {4} {5}".format(
                plnt.name().capitalize(),
                plnt.right_ascension(formatted=True),
                plnt.declination(formatted=True),
                plnt.distance(),
                plnt.right_ascension(zodiac=True),
                plnt.zodiac_sign()
        ))

    print("\n* Distance for the moon given in Earth radii.")
