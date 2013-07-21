"""Astronomical module.

Some of the moon astronomical functions here adapted from
Brett Hamilton's Astro-MoonPhase Perl module, available at the time
of writing at:
  http://cpansearch.perl.org/src/BRETT/Astro-MoonPhase-0.60/MoonPhase.pm

Other planetary data and mathematical formulae obtained from NASA's
Jet Propulsion Laboratory at:
  http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf

Much useful information and assistance also gained from:
  http://www.stargazing.net/kepler/ellipse.html

"""


from __future__ import division

from datetime import datetime, tzinfo, timedelta
from math import radians, degrees, sin, cos, atan, sqrt, floor
from collections import namedtuple


# Non-public named tuples

# Disable pylint message about constant names, we use named tuples
# like classes, so we'll follow the normal naming convention for those.
#
# pylint: disable=C0103

_ZodiacInfo = namedtuple("_ZodiacInfo", ["longitude",
                                         "sign_index",
                                         "sign_name",
                                         "sign_short_name",
                                         "degrees",
                                         "minutes",
                                         "seconds"])

_OrbitalElements = namedtuple("_OrbitalElements", ["semi_major_axis",
                                                   "eccentricity",
                                                   "inclination",
                                                   "mean_longitude",
                                                   "longitude_perihelion",
                                                   "longitude_asc_node"])

_HMS = namedtuple("_HMS", ["hours", "minutes", "seconds"])

_DMS = namedtuple("_DMS", ["degrees", "minutes", "seconds"])

SphCoords = namedtuple("SphCoords", ["right_ascension",
                                     "declination",
                                     "distance"])

RectCoords = namedtuple("RectCoords", ["x", "y", "z"])

# pylint: enable=C0103


# Non-public constants

# NOTE: Data for Keplerian elements obtained from:
#   http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf

_OE_J2000 = {
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

_ZODIAC_SIGNS = [
    "Aries", "Taurus", "Gemini", "Cancer",
    "Leo", "Virgo", "Libra", "Scorpio",
    "Sagittarius", "Capricorn", "Aquarius", "Pisces"
]

_ZODIAC_SIGNS_SHORT = [
    "AR", "TA", "GE", "CN", "LE", "VI",
    "LI", "SC", "SG", "CP", "AQ", "PI"
]

_J2000 = 2451545.0
_SECONDS_IN_A_DAY = 86400.0


# UTC class needed for datetime functions

class UTC(tzinfo):

    """UTC tzinfo class, needed for constructing datetime objects."""

    # Disable pylint message about unused arguments
    # pylint: disable=W0613

    def utcoffset(self, dtime):
        return timedelta(0)

    def tzname(self, dtime):
        return "UTC"

    def dst(self, dtime):
        return timedelta(0)

    # pylint: enable=W0613


# Non-public functions

def _normalize_degrees(angle):

    """Normalizes a given angle in degrees so that
    it returns 0 <= angle < 360.

    Arguments:
    angle -- the angle, in degrees, to normalize

    Returns: the normalized angle.

    """

    return angle - 360 * floor(angle / 360.0)


def _deg_to_hms(degs):

    """Converts the supplied angle in degrees to
    a string formatted as hours, minutes and seconds.

    Arguments:
    degs -- the degrees to format

    """

    degs = _normalize_degrees(degs)
    total_seconds = degs / 360.0 * _SECONDS_IN_A_DAY
    hours = int(total_seconds / 3600)
    minutes = int((total_seconds - hours * 3600) / 60)
    seconds = int(round(total_seconds - hours * 3600 - minutes * 60))

    return _HMS(hours, minutes, seconds)


def _deg_to_dms(degs):

    """Converts the supplied angle in degrees to a
    a string formatted as degrees, minutes and seconds.

    Arguments:
    degs -- the degrees to format

    """

    total_seconds = degs * 3600
    degs = int(total_seconds / 3600)
    minutes = int((total_seconds - degs * 3600) / 60)
    seconds = int(round(total_seconds - degs * 3600 - minutes * 60))

    return _DMS(degs, minutes, seconds)


def _get_zodiac_info(rasc):

    """Returns a populated _ZodiacInfo named tuple for
    a given right ascension, in degrees.

    Arguments:
    rasc -- a right ascension, in degrees

    Returns: an appropriate _ZodiacInfo named tuple

    """

    dms = _deg_to_dms(rasc)
    sign_index = dms.degrees // 30
    sign_short_name = _ZODIAC_SIGNS_SHORT[sign_index]
    sign_name = _ZODIAC_SIGNS[sign_index]

    return _ZodiacInfo(rasc, sign_index, sign_name, sign_short_name,
                       dms.degrees % 30, dms.minutes, dms.seconds)


def _get_current_utc_datetime():

    """Returns an aware datetime object for the current UTC time."""

    dti = datetime.utcnow()
    return datetime(dti.year, dti.month, dti.day,
                    dti.hour, dti.minute, dti.second,
                    dti.microsecond, UTC())


# Public functions

def julian_date(dtime=None):

    """Returns a Julian Date for a given datetime object,
    or for the current time if one is not provided.

    Arguments:
    dtime -- a datetime object containing the desired time

    Returns: the corresponding Julian Date

    """

    # Get datetime object for our selected Julian Date epoch,
    # and for the current datetime if one has not been specified.

    epoch_dt = datetime(2000, 1, 1, 12, 0, 0, 0, UTC())
    if dtime is None:
        dtime = _get_current_utc_datetime()

    # Get a timedelta object between the desired time
    # and our selected Julian Date epoch, and express
    # the timedelta in fractional days.

    timediff = dtime - epoch_dt
    days_diff = timediff.days + timediff.seconds / _SECONDS_IN_A_DAY

    # Return the fractional days plus the epoch

    return days_diff + _J2000


def kepler(m_anom, eccentricity):

    """Solves Kepler's equation for the eccentric anomaly
    using Newton's method.

    Arguments:
    m_anom -- mean anomaly in radians
    eccentricity -- eccentricity of the ellipse

    Returns: eccentric anomaly in radians.

    """

    desired_accuracy = 1e-6
    e_anom = m_anom

    while True:
        diff = e_anom - eccentricity * sin(e_anom) - m_anom
        e_anom -= diff / (1 - eccentricity * cos(e_anom))
        if abs(diff) <= desired_accuracy:
            break
    return e_anom


def rec_to_sph(rcds):

    """Transforms the provided rectangular coordinates into spherical
    coordinates.

    Arguments:
    rcds -- a RectCoords named tuple containing the rectangular coordinates:
      x -- rectangular x-axis coordinate (to vernal equinox)
      y -- rectangular y-axis coordinate (90 degrees east in
           the plane of the equator)
      z -- rectangular z-axis coodinrate (to north pole)

    Returns: a SphCoords named tuple consisting of:
      rasc -- right ascension, in degrees
      decl -- declination, in degrees
      dist -- distance, in astronomical units

    """

    rasc = degrees(atan(rcds.y / rcds.x))
    if rcds.x < 0:
        rasc += 180
    elif rcds.y < 0:
        rasc += 360
    decl = degrees(atan(rcds.z / sqrt(rcds.x ** 2 + rcds.y ** 2)))
    dist = sqrt(rcds.x ** 2 + rcds.y ** 2 + rcds.z ** 2)

    return SphCoords(rasc, decl, dist)


def zodiac_sign(rasc, short=False):

    """Returns a string containing the zodiac sign of
    a specified right ascension given in degrees.

    Arguments:
    rasc -- the right ascension to convert
    short -- set to True to show short two letter abbreviations
    rather than the full sign name, e.g. AR instead of Aries

    """

    zinfo = _get_zodiac_info(rasc)
    if short:
        return zinfo.sign_short_name
    else:
        return zinfo.sign_name


def rasc_to_zodiac(rasc, seconds=False):

    """Returns a string showing a right ascension in terms
    of the zodiac, e.g. 20CN15, meaning the 15th minute of the
    20th degree in the sign of Cancer.

    Arguments:
    rasc -- the right ascension to convert
    seconds -- set to True to also show the seconds, e.g. 20CN15'33

    """

    zinfo = _get_zodiac_info(rasc)
    if seconds:
        second_string = "'{0:02}".format(zinfo.seconds)
    else:
        second_string = ""
    return "{0:02}{1}{2:02}{3}".format(zinfo.degrees,
                                       zinfo.sign_short_name,
                                       zinfo.minutes,
                                       second_string)


def rasc_string(rasc):

    """Converts the supplied right ascension, in degrees, to
    a string formatted as hours, minutes and seconds.

    Arguments:
    rasc -- the right ascension to format

    """

    hms = _deg_to_hms(rasc)

    return "{0:02}h {1:02}m {2:02}s".format(hms.hours,
                                            hms.minutes, hms.seconds)


def decl_string(decl):

    """Converts the supplied declination, in degrees, to
    a string formatted as degrees, minutes and seconds.

    Arguments:
    decl -- the declination to format

    """

    dms = _deg_to_dms(decl)

    return "{0:+03}d {1:02}m {2:02}s".format(dms.degrees,
                                             abs(dms.minutes),
                                             abs(dms.seconds))


# Classes

class _Planet(object):

    """Abstract planet class."""

    def __init__(self, dtime=None):

        """Initializer."""

        # Get current datetime if one not provided

        self._dtime = dtime
        if self._dtime is None:
            self._dtime = _get_current_utc_datetime()

        # Get and store Keplerian elements

        self._oes = self._get_orbital_elements()

    def _get_orbital_elements(self):

        """Returns the Keplerian elements for the datetime
        specified at object construction.

        """

        # Calculate the number of Julian centuries since J2000

        jdc = (julian_date(self._dtime) - _J2000) / 36525

        # Calculate and return Keplerian elements

        # Disable pylint warning for abstract _Planet class
        # having no _pname member.
        #
        # pylint: disable=E1101

        tmp = []
        for elem in range(6):
            tmp.append(_OE_J2000[self._pname][elem] +
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
        any heliocentric coordinates for Earth will be zero.

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


class Moon(object):

    """Moon class."""

    def __init__(self):

        """Initializer."""

        self._m_long_epoch = 64.975464
        self._m_long_perigree = 349.383063
        self._m_long_node = 151.950429
        self._inclination = 5.145396
        self._eccentricity = 0.054900
        self._ang_size = 0.5181
        self._sm_axis = 384401.0
        self._parallax = 0.9507


# main function

def main():

    """Main function."""

    print "Current planetary data:\n"
    print "PLANET    R.ASCENSION   DECLINATION  DIST (AU) ZODIAC ZODIAC SIGN"
    print "=======   ===========  ============= ========= ====== ==========="

    for planet in [Sun, Mercury, Venus, Mars, Jupiter,
                   Saturn, Uranus, Neptune, Pluto]:
        plnt = planet()
        print "{0:8}: {1}, {2}, {3:9.6f} {4} {5}".format(
                plnt.name().capitalize(),
                plnt.right_ascension(formatted=True),
                plnt.declination(formatted=True),
                plnt.distance(),
                plnt.right_ascension(zodiac=True),
                plnt.zodiac_sign()
        )


# Direct entry point

if __name__ == "__main__":
    main()
