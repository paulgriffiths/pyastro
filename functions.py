"""General functions for astronomical module.

Library Release 1.0

Copyright 2013 Paul Griffiths
Email: mail@paulgriffiths.net

Distributed under the terms of the GNU General Public License.
http://www.gnu.org/licenses/

"""


from __future__ import division
from __future__ import print_function

from datetime import datetime, tzinfo, timedelta
from math import degrees, sin, cos, atan, sqrt, floor, atan2, hypot
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

_HMS = namedtuple("_HMS", ["hours", "minutes", "seconds"])

_DMS = namedtuple("_DMS", ["degrees", "minutes", "seconds"])

SphCoords = namedtuple("SphCoords", ["right_ascension",
                                     "declination",
                                     "distance"])

RectCoords = namedtuple("RectCoords", ["x", "y", "z"])

# pylint: enable=C0103


# Non-public constants

_ZODIAC_SIGNS = [
    "Aries", "Taurus", "Gemini", "Cancer",
    "Leo", "Virgo", "Libra", "Scorpio",
    "Sagittarius", "Capricorn", "Aquarius", "Pisces"
]

_ZODIAC_SIGNS_SHORT = [
    "AR", "TA", "GE", "CN", "LE", "VI",
    "LI", "SC", "SG", "CP", "AQ", "PI"
]

J2000 = 2451545.0
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

    dms = _deg_to_dms(_normalize_degrees(rasc))
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

    return days_diff + J2000


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

    rasc = degrees(atan2(rcds.y, rcds.x))
    decl = degrees(atan(rcds.z / hypot(rcds.x, rcds.y)))
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
