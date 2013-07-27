"""Unit test module for coordinates functions."""

# Disable the following pylint warnings:
#  - long identifier names, these are deliberate for unittests
#  - too many public methods, supposed to be lots of tests here
#
# pylint: disable=C0103
# pylint: disable=R0904


import unittest
from datetime import datetime

from ..core import Mars, Earth, Venus, Mercury, Jupiter, Saturn
from ..core import Uranus, Neptune, Pluto, Sun, Moon, UTC
from ..functions import rasc_string, decl_string, rec_to_sph


class TestSequenceFunctions(unittest.TestCase):

    """Test sequence class for coordintes functions."""

    def setUp(self):
        pass

    def test_mars_helio_ecl_coords(self):

        """Test helio_ecl_coords() method for Mars."""

        # Test numbers taken from:
        #   http://www.stargazing.net/kepler/ellipse.html

        dtm = datetime(1997, 6, 21, 0, 0, 0, 0, UTC())

        test_xec = -1.18659376
        test_yec = -1.03200869
        test_zec = 0.00755328

        mars = Mars(dtm)
        xec, yec, zec = mars.helio_ecl_coords()

        self.assertAlmostEqual(xec, test_xec, places=3)
        self.assertAlmostEqual(yec, test_yec, places=3)
        self.assertAlmostEqual(zec, test_zec, places=3)

    def test_earth_helio_ecl_coords(self):

        """Test helio_ecl_coords() method for Earth."""

        # Test numbers taken from:
        #   http://www.stargazing.net/kepler/ellipse.html

        dtm = datetime(1997, 6, 21, 0, 0, 0, 0, UTC())

        test_xec = -0.00515816
        test_yec = -1.01625196
        test_zec = 0

        earth = Earth(dtm)
        xec, yec, zec = earth.helio_ecl_coords()

        self.assertAlmostEqual(xec, test_xec, places=3)
        self.assertAlmostEqual(yec, test_yec, places=3)
        self.assertAlmostEqual(zec, test_zec, places=3)

    def test_mars_geo_ecl_coords(self):

        """Test geo_ecl_coords() method for Mars."""

        # Test numbers taken from:
        #   http://www.stargazing.net/kepler/ellipse.html

        dtm = datetime(1997, 6, 21, 0, 0, 0, 0, UTC())

        test_xec = -1.18143560
        test_yec = -0.01575673
        test_zec = 0.00755328

        mars = Mars(dtm)
        xec, yec, zec = mars.geo_ecl_coords()

        self.assertAlmostEqual(xec, test_xec, places=3)
        self.assertAlmostEqual(yec, test_yec, places=3)
        self.assertAlmostEqual(zec, test_zec, places=3)

    def test_mars_geo_equ_coords(self):

        """Test geo_equ_coords() method for Mars."""

        # Test numbers taken from:
        #   http://www.stargazing.net/kepler/ellipse.html

        dtm = datetime(1997, 6, 21, 0, 0, 0, 0, UTC())

        test_xec = -1.18143560
        test_yec = -0.01746127
        test_zec = 0.00066224

        mars = Mars(dtm)
        xec, yec, zec = mars.geo_equ_coords()

        self.assertAlmostEqual(xec, test_xec, places=3)
        self.assertAlmostEqual(yec, test_yec, places=3)
        self.assertAlmostEqual(zec, test_zec, places=3)

    def test_mercury_geo_equ_coords(self):

        """Test geo_equ_coords() method for Mercury."""

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(1980, 5, 5, 20, 23, 0, 0, UTC())

        # Right ascension accurate to within four seconds
        # Declination accurate to within one minute
        # Distance accurate to three decimal places

        test_rasc = ["02h 20m 15s", "02h 20m 16s", "02h 20m 17s",
                     "02h 20m 18s", "02h 20m 19s", "02h 20m 20s",
                     "02h 20m 21s"]
        test_decl = "+12d 47m"
        test_dist = 1.30377991332952

        mercury = Mercury(dtm)
        rasc, decl, dist = rec_to_sph(mercury.geo_equ_coords())

        rasc_s = rasc_string(rasc)
        decl_s = decl_string(decl)

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)
        self.assertAlmostEqual(dist, test_dist, places=3)

    def test_venus_geo_equ_coords(self):

        """Test geo_equ_coords() method for Venus."""

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(1982, 6, 14, 8, 30, 0, 0, UTC())

        # Right ascension accurate to within two seconds
        # Declination accurate to within one minute
        # Distance accurate to three decimal places

        test_rasc = ["03h 00m 56s", "03h 00m 57s", "03h 00m 58s"]
        test_decl = "+15d 02m"
        test_dist = 1.23214680640560

        venus = Venus(dtm)
        rasc, decl, dist = rec_to_sph(venus.geo_equ_coords())

        rasc_s = rasc_string(rasc)
        decl_s = decl_string(decl)

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)
        self.assertAlmostEqual(dist, test_dist, places=3)

    def test_jupiter_geo_equ_coords(self):

        """Test geo_equ_coords() method for Jupiter."""

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(1991, 6, 17, 14, 00, 0, 0, UTC())

        # Right ascension accurate to within one second
        # Declination accurate to within one minute
        # Distance accurate to two decimal places

        test_rasc = ["08h 58m 28s"]
        test_decl = "+17d 55m"
        test_dist = 5.99194574022366

        jupiter = Jupiter(dtm)
        rasc, decl, dist = rec_to_sph(jupiter.geo_equ_coords())

        rasc_s = rasc_string(rasc)
        decl_s = decl_string(decl)

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)
        self.assertAlmostEqual(dist, test_dist, places=2)

    def test_saturn_geo_equ_coords(self):

        """Test geo_equ_coords() method for Saturn."""

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(1947, 12, 1, 12, 00, 0, 0, UTC())

        # Right ascension accurate to within four seconds
        # Declination accurate to within one minute
        # Distance accurate to two decimal places

        test_rasc = ["09h 44m 15s", "09h 44m 16s", "09h 44m 17s",
                     "09h 44m 18s", "09h 44m 19s", "09h 44m 20s",
                     "09h 44m 21s"]
        test_decl = "+14d 41m"
        test_dist = 8.84280809502593

        saturn = Saturn(dtm)
        rasc, decl, dist = rec_to_sph(saturn.geo_equ_coords())

        rasc_s = rasc_string(rasc)
        decl_s = decl_string(decl)

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)
        self.assertAlmostEqual(dist, test_dist, places=2)

    def test_uranus_geo_equ_coords(self):

        """Test geo_equ_coords() method for Uranus."""

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(1975, 10, 31, 8, 0, 0, 0, UTC())

        # Right ascension accurate to within four seconds
        # Declination accurate to within one minute
        # Distance accurate to two decimal places

        test_rasc = ["14h 05m 28s", "14h 05m 29s", "14h 05m 30s",
                     "14h 05m 31s", "14h 05m 32s", "14h 05m 33s",
                     "14h 05m 34s"]
        test_decl = "-12d 12m"
        test_dist = 19.487343031752

        uranus = Uranus(dtm)
        rasc, decl, dist = rec_to_sph(uranus.geo_equ_coords())

        rasc_s = rasc_string(rasc)
        decl_s = decl_string(decl)

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)
        self.assertAlmostEqual(dist, test_dist, places=2)

    def test_neptune_geo_equ_coords(self):

        """Test geo_equ_coords() method for Neptune."""

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(1966, 9, 15, 9, 0, 0, 0, UTC())

        # Right ascension accurate to within four seconds
        # Declination accurate to within one minute
        # Distance accurate to two decimal places

        test_rasc = ["15h 13m 44s", "15h 13m 45s", "15h 13m 46s",
                     "15h 13m 47s", "15h 13m 48s", "15h 13m 49s",
                     "15h 13m 50s"]
        test_decl = "-16d 12m"
        test_dist = 30.8456979309321

        neptune = Neptune(dtm)
        rasc, decl, dist = rec_to_sph(neptune.geo_equ_coords())

        rasc_s = rasc_string(rasc)
        decl_s = decl_string(decl)

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)
        self.assertAlmostEqual(dist, test_dist, places=2)

    def test_pluto_geo_equ_coords(self):

        """Test geo_equ_coords() method for Pluto."""

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(2001, 9, 11, 14, 0, 0, 0, UTC())

        # Right ascension accurate to within four seconds
        # Declination accurate to within one minute
        # Distance accurate to two decimal places

        test_rasc = ["16h 49m 55s", "16h 49m 56s", "16h 49m 57s",
                     "16h 49m 58s", "16h 49m 59s", "16h 50m 00s",
                     "16h 50m 01s"]
        test_decl = "-12d 12m"      # Rounded to nearest minute from NASA data
        test_dist = 30.5130064148932

        pluto = Pluto(dtm)
        rasc, decl, dist = rec_to_sph(pluto.geo_equ_coords())

        rasc_s = rasc_string(rasc)
        decl_s = decl_string(decl)

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)
        self.assertAlmostEqual(dist, test_dist, places=1)

    def test_sun_geo_equ_coords(self):

        """Test geo_equ_coords() method for the Sun."""

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(2013, 6, 4, 1, 15, 0, 0, UTC())

        # Right ascension accurate to within four seconds
        # Declination accurate to within one minute
        # Distance accurate to two decimal places

        test_rasc = ["04h 48m 00s", "04h 48m 01s", "04h 48m 02s",
                     "04h 48m 03s", "04h 48m 04s", "04h 48m 05s",
                     "04h 48m 06s"]
        test_decl = "+22d 24m"
        test_dist = 1.01446426179846

        sun = Sun(dtm)
        rasc_s = sun.right_ascension(formatted=True)
        decl_s = sun.declination(formatted=True)
        dist = sun.distance()

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)
        self.assertAlmostEqual(dist, test_dist, places=1)

    def test_moon_equ_geo_coords(self):

        """Test geo_equ_coords() method for the Moon."""

        # Test numbers taken from:
        #   http://www.stargazing.net/kepler/moon2.html

        dtm = datetime(1998, 8, 10, 0, 0, 0, 0, UTC())

        # Right ascension accurate to within two seconds
        # Declination accurate to within one minute
        # Distance not assessed in this test

        test_rasc = ["22h 56m 47s", "22h 56m 48s", "22h 56m 49s"]
        test_decl = "-07d 47m"

        moon = Moon(dtm)
        rasc_s = moon.right_ascension(formatted=True)
        decl_s = moon.declination(formatted=True)

        self.assertTrue(rasc_s in test_rasc)
        self.assertEqual(decl_s[:8], test_decl)

        # Test numbers taken from:
        #   http://ssd.jpl.nasa.gov/horizons

        dtm = datetime(2013, 6, 8, 17, 45, 0, 0, UTC())

        # Right ascension accurate to within two minutes
        # Declination accurate to within two minutes
        # Distance not assessed in this test

        test_rasc = ["05h 10m", "05h 11m", "05h 12m"]
        test_decl = ["+20d 10m", "+20d 11m", "+20d 12m"]

        moon = Moon(dtm)
        rasc_s = moon.right_ascension(formatted=True)
        decl_s = moon.declination(formatted=True)

        self.assertTrue(rasc_s[:7] in test_rasc)
        self.assertTrue(decl_s[:8] in test_decl)
