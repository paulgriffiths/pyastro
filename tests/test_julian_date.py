"""Unit test module for julian_date() function."""

# Disable the following pylint warnings:
#  - long identifier names, these are deliberate for unittests
#  - too many public methods, supposed to be lots of tests here
#
# pylint: disable=C0103
# pylint: disable=R0904


import unittest
from datetime import datetime

from ..functions import julian_date, UTC


class TestSequenceFunctions(unittest.TestCase):

    """Test sequence class for julian_date function."""

    def setUp(self):
        pass

    def test_julian_date(self):

        """Test julian_date() function.

        Test case dates are taken from:
            http://aa.usno.navy.mil/data/docs/JulianDate.php

        """

        dt = datetime(2013, 6, 2, 0, 0, 0, 0, UTC())
        self.assertAlmostEqual(julian_date(dt), 2456445.5)

        dt = datetime(1980, 1, 1, 0, 0, 0, 0, UTC())
        self.assertAlmostEqual(julian_date(dt), 2444239.5)

        dt = datetime(1918, 11, 11, 11, 11, 11, 0, UTC())
        self.assertAlmostEqual(julian_date(dt), 2421908.9661, places=4)
