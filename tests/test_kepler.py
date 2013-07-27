"""Unit test module for kepler() function."""

# Disable the following pylint warnings:
#  - long identifier names, these are deliberate for unittests
#  - too many public methods, supposed to be lots of tests here
#
# pylint: disable=C0103
# pylint: disable=R0904


import unittest
from math import radians

from ..functions import kepler


class TestSequenceFunctions(unittest.TestCase):

    """Test sequence class for kepler function."""

    def setUp(self):
        pass

    def test_kepler(self):

        """Test kepler() function.

        Test case solutions are taken from:
            http://www.jgiesen.de/kepler/kepler.html

        """

        self.assertAlmostEqual(kepler(radians(20), 0.5),
                               radians(37.40006), places=5)
        self.assertAlmostEqual(kepler(radians(27), 0.5),
                               radians(48.43418), places=5)
        self.assertAlmostEqual(kepler(radians(235), 0.2),
                               radians(226.66512), places=5)
        self.assertAlmostEqual(kepler(radians(0), 0),
                               radians(0), places=5)
        self.assertAlmostEqual(kepler(radians(360), 0),
                               radians(360), places=5)
        self.assertAlmostEqual(kepler(radians(300), 0),
                               radians(300), places=5)
        self.assertAlmostEqual(kepler(radians(45), 0.9),
                               radians(96.25884), places=5)
