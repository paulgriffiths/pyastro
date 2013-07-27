"""Unit test module for normalize_degrees() function."""

# Disable the following pylint warnings:
#  - long identifier names, these are deliberate for unittests
#  - too many public methods, supposed to be lots of tests here
#
# pylint: disable=C0103
# pylint: disable=R0904


import unittest

from ..functions import _normalize_degrees


class TestSequenceFunctions(unittest.TestCase):

    """Test sequence class for normalize_degrees function."""

    def setUp(self):
        pass

    def test_normalize_degrees(self):

        """Test normalize_degrees() function."""

        self.assertEqual(_normalize_degrees(50), 50)
        self.assertEqual(_normalize_degrees(400), 40)
        self.assertEqual(_normalize_degrees(-60), 300)
        self.assertEqual(_normalize_degrees(-460), 260)
        self.assertEqual(_normalize_degrees(500), 140)
        self.assertEqual(_normalize_degrees(360), 0)
        self.assertEqual(_normalize_degrees(0), 0)
