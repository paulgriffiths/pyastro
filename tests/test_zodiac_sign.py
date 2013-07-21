"""Unit test module for zodiac_sign() function."""

# Disable the following pylint warnings:
#  - long identifier names, these are deliberate for unittests
#  - too many public methods, supposed to be lots of tests here
#
# pylint: disable=C0103
# pylint: disable=R0904


import unittest

from astro import zodiac_sign, rasc_to_zodiac


class TestSequenceFunctions(unittest.TestCase):

    """Test sequence class for zodiac_sign() function."""

    def setUp(self):
        pass

    def test_zodiac_sign(self):

        """Test zodiac_sign() function."""

        self.assertEqual(zodiac_sign(0), "Aries")
        self.assertEqual(zodiac_sign(15), "Aries")
        self.assertEqual(zodiac_sign(29.9), "Aries")
        self.assertEqual(zodiac_sign(30), "Taurus")
        self.assertEqual(zodiac_sign(45), "Taurus")
        self.assertEqual(zodiac_sign(59.9), "Taurus")
        self.assertEqual(zodiac_sign(60), "Gemini")
        self.assertEqual(zodiac_sign(75), "Gemini")
        self.assertEqual(zodiac_sign(89.9), "Gemini")
        self.assertEqual(zodiac_sign(90), "Cancer")
        self.assertEqual(zodiac_sign(105), "Cancer")
        self.assertEqual(zodiac_sign(119.9), "Cancer")
        self.assertEqual(zodiac_sign(120), "Leo")
        self.assertEqual(zodiac_sign(135), "Leo")
        self.assertEqual(zodiac_sign(149.9), "Leo")
        self.assertEqual(zodiac_sign(150), "Virgo")
        self.assertEqual(zodiac_sign(165), "Virgo")
        self.assertEqual(zodiac_sign(179.9), "Virgo")
        self.assertEqual(zodiac_sign(180), "Libra")
        self.assertEqual(zodiac_sign(195), "Libra")
        self.assertEqual(zodiac_sign(209.9), "Libra")
        self.assertEqual(zodiac_sign(210), "Scorpio")
        self.assertEqual(zodiac_sign(225), "Scorpio")
        self.assertEqual(zodiac_sign(239.9), "Scorpio")
        self.assertEqual(zodiac_sign(240), "Sagittarius")
        self.assertEqual(zodiac_sign(255), "Sagittarius")
        self.assertEqual(zodiac_sign(269.9), "Sagittarius")
        self.assertEqual(zodiac_sign(270), "Capricorn")
        self.assertEqual(zodiac_sign(285), "Capricorn")
        self.assertEqual(zodiac_sign(299.9), "Capricorn")
        self.assertEqual(zodiac_sign(300), "Aquarius")
        self.assertEqual(zodiac_sign(315), "Aquarius")
        self.assertEqual(zodiac_sign(329.9), "Aquarius")
        self.assertEqual(zodiac_sign(330), "Pisces")
        self.assertEqual(zodiac_sign(345), "Pisces")
        self.assertEqual(zodiac_sign(359.9), "Pisces")

    def test_zodiac_sign_short(self):

        """Test zodiac_sign() function with short option set."""

        self.assertEqual(zodiac_sign(0, short=True), "AR")
        self.assertEqual(zodiac_sign(15, short=True), "AR")
        self.assertEqual(zodiac_sign(29.9, short=True), "AR")
        self.assertEqual(zodiac_sign(30, short=True), "TA")
        self.assertEqual(zodiac_sign(45, short=True), "TA")
        self.assertEqual(zodiac_sign(59.9, short=True), "TA")
        self.assertEqual(zodiac_sign(60, short=True), "GE")
        self.assertEqual(zodiac_sign(75, short=True), "GE")
        self.assertEqual(zodiac_sign(89.9, short=True), "GE")
        self.assertEqual(zodiac_sign(90, short=True), "CN")
        self.assertEqual(zodiac_sign(105, short=True), "CN")
        self.assertEqual(zodiac_sign(119.9, short=True), "CN")
        self.assertEqual(zodiac_sign(120, short=True), "LE")
        self.assertEqual(zodiac_sign(135, short=True), "LE")
        self.assertEqual(zodiac_sign(149.9, short=True), "LE")
        self.assertEqual(zodiac_sign(150, short=True), "VI")
        self.assertEqual(zodiac_sign(165, short=True), "VI")
        self.assertEqual(zodiac_sign(179.9, short=True), "VI")
        self.assertEqual(zodiac_sign(180, short=True), "LI")
        self.assertEqual(zodiac_sign(195, short=True), "LI")
        self.assertEqual(zodiac_sign(209.9, short=True), "LI")
        self.assertEqual(zodiac_sign(210, short=True), "SC")
        self.assertEqual(zodiac_sign(225, short=True), "SC")
        self.assertEqual(zodiac_sign(239.9, short=True), "SC")
        self.assertEqual(zodiac_sign(240, short=True), "SG")
        self.assertEqual(zodiac_sign(255, short=True), "SG")
        self.assertEqual(zodiac_sign(269.9, short=True), "SG")
        self.assertEqual(zodiac_sign(270, short=True), "CP")
        self.assertEqual(zodiac_sign(285, short=True), "CP")
        self.assertEqual(zodiac_sign(299.9, short=True), "CP")
        self.assertEqual(zodiac_sign(300, short=True), "AQ")
        self.assertEqual(zodiac_sign(315, short=True), "AQ")
        self.assertEqual(zodiac_sign(329.9, short=True), "AQ")
        self.assertEqual(zodiac_sign(330, short=True), "PI")
        self.assertEqual(zodiac_sign(345, short=True), "PI")
        self.assertEqual(zodiac_sign(359.9, short=True), "PI")

    def test_rasc_to_zodiac(self):

        """Test rasc_to_zodiac() function."""

        self.assertEqual(rasc_to_zodiac(0), "00AR00")
        self.assertEqual(rasc_to_zodiac(30), "00TA00")
        self.assertEqual(rasc_to_zodiac(65.5), "05GE30")
        self.assertEqual(rasc_to_zodiac(145.7), "25LE42")
        self.assertEqual(rasc_to_zodiac(325.2), "25AQ12")

        self.assertEqual(rasc_to_zodiac(325.2, seconds=True),
                                         "25AQ12'00")
        self.assertEqual(rasc_to_zodiac(325.21, seconds=True),
                                         "25AQ12'36")
        self.assertEqual(rasc_to_zodiac(95.48, seconds=True),
                                         "05CN28'48")
