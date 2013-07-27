#!/usr/bin/python

from distutils.core import setup

long_desc='''An astronomical library. Contains classes for the major planetary bodies enabling the right ascension, declination, distance from earth, and zodiac sign to be calculated. Also contains functions for solving Kepler's equation, transforming rectangular to spherical coordinates, and calculating Julian dates.'''

setup(name='pyastro',
      version='1.1',
      description='Astronomical Library',
      long_description=long_desc,
      author='Paul Griffiths',
      author_email='mail@paulgriffiths.net',
      url='https://github.com/paulgriffiths/pyastro',
      packages=['pyastro'],
      license='GNU General Public License version 3'
     )

