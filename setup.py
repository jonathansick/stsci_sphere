#!/usr/bin/env python

CONTACT = "Michael Droettboom"
EMAIL = "help@stsci.edu"
VERSION = "0.1"

from distutils.core import setup, Extension
import sys

try:
    import numpy
except ImportError:
    print("numpy must be installed to build sphere.")
    print("ABORTING.")
    raise

major, minor, rest = numpy.__version__.split(".", 2)
if (int(major), int(minor)) < (1, 4):
    print("numpy version 1.4 or later must be installed to build pywcs.")
    print("ABORTING.")
    raise ImportError

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

extensions = [
    Extension('sphere.math_util',
              ['src/math_util.c'],
              include_dirs = [numpy_include],
              libraries = ['m'])
    ]

setup(
    name =           'sphere',
    version =        VERSION,
    description =    "Python based tools for spherical geometry",
    author =         CONTACT,
    author_email =   EMAIL,
    packages =       ['sphere', 'sphere.test'],
    package_dir =    {'sphere': 'lib', 'sphere.test': 'lib/test'},
    package_data =   {'sphere.test': ['data/*.fits', 'data/*.fits.gz']},
    ext_modules =    extensions
    )
