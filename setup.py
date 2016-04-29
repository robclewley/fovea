#!/usr/bin/env python
"""
Setup script for fovea dynamic model diagnostic tool
"""

from setuptools import setup, os, find_packages
from setuptools import Command
from codecs import open
import sys

MAJOR = 0
MINOR = 1
MICRO = 0
__version__ = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8').read()


if sys.version_info[:2] < (2, 7) or (3, 0) <= sys.version_info[0:2] < (3, 3):
    raise RuntimeError("Python version 2.7 or >= 3.3 required.")

setup(
    name="fovea",
    version=__version__,
    packages=find_packages(),
    dependency_links = ["https://github.com/robclewley/pyeuclid/tarball/master#egg=euclid-0.01"],
    install_requires=[
        "pydstool>=0.90",
        "shapely>=1.2",
        "descartes>=1.0",
        "euclid>=0.01",
        "pyyaml>=3.11",
        "structlog>=15.1",
        "tinydb>=2.0",
        "matplotlib>=1.5.0",
    ],
    author="Rob Clewley",
    author_email="rob.clewley@gmail.com",
    maintainer="Rob Clewley",
    maintainer_email="rob.clewley@gmail.com",
    description="Dynamic modeling diagnostic and visualization tools",
    long_description=read('README.md'), # + '\n\n' + read('WHATS_NEW.txt'),
    license="BSD",
    keywords="dynamical systems, bioinformatics, modeling, diagnostics, " + \
                "visualization, literate modeling",
    url="https://github.com/robclewley/fovea",
    #download_url="https://github.com/robclewley/fovea/tarball/v%s" % __version__,
    include_package_data=True,
    platforms=["any"],
    package_data={
        '': ['*.txt', '*.rst', '*.md'],
    },
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Human Machine Interfaces",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: BSD :: FreeBSD",
        "Operating System :: POSIX :: Linux",
    ],
)
