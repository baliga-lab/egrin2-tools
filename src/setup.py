import os
import re
import codecs
from setuptools import setup, find_packages

NAME = 'egrin2-tools'
PACKAGES = find_packages()
DESCRIPTION = 'egrin2-tools is a tool suite to build EGRIN2 networks'
LICENSE = 'LGPL V3'
URI = 'https://github.com/baliga-lab/egrin2-tools'
AUTHOR = 'Baliga Lab, Institute for Systems Biology'
VERSION = '0.8.0'

KEYWORDS = ['class', 'egrin2']

# See trove classifiers
# https://testpypi.python.org/pypi?%3Aaction=list_classifiers

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules"
]
INSTALL_REQUIRES = ['numpy', 'scipy', 'pandas', 'biopython', 'statsmodels',
                    'plotly', 'matplotlib', 'pymongo', 'cmonkey2']

if __name__ == '__main__':
    setup(name=NAME, description=DESCRIPTION,
          license=LICENSE,
          url=URI,
          version=VERSION,
          author=AUTHOR,
          author_email='wwu@systemsbiology.net',
          maintainer=AUTHOR,
          maintainer_email='wwu@systemsbiology.net',
          keywords=KEYWORDS,
          packages=PACKAGES,
          zip_safe=False,
          classifiers=CLASSIFIERS,
          install_requires=INSTALL_REQUIRES,
          scripts=['bin/egrin2-make_ensemble', 'bin/egrin2-fimojobs'])
