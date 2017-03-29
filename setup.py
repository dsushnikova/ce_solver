#!/usr/bin/env python
"""CE: Python (with Fortran core) implementation of the "Compress and
Eliminate" algorithm. Paper: https://arxiv.org/pdf/1603.09133v3.pdf"""

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

DOCLINES = (__doc__ or '').split('\n')

PLATFORMS = [
             'Windows',
             'Linux',
             'Solaris',
             'Mac OS-X',
             'Unix',
             ]

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT License
Programming Language :: Fortran
Programming Language :: Python
Programming Language :: Python :: 2
Programming Language :: Python :: 2.6
Programming Language :: Python :: 2.7
Programming Language :: Python :: Implementation :: CPython
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: Unix
Operating System :: MacOS
"""

MAJOR = 1
MINOR = 0
MICRO = 0
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


def configuration(parent_package='', top_path=None):
    config = Configuration(None, parent_package, top_path)
    config.set_options(
                       ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True,
                       )
    config.add_subpackage('ce')
    return config


def setup_package():
    metadata = dict(
                    name='cepy',
                    version=VERSION,
                    description=DOCLINES[0],
                    long_description='\n'.join(DOCLINES[2:]),
                    url='https://bitbucket.org/dsushnikova/fdirsparsefmm',
                    author='Daria Sushnikova, Ivan Oseledets',
                    maintainer='Daria Sushnikova',
                    platforms=PLATFORMS,
                    classifiers=[line for line
                                 in CLASSIFIERS.split('\n') if line],
                    configuration=configuration,
                    )
    setup(**metadata)


if __name__ == '__main__':
    setup_package()
