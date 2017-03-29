from __future__ import print_function, absolute_import
from distutils.util import get_platform
from numpy.distutils.misc_util import Configuration, get_info
from numpy.distutils.core import setup
from os.path import join
import sys


def configuration(parent_package='', top_path=None):
    plat_specifier = ".%s-%s" % (get_platform(), sys.version[0:3])
    config = Configuration('ce', parent_package, top_path)
    config.set_options(
                       ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=False,
                       )

    config.add_subpackage('fortran_core')
    return config


if __name__ == '__main__':
    print('This is the wrong setup.py to run')
