from __future__ import print_function, absolute_import
from distutils.util import get_platform
from numpy.distutils.misc_util import Configuration, get_info
from numpy.distutils.core import setup
from os.path import join
import sys

dir_ce = 'fortran_core'
src = ['sparsekit.f90', 'dlauc1.f', 'dgeqpw.f', 'dgeqpc.f', 'dtrrnk.f',
       'dlasmx.f', 'dtrqxc.f', 'dgeqpx.f', 'dtrqpx.f', 'dgeqpb.f',
       'dsort.f', 'blassm.f', 'matvec.f', 'unary.f', 'formats.f']
print_src = ['putstrmodule.F90', 'dispmodule.f90', 'fast_direct.f90']

def configuration(parent_package='', top_path=None):
    #sys.argv.extend(['config_fc', '--fcompiler=gnu95'])
    #ce_src = [join(dir_ce, x) for x in src]
    #print_ce_src = [join(dir_ce, x) for x in print_src]
    plat_specifier = ".%s-%s" % (get_platform(), sys.version[0:3])
    inc_dir = ['build/temp%s' % plat_specifier]
    config = Configuration('ce', parent_package, top_path)
    config.add_include_dirs(inc_dir)
    config.set_options(
                       ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=False,
                       )
    print("make mylib")
    config.add_library('mylib', sources=[join(dir_ce, x) for x in src])
    config.add_library('print_lib', sources=[join(dir_ce, x) for x in print_src])
    
    print("add_subpackage('fortran_core')")
    config.add_subpackage('fortran_core')
    print("!!!!!")
    return config


if __name__ == '__main__':
    print('This is the wrong setup.py to run')
