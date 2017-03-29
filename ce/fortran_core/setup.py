# This script will build the main subpackages
from distutils.util import get_platform
import sys


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_info
    sys.argv.extend(['config_fc', '--fcompiler=gnu95'])
    config = Configuration('fortran_core', parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=False)
    plat_specifier = ".%s-%s" % (get_platform(), sys.version[0:3])
    inc_dir = ['build/temp%s' % plat_specifier]

    src = ['sparsekit.f90', 'dlauc1.f', 'dgeqpw.f', 'dgeqpc.f', 'dtrrnk.f',
           'dlasmx.f', 'dtrqxc.f', 'dgeqpx.f', 'dtrqpx.f', 'dgeqpb.f',
           'dsort.f', 'blassm.f', 'matvec.f', 'unary.f', 'formats.f']
    config.add_include_dirs(inc_dir)
    config.add_library('mylib', sources=src)
    print_src = ['putstrmodule.F90', 'dispmodule.f90', 'fast_direct.f90']
    config.add_library('print_lib', sources=print_src)
    local_src = ['prec.f90']
    config.add_extension('core', sources=local_src,
                         depends=['mylib', 'print_lib'],
                         libraries=['mylib', 'print_lib', 'lapack', 'blas'])
    return config

if __name__ == '__main__':
    print('This is the wrong setup.py to run')
