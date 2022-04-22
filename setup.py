# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:11:18 2019

@author: ciaran
"""

from distutils import util
import os

import numpy
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

cdflibnames=['biomath_constants_mod', 'biomath_sort_mod', 'biomath_strings_mod', 'biomath_interface_mod', 'biomath_mathlib_mod', 'zero_finder', 'cdf_aux_mod', 'cdf_beta_mod', 'cdf_binomial_mod', 'cdf_gamma_mod', 'cdf_chisq_mod', 'cdf_f_mod', 'cdf_nc_chisq_mod', 'cdf_nc_f_mod', 'cdf_normal_mod', 'cdf_t_mod', 'cdf_nc_t_mod', 'cdf_neg_binomial_mod', 'cdf_poisson_mod']

def add_mod_dir(ext,build_dir):
    build_path = build_dir.split(os.sep)
    mod_dir = os.path.join(build_path[0], build_path[1].replace("src", "temp"))
    ext.include_dirs += [mod_dir]
    return None

config = {
    'description': 'MESAS - Multiscale Estimation of StorAge Selection',
    'author': 'Ciaran J. Harman',
    'url': '',
    'download_url': '',
    'author_email': 'charman1@jhu.edu',
    'version': '0.2022.0421',
    'packages': ['mesas', 'mesas.sas', 'mesas.me', 'mesas.utils'],
    'package_dir': {'mesas': 'mesas',
                    'mesas.sas': util.convert_path('mesas/sas'),
                    'mesas.me': util.convert_path('mesas/me'),
                    'mesas.utils': util.convert_path('mesas/utils')
                    },
    'scripts': [],
    'name': 'mesas',
    'ext_modules': [Extension(name='solve', sources=[add_mod_dir, util.convert_path('mesas/sas/solve.f90')],
                              include_dirs=[numpy.get_include()],
                              extra_f90_compile_args=["-Ofast", '-fno-stack-arrays'],
                              libraries=reversed(cdflibnames))],
}

config['libraries'] = []
for src in cdflibnames:
    config['libraries'] += [(
        src,
        {'sources': [f'mesas/sas/cdflib90/{src}.f90'], 'extra_f90_compile_args': ["-Ofast"]}
    )]
#extra_f90_compile_args=["-Ofast", '-fno-stack-arrays'],
#extra_f90_compile_args=["-fast", '-acc', '-Minfo', '-Mvect=levels:10', '-ta=tesla:cc35'],
#extra_f90_compile_args=["-fast", '-acc', '-Minfo', '-ta=multicore'],
#extra_link_args=['-acc'],
#extra_f90_compile_args=["-fbacktrace", '-fcheck=all'],


setup(**config, requires=['pandas', 'numpy', 'scipy', 'matplotlib'])
