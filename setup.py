# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:11:18 2019

@author: ciaran
"""

from distutils import util

import numpy
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

config = {
    'description': 'MESAS - Multiscale Estimation of StorAge Selection',
    'author': 'Ciaran J. Harman',
    'url': '',
    'download_url': '',
    'author_email': 'charman1@jhu.edu',
    'version': '0.2',
    'packages': ['mesas', 'mesas.sas', 'mesas.me', 'mesas.utils'],
    'package_dir': {'mesas': 'mesas',
                    'mesas.sas': util.convert_path('mesas/sas'),
                    'mesas.me': util.convert_path('mesas/me'),
                    'mesas.utils': util.convert_path('mesas/utils')
                    },
    'scripts': [],
    'name': 'mesas',
    'ext_modules': [Extension(name='solve', sources=[util.convert_path('./mesas/sas/solve.f90')],
                              include_dirs=[numpy.get_include()],
                              extra_f90_compile_args=["-Ofast", '-fno-stack-arrays'],
                              libraries=None)],
}
#extra_f90_compile_args=["-Ofast", '-fno-stack-arrays'],
#extra_f90_compile_args=["-fbacktrace", '-fcheck=all'],


setup(**config, requires=['pandas', 'numpy', 'scipy', 'matplotlib'])
