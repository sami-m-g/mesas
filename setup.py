# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:11:18 2019

@author: ciaran
"""

from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import numpy

config = {
    'description': 'MESAS - Multiscale Estimation of StorAge Selection',
    'author': 'Ciaran J. Harman',
    'url': '',
    'download_url': '',
    'author_email': 'charman1@jhu.edu',
    'version': '0.01',
    'install_requires': ['numpy', 'scipy'],
    'packages' : ['mesas'],
    'scripts': [],
    'name': 'mesas',
    'ext_modules': [Extension(name='f_solve', sources=['./mesas/sas/solve.f90'],
                              include_dirs=[numpy.get_include()],
                              extra_f90_compile_args = ["-fcheck=all"],
                              libraries=None),
        Extension(name='f_convolve', sources=['./mesas/sas/convolve.f90'],
                              include_dirs=[numpy.get_include()],
                              extra_f90_compile_args = ["-fcheck=all"],
                              libraries=None)],
}

setup(**config)

