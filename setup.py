# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:11:18 2019

@author: ciaran
"""
from distutils import util
from numpy.distutils.misc_util import Configuration
import os

cdflibnames=[ 'biomath_constants_mod', 'biomath_sort_mod', 'biomath_strings_mod', 'biomath_interface_mod', 'biomath_mathlib_mod', 'zero_finder', 'cdf_aux_mod', 'cdf_beta_mod', 'cdf_binomial_mod', 'cdf_gamma_mod', 'cdf_chisq_mod', 'cdf_f_mod', 'cdf_nc_chisq_mod', 'cdf_nc_f_mod', 'cdf_normal_mod', 'cdf_t_mod', 'cdf_nc_t_mod', 'cdf_neg_binomial_mod', 'cdf_poisson_mod']

def add_mod_dir(ext,build_dir):
    mod_dir = f'{os.path.split(build_dir)[0].replace("src", "temp")}/'
    ext.include_dirs += [mod_dir]
    return None

def configuration(parent_package='',top_path=None):
    config = Configuration('mesas',parent_package,top_path)
    config.description = 'MESAS - Multiscale Estimation of StorAge Selection'
    config.author = 'Ciaran J. Harman'
    config.url = ''
    config.download_url = ''
    config.author_email = 'charman1@jhu.edu'
    config.version = '0.2021.0221'
    config.packages = ['mesas', 'mesas.sas', 'mesas.me', 'mesas.utils']
    config.package_dir = {'mesas': 'mesas',
                   'mesas.sas': util.convert_path('mesas/sas'),
                   'mesas.me': util.convert_path('mesas/me'),
                   'mesas.utils': util.convert_path('mesas/utils')
                   }
    config.scripts = []
    config.name = 'mesas'
    for src in cdflibnames:
        config.add_library(name=src,
                           sources=[f'./mesas/sas/cdflib90/{src}.f90'],
                           extra_f90_compile_args=["-Ofast"])
    config.add_extension(name='solver',
                         sources=[add_mod_dir, util.convert_path('mesas/sas/solver.f90')],
                         extra_f90_compile_args=["-Ofast", '-fno-stack-arrays']
                         )
    config.requires = ['pandas', 'numpy', 'scipy', 'matplotlib']
    config.make_config_py(name='__config__')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)



#cdflibs =  [Extension(name=name, sources=[util.convert_path(f'./mesas/sas/cdflib90/{name}.f90')], extra_f90_compile_args=["-Ofast", '-fno-stack-arrays'], libraries=None)
            #for name in cdflibnames]

#extra_f90_compile_args=["-Ofast", '-fno-stack-arrays'],
#extra_f90_compile_args=["-fast", '-acc', '-Minfo', '-Mvect=levels:10', '-ta=tesla:cc35'],
#extra_f90_compile_args=["-fast", '-acc', '-Minfo', '-ta=multicore'],
#extra_link_args=['-acc'],
#extra_f90_compile_args=["-fbacktrace", '-fcheck=all'],


#setup(**config, requires=['pandas', 'numpy', 'scipy', 'matplotlib'])
