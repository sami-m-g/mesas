# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:11:18 2019

@author: ciaran
"""
from pathlib import Path

import setuptools
from numpy.distutils.core import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="mesas.sas.solve",
            sources=["mesas/sas/solve.f90"],
            extra_f90_compile_args=["-Ofast", "-fPIC", "-fno-stack-arrays"],
        )
    ]
)
