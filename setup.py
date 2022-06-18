# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:11:18 2019

@author: ciaran
"""
from pathlib import Path

import numpy
from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="mesas.sas.solve",
            sources=["mesas/sas/solve.f90"],
            include_dirs=[numpy.get_include(), "mesas/sas/cdflib90/_build/mod"],
            library_dirs=["mesas/sas/cdflib90/_build/mod", "mesas/sas/cdflib90/_build"],
            extra_f90_compile_args=["-Ofast", "-fno-stack-arrays"],
            libraries=[p.stem for p in Path("mesas/sas/cdflib90").glob("*_mod.f90")],
        )
    ]
)
