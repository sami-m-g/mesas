# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:11:18 2019

@author: ciaran
"""
from pathlib import Path

import setuptools
from numpy.distutils.core import Extension, setup

CDFLIB_SRC_DIR = Path("mesas/sas/cdflib90")
CDFLIB_BUILD_DIR = CDFLIB_SRC_DIR / "_build"
CDFLIB_MODULE_DIR = CDFLIB_BUILD_DIR / "mod"
CDFLIB_MODULES = [p.stem for p in CDFLIB_SRC_DIR.glob("*_mod.f90")]


setup(
    ext_modules=[
        Extension(
            name="mesas.sas.solve",
            sources=["mesas/sas/solve.f90"],
            include_dirs=[str(CDFLIB_MODULE_DIR)],
            library_dirs=[str(CDFLIB_MODULE_DIR), str(CDFLIB_BUILD_DIR)],
            libraries=CDFLIB_MODULES
        )
    ]
)

            #extra_f90_compile_args=["-Ofast", "-fno-stack-arrays", "-pg","-O"], extra_link_args=["-pg", "-O"],f2py_options=['-pg','-O']