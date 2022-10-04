# mesas

StorAge Selection is a theoretical framework for modeling transport through control volumes. It is appropriate if you are interested a system that can be treated as a single control volume (or a collection of such), and wish to make minimal assumptions about the internal organization of the transport. SAS assumes that the material leaving a system is some combination of the material that entered at earlier times. This can be useful for constructing simple models of very complicated flow systems, and for inferring the emergent transport properties of a system from tracer data.

For more information see the free HydroLearn course: [Tracers and transit times in time-variable hydrologic systems: A gentle introduction to the StorAge Selection (SAS) approach](https://edx.hydrolearn.org/courses/course-v1:JHU+570.412+Sp2020)

## Installation

The current version of mesas.py can be installed using Conda with::

    conda install -c conda-forge mesas

This will install any additional dependencies at the same time.

Alternatively, the code can be obtained from GitHub: https://github.com/charman2/mesas. Note that a fortran compiler is required to build from source (but is not required to install through Conda).

Clone the git repo and open a command prompt in the mesas directory (where the setup.py file sits). Make and install with::

    python setup.py install

## Documentation

Documentation for the code is available here: https://mesas.readthedocs.io/en/latest/

## Citation

A paper describing the code and presenting validation results is in preparation

[![DOI](https://zenodo.org/badge/183813641.svg)](https://zenodo.org/badge/latestdoi/183813641)
