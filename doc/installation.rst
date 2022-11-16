
============
Installation
============

The current version of mesas.py can be installed using Conda with::

    conda install -c conda-forge mesas

This will install any additional dependencies at the same time.

Alternatively, the code can be obtained from GitHub: https://github.com/charman2/mesas. Note that a fortran compiler is required to build from source (but is not required to install through Conda).

Clone the git repo and open a command prompt in the mesas directory (where the setup.py file sits). Make and install with::

    python setup.py install

Build from source
-----------------

Building MESAS from source is a two step process: (1) compile the *cdflib90* Fortran
source code, and (2) install the *mesas* python package.

Compile *cdflib90*
``````````````````

To compile the *cdflib90* source code (located under ``mesas/sas/cdflib90/``) you
will need a Fortran compiler and *cmake*. Both of which can be installed
through *conda*,

.. code:: bash

  $ conda install fortran-compiler cmake -c conda-forge

To compile *cdflib90*,

* create a build directory into which library and module files will be placed
* run *cmake* to setup the build folder
* run *cmake* to compile the code

.. code:: bash

  $ mkdir mesas/sas/cdflib90/_build
  $ cmake -S mesas/sas/cdflib90 -B mesas/sas/cdflib90/_build
  $ cmake --build mesas/sas/cdflib90/_build

.. note::

  The location of the build folder (``mesas/sas/cdflib90/_build/``) is important.
  The location of this folder MUST match the value that is in *setup.py*.


Build and install *mesas*
`````````````````````````

With *cdflib90* built, *mesas* can be *pip*-installed in the usual way,

.. code:: bash

  $ pip install -e .
