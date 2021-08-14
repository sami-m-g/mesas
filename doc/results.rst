.. _results:

==================
Extracting results
==================

Once the model has been run, timeseries of solute outputs can be found in ``<Model object>.data_df``. To save these as a ``.csv`` use::

    my_model.data_df.to_csv('my_model_outputs.csv')

where ``my_model`` is a Model object.

Alternatively, raw results in the form of output arrays can be accessed through the ``<Model object>.results`` property, which returns a dict with the following keys:

    ``sT`` : m x n+1 numpy float64 2D array
        Array of instantaneous age-ranked storage for n+1 times, m ages. First column is initial condition (given by the ``sT_init`` option if provided, or all zeros otherwise)
    ``pQ`` : m x n x q numpy float64 2D array
        Array of timestep-averaged time-varying backward transit time distributions over m ages, at n times, for q fluxes.
    ``WaterBalance`` : m x n numpy float64 2D array
        Should always be within tolerances of zero, unless something is very wrong.
    ``C_Q`` : n x q x s float64 ndarray
        If input concetrations are provided (see :ref:`solspec`), this gives the timeseries of timestep-averaged outflow concentration
    ``mT`` : m x n+1 x s float64 ndarray
        Array of instantaneous age-ranked solute mass over m ages, at n times, for s solutes. First column is initial condition
    ``mQ`` : m x n x q x s float64 ndarray
        Array of timestep-averaged age-ranked solute mass flux over m ages, at n times, for q fluxes and s solutes.
    ``mR`` : m x n x s float64 ndarray
        Array of timestep-averaged age-ranked solute reaction flux over m ages, at n times, for s solutes.
    ``SoluteBalance`` : m x n x s float64 ndarray
        Should always be within tolerances of zero, unless something is very wrong.

The order that the fluxes ``q`` and solutes ``s`` appear in these arrays is given by the properties ``my_model.fluxorder`` and ``my_model.solorder``. These provide lists of the column names in the order they are given in ``results``.

