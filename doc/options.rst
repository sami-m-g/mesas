.. _options:

===================
Optional parameters
===================

A number of optional parameters can be set. These can be specified in the ``config.json`` file under the ``"options"`` key:

.. code-block:: json

    {
    
    "...": "...",

    "options":{
        "dt":{
            "dt": 3600,
            "n_substeps": 5
            }
        }
    }

Options can be set as keywords when the instance is created::

    my_model = Model(..., option1='value', option2='another')

or by assigning a ``dict`` of valid key-value pairs::

    my_model.options = {option1='value', option2='another'}

In the latter case, only the options specifically referenced will be changed. The others will retain their previous values (i.e. the defaults unless they have been previously changed).

Default options are

``dt``: (scalar, default=1.0)
  Timestep in appropriate units, such that ``dt`` multiplied by any of the fluxes in ``data_df`` gives the total volume of flux over the timestep

``verbose``: (bool, default=False)
  Print information about the calculation progress

``debug``: (bool, default=False)
  Print information useful for debugging. Warning: do not use.

``warning``: (bool, default=True)
  Turn off and on warnings about calculation issues

``n_substeps``: (int, default=1)
  Number of substeps in the calculation. Each timestep can be subdivided to increase the numerical accuracy of the solution and address some numerical issues, at the cost of longer run times. Note that the substep calculations are not retained in the output. The substeps are aggregated back to the full timestep first

``max_age``: (int, default=len(data_df))
  The maximum age that will be calculated. This controls the number of rows in the output matricies. Set to a smaller value than the default to reduce calculation time (at the cost of replacing calculated concentrations of old water with ``C_old``)

``sT_init``: (1-D numpy.array, default=np.zeros(len(data_df)))
  Initial distribution of age-ranked storage (in density form). Useful for starting a run using output from another model run, e.g. for spin up. If the length of this array is less than the length of ``data_df``, then ``max_age`` will be set to ``len(sT_init)``

``influx``: (str, default=``J``)
  The name of the column in ``data_df`` containing the inflow rate

