=======================
SAS Specification Class
=======================

A SAS specification is an object that can be queried to get the sas function at each timestep.
This allows us to specify a time-varying sas function in a couple of different ways.

The :class:`Fixed` spec is the simplest, producing a time-invariant function

The :class:`StateSwitch` spec switches between several sas functions, depending on the state
of the system at a given timestep

The :class:`Weighted` spec produces a weighted average of several sas functions. The weights can
vary in time

mesas.sas.specs module
----------------------

.. automodule:: mesas.sas.specs
    :members:
