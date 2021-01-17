.. _results:

=======
Results
=======
        Returns a dict with the following keys:

            'sT' : m x n+1 numpy float64 2D array
                Array of instantaneous age-ranked storage for n times, m ages. First column is initial condition
            'pQ' : m x n x q numpy float64 2D array
                Array of timestep-averaged time-varying cumulative transit time distributions for n times,
                m ages, and q fluxes.
            'WaterBalance' : m x n numpy float64 2D array
                Should always be within tolerances of zero, unless something is very
                wrong.
            'C_Q' : n x q x s float64 ndarray
                If C_J is supplied, C_Q is the timeseries of timestep-averaged outflow concentration
            'mT' : m x n+1 x s float64 ndarray
                Array of instantaneous age-ranked solute mass for n times, m ages, and s solutes. First column is initial condition
            'mQ' : m x n x q x s float64 ndarray
                Array of timestep-averaged age-ranked solute mass flux for n times, m ages, q fluxes and s
                solutes.
            'mR' : m x n x s float64 ndarray
                Array of timestep-averaged age-ranked solute reaction flux for n times, m ages, and s
                solutes.
            'SoluteBalance' : m x n x s float64 ndarray
                Should always be within tolerances of zero, unless something is very
                wrong.

        For each of the arrays in the full outputs each row represents an age, and each
        column is a timestep. For n timesteps and m ages, sT will have dimensions
        (m) x (n+1), with the first row representing age T = dt and the first
        column derived from the initial condition.
