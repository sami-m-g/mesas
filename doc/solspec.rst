.. _solspec:

============================
Specifying solute properties
============================
        This is a dictionary whose keys correspond to columns in data_df. Each entry in the dictionary
        is itself a dict with the following default key:value pairs

            'mT_init': 0.   # initial age-ranked mass in the system
            'C_old': 0.   # old water concentration
            'k1': 0.   # reaction rate constant
            'C_eq': 0.   # equilibrium concentration
            'alpha': {'Q': 1., ...}   # Partitioning coefficient for flux 'Q'
            'observations': {'Q': 'obs C in Q', ...}   # column name in data_df of observations

        note that 'alpha' and 'observations' are both dictionaries with keys corresponding to the names of fluxes

