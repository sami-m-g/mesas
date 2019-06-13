import numpy as np
from f_solve import f_solve
from copy import deepcopy
dtype = np.float64


class Model:
    def __init__(self, data_df, sas_blends, solute_parameters=None, **kwargs):
        # defaults
        self._default_options = {
            'dt': 1,
            'verbose': False,
            'debug': False,
            'full_outputs': True,
            'n_substeps': 1,
            'max_age': None,
            'ST_init': None,
            'influx': 'J',
        }
        self._options = self._default_options
        # do input Checking
        self.data_df = data_df
        self.set_options(**kwargs)
        self.sas_blends = sas_blends
        # defaults for solute transport
        self._default_parameters = {
            'CS_init': 0.,
            'C_old': 0.,
            'k1': 0.,
            'C_eq': 0.,
            'alpha': dict((flux, 1.) for flux in self._fluxorder)
        }
        self.solute_parameters = solute_parameters

    @property
    def data_df(self):
        return self._data_df

    @data_df.setter
    def data_df(self, new_data_df):
        self.set_data_df(new_data_df)

    def set_data_df(self, new_data_df):
        self._data_df = new_data_df
        self._timeseries_length = len(self._data_df)

    @property
    def options(self):
        return self._options

    @options.setter
    def options(self, new_options):
        self._options = self._default_options
        self.set_options(**new_options)

    def set_options(self, **new_options):
        invalid_options = [optkey for optkey in new_options.keys()
                           if optkey not in self._default_options.keys()]
        if any(invalid_options):
            raise KeyError("Invalid options: {}".format(invalid_options))
        self._options.update(new_options)
        if self._options['max_age'] is None:
            self._options['max_age'] = self._timeseries_length
        if self._options['ST_init'] is None:
            self._options['ST_init'] = np.zeros(self._options['max_age']+1)
        else:
            self._options['max_age'] = len(self._options['ST_init'])-1
        self._max_age = self._options['max_age']

    @property
    def sas_blends(self):
        return self._sas_blends

    @sas_blends.setter
    def sas_blends(self, new_sas_blends):
        self._sas_blends = new_sas_blends
        self._numflux = len(self._sas_blends)
        self._fluxorder = self._sas_blends.keys()

    def set_sas_blend(self, flux, sas_blends):
        self._sas_blends[flux] = sas_blends

    @property
    def solute_parameters(self):
        return self._solute_parameters

    @solute_parameters.setter
    def solute_parameters(self, new_solute_parameters):
        # set parameters for solute transport
        # provide defaults if absent
        if new_solute_parameters is not None:
            self._solute_parameters = {}
            for sol, params in new_solute_parameters.items():
                self._solute_parameters[sol] = deepcopy(self._default_parameters)
                self.set_solute_parameters(sol, params)
            self._numsol = len(self._solute_parameters)
            self._solorder = self._solute_parameters.keys()
        else:
            self._numsol = 0

    def set_solute_parameters(self, sol, params):
        invalid_parameters = [paramkey for paramkey in params.keys()
                              if paramkey not in self._default_parameters.keys()]
        if any(invalid_parameters):
            raise KeyError("invalid parameters for {}: {}".format(sol, invalid_parameters))
        self._solute_parameters[sol].update(params)

    def _create_solute_inputs(self):
        C_J = np.zeros((self._timeseries_length, self._numsol), dtype=dtype)
        CS_init = np.zeros((self._max_age, self._numsol), dtype=dtype)
        C_old = np.zeros(self._numsol, dtype=dtype)
        k1 = np.zeros((self._timeseries_length, self._numsol), dtype=dtype)
        C_eq = np.zeros((self._timeseries_length, self._numsol), dtype=dtype)
        alpha = np.ones((self._timeseries_length, self._numflux, self._numsol), dtype=dtype)
        if self.solute_parameters is not None:
            def _get_array(param, N):
                if param in self.data_df:
                    return self.data_df[param].values
                else:
                    return param * np.ones(N)
            for isol, sol in enumerate(self._solorder):
                C_J[:, isol] = self.data_df[sol].values
                C_old[isol] = self.solute_parameters[sol]['C_old']
                CS_init[:, isol] = _get_array(self.solute_parameters[sol]['CS_init'], self._max_age)
                k1[:, isol] = _get_array(self.solute_parameters[sol]['k1'], self._timeseries_length)
                C_eq[:, isol] = _get_array(self.solute_parameters[sol]['C_eq'],
                                           self._timeseries_length)
                for iflux, flux in enumerate(self._fluxorder):
                    alpha[:, iflux, isol] = _get_array(
                        self.solute_parameters[sol]['alpha'][flux], self._timeseries_length)
        return C_J, CS_init, C_old, alpha, k1, C_eq

    def _create_sas_lookup(self):
        nP_list = np.array([len(self.sas_blends[flux].P) for flux in self._fluxorder])
        nP_total = np.sum(nP_list)
        P_list = np.concatenate([self.sas_blends[flux].P for flux in self._fluxorder], axis=0)
        SAS_lookup = np.concatenate([self.sas_blends[flux].ST for flux in self._fluxorder], axis=0)
        return SAS_lookup, P_list, nP_list, nP_total

    def run(self):
        # water fluxes
        J = self.data_df[self.options['influx']].values
        Q = self.data_df[self._fluxorder].values
        ST_init = self.options['ST_init']
        timeseries_length = self._timeseries_length
        numflux = self._numflux
        #
        # SAS lookup table
        SAS_lookup, P_list, nP_list, nP_total = self._create_sas_lookup()
        SAS_lookup = np.asfortranarray(SAS_lookup)
        P_list = np.asfortranarray(P_list)
        #
        # Solutes
        C_J, CS_init, C_old, alpha, k1, C_eq = self._create_solute_inputs()
        numsol = self._numsol
        #
        # options
        dt = self.options['dt']
        verbose = self.options['verbose']
        debug = self.options['debug']
        full_outputs = self.options['full_outputs']
        n_substeps = self.options['n_substeps']
        max_age = self.options['max_age']

        # call the Fortran code
        fresult = f_solve(
            J, Q, SAS_lookup, P_list, ST_init, dt,
            verbose, debug, full_outputs,
            CS_init, C_J, alpha, k1, C_eq, C_old,
            n_substeps,  nP_list, numflux, numsol, max_age, timeseries_length, nP_total)
        ST, PQ, WaterBalance, MS, MQ, MR, C_Q, SoluteBalance = fresult

        if numsol > 0:
            self.result = {'C_Q': C_Q}
        else:
            self.result = {}
        if full_outputs:
            self.result.update({'ST': ST, 'PQ': PQ, 'WaterBalance': WaterBalance})
            if numsol > 0:
                self.result.update({'MS': MS, 'MQ': MQ, 'MR': MR, 'SoluteBalance': SoluteBalance})
