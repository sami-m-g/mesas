"""

    Module Models
    =============

    Text here
"""
import copy
from collections import OrderedDict
from copy import deepcopy

import numpy as np
from f_solve import f_solve

dtype = np.float64


class Model:
    '''
    Class for building and running SAS models

    To use mesas, an instance of this class must be constructed, populated with
    parameters (held in the dicts `sas_blends` and optionally `solute_parameters`),j
    associated with a dataset (`data_df`) and run using the `run` method.
    '''

    def __init__(self, data_df, sas_blends, solute_parameters=None, components_to_learn=None, **kwargs):
        # defaults
        self._result = None
        self.jacobian = {}
        self._default_options = {
            'dt': 1,
            'verbose': False,
            'debug': False,
            'warning': True,
            'full_outputs': True,
            'n_substeps': 1,
            'max_age': None,
            'ST_init': None,
            'influx': 'J',
            'ST_smallest_segment': 1./100,
            'ST_largest_segment': np.inf,
        }
        self._options = self._default_options
        # do input Checking
        self.data_df = data_df
        self.options = kwargs
        self.sas_blends = sas_blends
        # defaults for solute transport
        self._default_parameters = {
            'CS_init': 0.,
            'C_old': 0.,
            'k1': 0.,
            'C_eq': 0.,
            'alpha': dict((flux, 1.) for flux in self._fluxorder),
            'observations': {}
        }
        self.solute_parameters = solute_parameters
        if components_to_learn is None:
            self.components_to_learn = self.get_component_labels()
        else:
            self.components_to_learn = components_to_learn

    def __repr__(self):
        """Creates a repr for the model"""
        repr = ''
        for flux, sas_blend in self.sas_blends.items():
            repr += f'flux = {flux}\n'
            repr += sas_blend.__repr__()
        return repr

    def copy_without_results(self):
        """creates a copy of the model without the results"""
        return Model(copy.deepcopy(self._data_df),
                     copy.deepcopy(self._sas_blends),
                     copy.deepcopy(self._solute_parameters),
                     copy.deepcopy(self._components_to_learn),
                     **copy.deepcopy(self._options))

    def subdivided_copy(self, flux, label, segment):
        """
        Creates a copy of the model with one segment of a sas function subdivided in two

        :param flux: name of the flux
        :param label: name of the component
        :param segment: segment (numbered from 0 for the youngest)
        :return:
        """

        # make a copy
        new_model = self.copy_without_results()

        # subdivide the component
        new_model.sas_blends[flux] = new_model.sas_blends[flux].subdivided_copy(label, segment)

        return new_model

    @property
    def result(self):
        """Results of running the sas model with the current parameters

        Returns a dict with the following keys:
            'ST' : m+1 x n+1 numpy float64 2D array
                Array of age-ranked storage for n times, m ages. (full_outputs=True only)
            'PQ' : m+1 x n+1 x q numpy float64 2D array
                List of time-varying cumulative transit time distributions for n times,
                m ages, and q fluxes. (full_outputs=True only)
            'WaterBalance' : m x n numpy float64 2D array
                Should always be within tolerances of zero, unless something is very
                wrong. (full_outputs=True only)
            'C_Q' : n x q x s float64 ndarray
                If C_J is supplied, C_Q is the timeseries of outflow concentration
            'MS' : m+1 x n+1 x s float64 ndarray
                Array of age-ranked solute mass for n times, m ages, and s solutes.
                (full_outputs=True only)
            'MQ' : m+1 x n+1 x q x s float64 ndarray
                Array of age-ranked solute mass flux for n times, m ages, q fluxes and s
                solutes. (full_outputs=True only)
            'MR' : m+1 x n+1 x s float64 ndarray
                Array of age-ranked solute reaction flux for n times, m ages, and s
                solutes. (full_outputs=True only)
            'SoluteBalance' : m x n x s float64 ndarray
                Should always be within tolerances of zero, unless something is very
                wrong. (full_outputs=True only)

        For each of the arrays in the full outputs each row represents an age, and each
        column is a timestep. For n timesteps and m ages, ST will have dimensions
        (n+1) x (m+1), with the first row representing age T = 0 and the first
        column derived from the initial condition.
        """
        if self._result is None:
            raise AttributeError(
                "results are only defined once the model is run. Use .run() method to generate results ")
        else:
            return self._result

    @result.setter
    def result(self, result):
        raise AttributeError("Model results are read-only. Use .run() method to generate results ")

    @property
    def data_df(self):
        """pandas dataframe holding the model inputs"""
        return self._data_df

    @data_df.setter
    def data_df(self, new_data_df):
        self._data_df = new_data_df
        self._timeseries_length = len(self._data_df)

    @property
    def options(self):
        """Options for running the model

        To modify, assign a dictionary of valid key-value pairs

        Default options are

            'dt': 1
            'verbose': False
            'debug': False
            'warning': True
            'full_outputs': True
            'n_substeps': 1
            'max_age': None
            'ST_init': None
            'influx': 'J'
            'ST_smallest_segment': 1./100
            'ST_largest_segment': 1000000.
        """
        return self._options

    @options.setter
    def options(self, new_options):
        invalid_options = [optkey for optkey in new_options.keys()
                           if optkey not in self._default_options.keys()]
        if any(invalid_options):
            raise KeyError("Invalid options: {}".format(invalid_options))
        self._options.update(new_options)
        if self._options['max_age'] is None:
            self._options['max_age'] = self._timeseries_length
        if self._options['ST_init'] is None:
            self._options['ST_init'] = np.zeros(self._options['max_age'] + 1)
        else:
            self._options['max_age'] = len(self._options['ST_init']) - 1
        self._max_age = self._options['max_age']

    @property
    def sas_blends(self):
        return self._sas_blends

    @sas_blends.setter
    def sas_blends(self, new_sas_blends):
        self._sas_blends = new_sas_blends
        self._numflux = len(self._sas_blends)
        self._fluxorder = list(self._sas_blends.keys())

    def set_sas_blend(self, flux, sas_blend):
        self._sas_blends[flux] = sas_blend
        self._sas_blends[flux].blend()

    def set_component(self, flux, component):
        label = component.label
        self._sas_blends[flux].components[label] = component
        self._sas_blends[flux].blend()

    def set_sas_fun(self, flux, label, sas_fun):
        self._sas_blends[flux].components[label].sas_fun = sas_fun
        self._sas_blends[flux].blend()

    def get_component_labels(self):
        component_labels = {}
        for flux in self._fluxorder:
            component_labels[flux] = list(self._sas_blends[flux].components.keys())
        return component_labels

    @property
    def components_to_learn(self):
        """
        A dictionary of lists giving the component labels that you want to train with MESAS

        The components in this dictionary can be set by assigning a dictionary. For example, this would
        include only the 'min' and 'max' components of the sas blender associated with the 'Q' flux.

            >>> model.components_to_learn = {'Q':['min', 'max']}

        The parameters associated with these components can be obtained as a 1-D array with:

            >>> seglist = model.get_segment_list()

        and modified using:

            >>> model.update_from_segment_list(seglst)

        :return:
        """
        return self._components_to_learn

    @components_to_learn.setter
    def components_to_learn(self, label_dict):
        self._components_to_learn = OrderedDict((flux,
                                                 [label for label in self._sas_blends[flux]._componentorder if
                                                  label in label_dict[flux]]
                                                 ) for flux in self._fluxorder if flux in label_dict.keys())
        self._comp2learn_fluxorder = list(self._components_to_learn.keys())
        for flux in self._comp2learn_fluxorder:
            self._sas_blends[flux]._comp2learn_componentorder = self._components_to_learn[flux]

    def trim_unused_ST(self):
        largest_observed_ST = self.result['ST'].max()
        for flux, labels in self.components_to_learn.items():
            for label in labels:
                self.sas_blends[flux].components[label].trim(largest_observed_ST)




    def get_segment_list(self):
        """
        Returns a concatenated list of segments from self.components_to_learn
        """
        return np.concatenate([self.sas_blends[flux].get_segment_list() for flux in self._comp2learn_fluxorder])

    def update_from_segment_list(self, segment_list):
        """
        Modifies the components in self.components_to_learn from a concatenated list of segments
        """
        starti = 0
        for flux in self._comp2learn_fluxorder:
            nparams = len(self.sas_blends[flux].get_segment_list())
            self.sas_blends[flux].update_from_segment_list(segment_list[starti: starti + nparams])
            starti += nparams

    @property
    def solute_parameters(self):
        """Parameters describing solute behavior

        This is a dictionary whose keys correspond to columns in data_df. Each entry in the dictionary
        is itself a dict with the following default key:value pairs

            'CS_init': 0.   # initial concentration in the system
            'C_old': 0.   # old water concentration
            'k1': 0.   # reaction rate constant
            'C_eq': 0.   # equilibrium concentration
            'alpha': {'Q': 1., ...}   # Partitioning coefficient for flux 'Q'
            'observations': {'Q': 'obs C in Q', ...}   # column name in data_df of observations

        note that 'alpha' and 'observations' are both dictionaries with keys corresponding to the names of fluxes

        }

        """
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
            self._solorder = list(self._solute_parameters.keys())
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
        warning = self.options['warning']
        full_outputs = self.options['full_outputs']
        n_substeps = self.options['n_substeps']
        max_age = self.options['max_age']

        # call the Fortran code
        fresult = f_solve(
            J, Q, SAS_lookup, P_list, ST_init, dt,
            verbose, debug, warning, full_outputs,
            CS_init, C_J, alpha, k1, C_eq, C_old,
            n_substeps, nP_list, numflux, numsol, max_age, timeseries_length, nP_total)
        ST, PQ, WaterBalance, MS, MQ, MR, C_Q, SoluteBalance = fresult

        if numsol > 0:
            self._result = {'C_Q': C_Q}
        else:
            self._result = {}
        if full_outputs:
            self._result.update({'ST': ST, 'PQ': PQ, 'WaterBalance': WaterBalance})
            if numsol > 0:
                self._result.update({'MS': MS, 'MQ': MQ, 'MR': MR, 'SoluteBalance': SoluteBalance})

    def get_jacobian(self, index=None, **kwargs):
        J = None
        if index is None:
            index = np.arange(self.N)
        self.jacobian = {}
        for isol, sol in enumerate(self._solorder):
            if 'observations' in self.solute_parameters[sol]:
                self.jacobian[sol] = {}
                for iflux, flux in enumerate(self._comp2learn_fluxorder):
                    if flux in self.solute_parameters[sol]['observations']:
                        MS = np.squeeze(self.result['MS'][:, :, isol])
                        ST = self.result['ST']
                        PQ = self.result['PQ'][:, :, iflux]
                        C_old = self.solute_parameters[sol]['C_old']
                        alpha = self.solute_parameters[sol]['alpha'][flux]
                        J_seg = self.sas_blends[flux].get_jacobian(ST, MS, C_old, alpha=alpha, index=index, **kwargs)
                        J_old = np.zeros(len(index))
                        observed_index = index[index < self._max_age]
                        unobserved_index = index[index >= self._max_age]
                        J_old[:len(observed_index)] = 1 - np.diag(PQ)[1:][observed_index]
                        J_old[len(observed_index):] = 1 - PQ[-1, unobserved_index]
                        J_old_sol = np.zeros((len(index), self._numsol))
                        J_old_sol[:, list(self._solorder).index(sol)] = J_old
                        J_sol = np.c_[J_seg, J_old_sol]
                        if J is None:
                            J = J_sol
                        else:
                            J = np.concatenate(J, J_sol, axis=0)
                        self.jacobian[sol][flux] = {}
                        self.jacobian[sol][flux]['seg'] = J_seg
                        self.jacobian[sol][flux]['C_old'] = J_old
        return J

    def get_residuals(self):
        ri = None
        for isol, sol in enumerate(self._solorder):
            if 'observations' in self.solute_parameters[sol]:
                for iflux, flux in enumerate(self._comp2learn_fluxorder):
                    if flux in self.solute_parameters[sol]['observations']:
                        obs = self.solute_parameters[sol]['observations'][flux]
                        C_obs = self.data_df[obs]
                        iflux = list(self._fluxorder).index(flux)
                        isol = list(self._solorder).index(sol)
                        this_ri = self.result['C_Q'][:, iflux, isol] - C_obs.values
                        if ri is None:
                            ri = this_ri
                        else:
                            ri = np.concatenate(ri, this_ri, axis=0)
                        self.data_df[f'residual {flux}, {sol}, {obs}'] = this_ri
        return ri

    def get_obs_index(self):
        index = None
        for isol, sol in enumerate(self._solorder):
            if 'observations' in self.solute_parameters[sol]:
                for iflux, flux in enumerate(self._comp2learn_fluxorder):
                    if flux in self.solute_parameters[sol]['observations']:
                        obs = self.solute_parameters[sol]['observations'][flux]
                        C_obs = self.data_df[obs]
                        this_index = ~np.isnan(C_obs.values)
                        if index is None:
                            index = this_index
                        else:
                            index = np.concatenate(index, this_index, axis=0)
        return np.nonzero(index)[0]

    # def get_residual_parts(self, flux, sol, trainingdata, **kwargs):
    # iflux = list(self._fluxorder).index(flux)
    # isol = list(self._solorder).index(sol)
    # C_obs = self.data_df[trainingdata]
    # MS = np.squeeze(self.result['MS'][:, :, isol])
    # ST = self.result['ST']
    # PQ = self.result['PQ'][:, :, iflux]
    # C_old = self.solute_parameters[sol]['C_old']
    # alpha = self.solute_parameters[sol]['alpha'][flux]
    # r_seg = self.sas_blends[flux].get_residual_parts(C_obs, ST, MS, alpha=alpha, **kwargs)
    # r_old = (C_old - C_obs) * (1 - np.diag(PQ)[1:])
    # return r_seg, r_old
