"""

    Module Models
    =============

    Text here
"""
import copy
from collections import OrderedDict
from copy import deepcopy

import numpy as np
from solve import solve

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
            'n_substeps': 1,
            'max_age': None,
            'sT_init': None,
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
            'mT_init': 0.,
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
            'n_substeps': 1
            'max_age': None
            'sT_init': None
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
        if self._options['sT_init'] is None:
            self._options['sT_init'] = np.zeros(self._options['max_age'])
        else:
            self._options['max_age'] = len(self._options['sT_init'])
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

    ## Taking this out because I don't think I want to do this ever
    #def trim_unused_ST(self):
    #    largest_observed_ST = self.result['sT'].sum()
    #    for flux, labels in self.components_to_learn.items():
    #        for label in labels:
    #            self.sas_blends[flux].components[label].trim(largest_observed_ST)

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

            'mT_init': 0.   # initial age-ranked mass in the system
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
        mT_init = np.zeros((self._max_age, self._numsol), dtype=dtype)
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
                mT_init[:, isol] = _get_array(self.solute_parameters[sol]['mT_init'], self._max_age)
                k1[:, isol] = _get_array(self.solute_parameters[sol]['k1'], self._timeseries_length)
                C_eq[:, isol] = _get_array(self.solute_parameters[sol]['C_eq'],
                                           self._timeseries_length)
                for iflux, flux in enumerate(self._fluxorder):
                    alpha[:, iflux, isol] = _get_array(
                        self.solute_parameters[sol]['alpha'][flux], self._timeseries_length)
        return C_J, mT_init, C_old, alpha, k1, C_eq

    def _create_sas_lookup(self):
        nC_list = []
        nP_list = []
        component_list = []
        for flux in self._fluxorder:
            nC_list.append(len(self._sas_blends[flux].components))
            for component in self._sas_blends[flux]._componentorder:
                component_list.append(self.sas_blends[flux].components[component])
                nP_list.append(len(component_list[-1].sas_fun.P))
        nC_total = np.sum(nC_list)
        nP_total = np.sum(nP_list)
        P_list = np.column_stack([component.P for component in component_list]).T
        SAS_lookup = np.column_stack([component.ST for component in component_list]).T
        weights = np.column_stack([component.weights for component in component_list])
        return SAS_lookup, P_list, weights, nC_list, nC_total, nP_list, nP_total

    def run(self):
        # water fluxes
        J = self.data_df[self.options['influx']].values
        Q = self.data_df[self._fluxorder].values
        sT_init = self.options['sT_init']
        timeseries_length = self._timeseries_length
        numflux = self._numflux
        #
        # SAS lookup table
        SAS_lookup, P_list, weights, nC_list, nC_total, nP_list, nP_total = self._create_sas_lookup()
        SAS_lookup = np.asfortranarray(SAS_lookup)
        P_list = np.asfortranarray(P_list)
        weights = np.asfortranarray(weights)
        #
        # Solutes
        C_J, mT_init, C_old, alpha, k1, C_eq = self._create_solute_inputs()
        numsol = self._numsol
        #
        # options
        dt = self.options['dt']
        verbose = self.options['verbose']
        debug = self.options['debug']
        warning = self.options['warning']
        n_substeps = self.options['n_substeps']
        max_age = self.options['max_age']

        # call the Fortran code
        fresult = solve(
            J, Q, SAS_lookup, P_list, weights, sT_init, dt,
            verbose, debug, warning,
            mT_init, C_J, alpha, k1, C_eq, C_old,
            n_substeps, nC_list, nP_list, numflux, numsol, max_age, timeseries_length, nC_total, nP_total)
        sT, pQ, WaterBalance, mT, mQ, mR, C_Q, dsTdSj, dmTdSj, dCdSj, SoluteBalance = fresult



        if numsol > 0:
            self._result = {'C_Q': C_Q}
        else:
            self._result = {}
        self._result.update({'sT': sT, 'pQ': pQ, 'WaterBalance': WaterBalance, 'dsTdSj':dsTdSj})
        self._result['last_T'] = np.ones(self._timeseries_length + 1, dtype=int) * (self._max_age - 1)
        if np.any(sT[:, 0] > 0):
            self._result['last_T'][0] = np.nonzero(sT_init > 0)[0].max()
            off_diag = np.arange(self._result['last_T'][0]+1, self._timeseries_length - 1)
        else:
            self._result['last_T'][0] = 0
            off_diag = np.arange(1, self._timeseries_length - 1)
        self._result['last_T'][1:len(off_diag)+1] = off_diag-1
        if numsol > 0:
            self._result.update({'mT': mT, 'mQ': mQ, 'mR': mR, 'SoluteBalance': SoluteBalance, 'dmTdSj':dmTdSj, 'dCdSj':dCdSj})

    def get_jacobian(self, mode='segment', logtransform=True):
        J = None
        self.jacobian = {}
        for isol, sol in enumerate(self._solorder):
            if 'observations' in self.solute_parameters[sol]:
                self.jacobian[sol] = {}
                iP = 0
                for isolflux, solflux in enumerate(self._fluxorder):
                    if solflux in self.solute_parameters[sol]['observations']:
                        J_seg = None
                        for iflux, flux in enumerate(self._comp2learn_fluxorder):
                            for label in self._components_to_learn[flux]:
                                nP = len(self.sas_blends[flux].components[label].sas_fun.P)
                                J_S = np.squeeze(self.result['dCdSj'][:, iP:iP+nP, isolflux, isol])
                                if mode == 'endpoint':
                                    pass
                                elif mode == 'segment':
                                    # To get the derivative with respect to the segment length, we add up the derivative w.r.t. the
                                    # endpoints that would be displaced by varying that segment
                                    A = np.triu(np.ones(nP), k=0)
                                    J_S = np.dot(A, J_S.T).T
                                    if logtransform:
                                        J_S = J_S * self.sas_blends[flux].components[label].sas_fun._segment_list
                                if J_seg is None:
                                    J_seg = J_S
                                else:
                                    J_seg = np.c_[J_seg, J_S]
                            PQ = np.sum(self.result['pQ'][:, :, iflux], axis=0)
                            J_old = 1 - PQ
                            J_old_sol = np.zeros((self._timeseries_length, self._numsol))
                            J_old_sol[:, list(self._solorder).index(sol)] = J_old.T
                            J_sol = np.c_[J_seg, J_old_sol]
                        if J is None:
                            J = J_sol
                        else:
                            J = np.concatenate((J, J_sol), axis=0)
                        self.jacobian[sol][flux] = {}
                        self.jacobian[sol][flux]['seg'] = J_seg
                        self.jacobian[sol][flux]['C_old'] = J_old
                    iP +=nP
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
                            ri = np.concatenate((ri, this_ri), axis=0)
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
                            index = np.concatenate((index, this_index), axis=0)
        return np.nonzero(index)[0]
