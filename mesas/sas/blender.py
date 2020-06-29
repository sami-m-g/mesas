"""
    Module blender
    ==============

    A sas blend is an object that can be queried to get the sas function at each timestep.
    This allows us to specify a time-varying sas function in a couple of different ways.

    The :class:`Fixed` blender is the simplest, producing a time-invariant function

    The :class:`StateSwitch` blender switches between several sas functions, depending on the state
    of the system at a given timestep

    The :class:`Weighted` blender produces a weighted average of several sas functions. The weights can
    vary in time
"""
import copy
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from collections import Iterable

from mesas.sas.functions import Continuous, Piecewise


class Blender:

    def __init__(self, spec, data_df):
        """
        Initializes a sas blender
        """

        self.components = OrderedDict()
        is_fixed = len(spec) == 1
        for label, component_spec in spec.items():
            self.components[label] = Component(label, component_spec, data_df, is_fixed)
        self._componentorder = list(self.components.keys())
        self._comp2learn_componentorder = self._componentorder

    def blend(self, data_df=None):
        """
        Create functions that can be called to return the CDF and inverse CDF
        """
        if data_df is None:
            data_df = self.model_data_df
        else:
            self.model_data_df = data_df

        # Trigger the method to create interpolators
        self.N = len(data_df)
        self._interp1d_inv = list(range(self.N))
        self._interp1d = list(range(self.N))
        self.ST_lists = list(range(self.N))
        self.P_lists = list(range(self.N))
        self.ST = None
        self.P = None
        for i in range(self.N):
            self.ST_lists[i] = np.sort(np.unique(np.concatenate(
                [component.sas_fun[i].ST for label, component in self.components.items()])))
            self.P_lists[i] = np.zeros_like(self.ST_lists[i])
            for label, component in self.components.items():
                self.P_lists[i] += component.sas_fun[i](self.ST_lists[i]) * component.weights.values[i]
            ST_min = self.ST_lists[i][0]
            ST_max = self.ST_lists[i][-1]
            self._interp1d_inv[i] = interp1d(self.P_lists[i], self.ST_lists[i],
                                             fill_value=(ST_min, ST_max),
                                             kind='linear', copy=False,
                                             bounds_error=False, assume_sorted=True)
            self._interp1d[i] = interp1d(self.ST_lists[i], self.P_lists[i],
                                         fill_value=(0., 1.),
                                         kind='linear', copy=False,
                                         bounds_error=False, assume_sorted=True)
        # create matricies of ST, P padded with extra values on the left
        maxj = np.max([len(ST) for ST in self.ST_lists])
        self.ST = np.zeros((maxj, self.N))
        self.P = np.zeros_like(self.ST)
        for i in range(self.N):
            Nj = len(self.ST_lists[i])
            self.ST[maxj - Nj:maxj, i] = self.ST_lists[i]
            self.P[maxj - Nj:maxj, i] = self.P_lists[i]
            if Nj < maxj:
                self.ST[:maxj - Nj, i] = -1 - np.arange(maxj - Nj)
                self.P[:maxj - Nj, i] = np.zeros(maxj - Nj)
        self.ST = self.ST.T
        self.P = self.P.T

    def __call__(self, ST, i):
        """
        The Cumulative Distribution Function of the SAS function at a timestep

        :param ST: Age-ranked Storage values (can be a number, list or array-like)
        :param i: index of the timestep
        :return: A numpy array of corresponding probabilities
        """
        return self._interp1d[i](ST)

    def inv(self, P, i):
        """
        The inverse Cumulative Distribution Function of the SAS function at a timestep

        :param P: Probabilities (can be a number, list or array-like)
        :param i: index of the timestep
        :return: A numpy array of corresponding ST values
        """
        return self._interp1d_inv[i](P)

    def __repr__(self):
        """produce a string representation of the blender"""
        repr = ''
        for label, component in self.components.items():
            repr += component.__repr__()
        return repr

    def get_parameter_list(self):
        """
        Returns a list of the parameters in each component sas function
        """
        return np.concatenate(
            [self.components[label].sas_fun.parameter_list for label in self._comp2learn_componentorder])

    def update_from_parameter_list(self, parameter_list):
        """Modify the component sas functions using a modified parameter_list"""
        starti = 0
        for label in self._comp2learn_componentorder:
            component = self.components[label]
            nparams = len(component.sas_fun.parameter_list)
            component.sas_fun.parameter_list = parameter_list[starti: starti + nparams]
            starti += nparams
        self.blend()

    def plot(self, ax=None, **kwargs):
        """Plot the component sas_functions"""
        if ax is None:
            fig, ax = plt.subplots()
        componentlines = dict([(label, component.plot(ax=ax, **kwargs))
                               for label, component in self.components.items()])
        ax.set_xlabel('$S_T$')
        ax.set_ylabel('$\Omega(S_T)$')
        ax.legend(frameon=False)
        return componentlines

    def get_jacobian(self, *args, index=None, **kwargs):
        """
        Get the jacobian of each component sas function

        Calls the get_jacobian method of each component sas function.
        Components are concatenated columnwise and returned as a single numpy array.
        """
        if index is None:
            index = np.arange(self.N)
        cat_me = [self.components[label].sas_fun.get_jacobian(*args, index=index, **kwargs).T *
                  self.components[label].weights.values[index]
                  for label in self._comp2learn_componentorder]
        return np.concatenate(cat_me, axis=0).T


class Component:
    """
    The Component class defines simple objects that package together sas
    functions with their weights timeseries and a label. The class also defines __repr__ and plot functions
    that call the corresponding functions of the underlying sas functions.
    """

    def __init__(self, label, spec, data_df, is_fixed):
        self.label = label
        self._data_df = data_df
        self.N = len(data_df)
        self.weights = pd.Series(index=data_df.index, data=1.) if is_fixed else data_df[label]
        ST = self.expand(spec['ST']) if 'ST' in spec else _NoneList()
        P = self.expand(spec['P']) if 'P' in spec else _NoneList()
        if 'args' in spec:
            assert isinstance(spec['args'], dict)
            self._args = OrderedDict()
            for arg, value in spec.pop('args').items():
                self._args[arg] = self.expand(value)
            self._args = [dict(zip(self._args, values)) for values in zip(*self._args.values())]
        if 'scipy.stats' in spec:
            fun = spec.pop('scipy.stats')
            self._sas_funs = [Continuous(fun, func_kwargs=self._args[i], P=P[i, :], **spec)
                              for i in range(self.N)]
        else:
            self._sas_funs = [Piecewise(ST=ST[i, :], P=P[i, :], **spec) for i in range(self.N)]


    def expand(self, value):
        if isinstance(value, Iterable) and not isinstance(value, str):
            return np.concatenate([[self.expand(x)] for x in value], axis=0).T
        else:
            if isinstance(value, str):
                return self._data_df[value].values
            else:
                return float(value) * np.ones(self.N)

    @property
    def ST(self):
        return np.concatenate([[sas_fun.ST] for sas_fun in self._sas_funs]).T

    @property
    def P(self):
        return np.concatenate([[sas_fun.P] for sas_fun in self._sas_funs]).T

    @property
    def sas_fun(self):
        return [sas_fun for sas_fun in self._sas_funs]

    def __get__(self, i):
        return self._sas_funs[i]

    def __repr__(self):
        repr = ''
        repr += '    component = {}'.format(self.label) + '\n'
        repr += self.sas_fun.__repr__()
        return repr

    def plot(self, *args, **kwargs):
        return self.sas_fun.plot(*args, label=self.label, **kwargs)


class _NoneList:
    def __init__(self):
        pass
    def __getitem__(self, *args, **kwargs):
        return None



