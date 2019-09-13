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


class Weighted:
    """
    Blender class for averaging several sas functions with time-varying weights

    The sas_functions are input in a dict, and the weights are input in a pandas dataframe.
    Each dict key must match a column of the dataframe.

    Example
    -------

    Here we create two timeseries of weights, one that rises in time, and one that falls. The rising one
    has a uniform distribution over the youngest 10 units of storage. The falling one has a uniform
    distribution over a much larger volume, 200 units. Note that when we evaluate it at a couple of ST values
    [5., 150.] the returned probabilities are different at timestep 2 and timestep 8

        >>>import pandas as pd
        >>>from mesas.sas.functions import Piecewise
        >>>weights_df=pd.DataFrame(data={'rising':np.linspace(0,1,10), 'falling':np.linspace(1,0,10)})
        >>>sas_fun_dict={'rising':Piecewise(ST=[0,10]), 'falling':Piecewise(ST=[0,200])}
        >>>sas_blend=Weighted(sas_fun_dict=sas_fun_dict, weights_df=weights_df)
        >>>sas_blend([5., 150.], 2)
        array([0.13055556, 0.80555556])
        >>>sas_blend([5., 150.], 8)
        array([0.44722222, 0.97222222])
    """

    def __init__(self, sas_fun_dict, weights_df):
        """
        Initializes a Weighted sas blender

        :param sas_fun_dict: A dict of sas functions
        :param weights_df: Timeseries of weights. Must contain columns with names matching the keys in `sas_fun_dict`
        """

        self.weights_df = weights_df
        self.N = len(self.weights_df)

        # The sas_function objects are stored in an ordered dictionary of :class:`Component` objects
        self.components = OrderedDict(
            (label, Component(label, sas_fun_dict[label], weights_df[label]))
            for label in sas_fun_dict.keys())

        # Trigger the method to create interpolators
        self._blend()

    def subdivided_copy(self, label, segment, **kwargs):
        """
        Returns a copy of itself with one segment in one componet sas function divided into two

        :param label: Which component sas function to subdivide
        :param segment: Which segment to subdivide
        :param kwargs: keyword arguments passed through to the function call
        :return: a new blender object
        """
        new_blend = copy.deepcopy(self)
        new_blend.components[label].sas_fun = self.components[label].sas_fun.subdivided_copy(segment, **kwargs)
        return new_blend

    def _blend(self):
        """
        Create functions that can be called to return the CDF and inverse CDF
        """
        _ST = np.sort(np.unique(np.concatenate(
            [component.sas_fun.ST for label, component in self.components.items()])))
        self.ST = np.broadcast_to(_ST, (self.N, len(_ST))).T
        self.P = np.zeros_like(self.ST)
        for label, component in self.components.items():
            _P = component.sas_fun(self.ST)
            self.P[:, :] += (_P * np.broadcast_to(component.weights, (len(_ST), self.N)))
        self._interp1d_inv = list(range(self.N))
        self._interp1d = list(range(self.N))
        for i in range(self.N):
            ST_min = self.ST[0, i]
            ST_max = self.ST[-1, i]
            self._interp1d_inv[i] = interp1d(self.P[:, i], self.ST[:, i],
                                             fill_value=(ST_min, ST_max),
                                             kind='linear', copy=False,
                                             bounds_error=False, assume_sorted=True)
            self._interp1d[i] = interp1d(self.ST[:, i], self.P[:, i],
                                         fill_value=(0., 1.),
                                         kind='linear', copy=False,
                                         bounds_error=False, assume_sorted=True)

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
        repr = '\n'
        for label, component in self.components.items():
            repr += component.__repr__()
            repr += ''+'\n'
        return repr

    def get_P_list(self):
        """
        Returns a list of the probabilities at each segment end in each component sas function
        """
        return np.concatenate([component.sas_fun.P for label, component in self.components.items()])

    def get_segment_list(self):
        """
        Returns a list of the segment lengths in each component sas function
        """
        return np.concatenate([component.sas_fun.segment_list for label, component in self.components.items()])

    def update_from_P_list(self, P_list):
        """Modify the component sas functions using a modified P_list"""
        starti = 0
        for label, component in self.components.items():
            nP = len(component.sas_fun.P)
            component.sas_fun.P = P_list[starti: starti + nP]
            starti += nP
        self._blend()

    def update_from_segment_list(self, segment_list):
        """Modify the component sas functions using a modified segment_list"""
        starti = 0
        for label, component in self.components.items():
            nparams = len(component.sas_fun.segment_list)
            component.sas_fun.segment_list = segment_list[starti: starti + nparams]
            starti += nparams
        self._blend()

    def update(self, P_list, segment_list):
        """
        Change the probabilities and segments of all component sas functions

        :param P_list:
        :param segment_list:
        :return:
        """
        starti = 0
        assert len(P_list) == len(segment_list)
        for label, component in self.components.items():
            nparams = nP = len(component.sas_fun.P)
            component.sas_fun.P = P_list[starti: starti + nP]
            component.sas_fun.segment_list = segment_list[starti: starti + nparams]
            starti += nP
        self._blend()

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

    def get_jacobian(self, *args, **kwargs):
        """
        Get the jacobian of each component sas function

        Calls the get_jacobian method of each component sas function.
        Components are concatenated columnwise and returned as a single numpy array.
        """
        cat_me = [component.sas_fun.get_jacobian(*args, **kwargs).T * component.weights.values
                  for label, component in self.components.items()]
        return np.concatenate(cat_me, axis=0).T

    def get_residual_parts(self, *args, **kwargs):
        """Call the get_residual_parts function on each component sas function and concatenate the results"""
        cat_me = [component.sas_fun.get_residual_parts(*args, **kwargs).T * component.weights.values
                  for label, component in self.components.items()]
        return np.concatenate(cat_me, axis=0).T


class Component:
    """
    The Component class defines simple objects that package together sas
    functions with their weights timeseries and a label. The class also defines __repr__ and plot functions
    that call the corresponding functions of the underlying sas functions.
    """
    def __init__(self, label, sas_fun, weights):
        self.label = label
        self.sas_fun = sas_fun
        self.weights = weights
        self.N = len(weights)

    def __repr__(self):
        repr = ''
        repr += 'component = {}'.format(self.label)+'\n'
        repr += ''+'\n'
        repr += self.sas_fun.__repr__()
        return repr

    def plot(self, *args, **kwargs):
        return self.sas_fun.plot(*args, label=self.label, **kwargs)


class Fixed(Weighted):
    """
    SAS blender class for time-invariant sas functions.

    This class subclasses the Weighted class, wrapping its __init__ method. Fixed.__init__ just takes a
    sas function and an integer N for the number of timesteps.

    Example
    -------

        >>>from mesas.sas.functions import Piecewise
        >>>sas_fun = Piecewise(ST=[0,20])
        >>>sas_blend = Fixed(sas_fun=sas_fun, N=10)
        >>>sas_blend([5., 150.], 2)
        array([0.25, 1. ])
        >>>sas_blend([5., 150.], 8)
        array([0.25, 1. ])
    """

    def __init__(self, sas_fun, N):
        """
        creates instance of a Fixed blender

        :param sas_fun: a sas function
        :param N: number of timesteps
        """
        weights_df = pd.DataFrame(index=range(N))
        weights_df['Fixed'] = 1.
        Weighted.__init__(self, {'Fixed': sas_fun}, weights_df)


class StateSwitch(Weighted):
    """
    Blender class for switching between several sas functions through time

    The sas_functions are input in a dict, and the state is input as an array-like timeseries
    (list, numpy array, pandas dataframe).
    Each entry in the timeseries must match a dict key.

    Example
    -------

    Here we create a random timeseries of strings 'a' and 'b', and two sas functions labeled 'a' and 'b'.

        >>>from mesas.sas.functions import Piecewise
        >>>state_ts = np.where(np.random.rand(10)>0.5, 'a', 'b')
        >>>sas_fun_dict = {'a':Piecewise(ST=[0,10]), 'b':Piecewise(ST=[0,200])}
        >>>sas_blend = StateSwitch(sas_fun_dict=sas_fun_dict, state_ts=state_ts)
        >>>sas_blend([5., 150.], 2)
        sas_blend([5., 150.], 2)
        >>>sas_blend([5., 150.], 8)
        array([0.025, 0.75 ])
    """

    def __init__(self, sas_fun_dict, state_ts):
        """
        creates instance of a StateSwitch blender

        :param sas_fun_dict: a dict of sas functions with keys corresponding to entries in state_ts
        :param state_ts: a timeseries whose entries are each a key in sas_fun_dict
        """
        weights_df = pd.DataFrame(index=range(len(state_ts)))
        for label in sas_fun_dict:
            weights_df[label] = 1. * (state_ts == label)
        Weighted.__init__(self, sas_fun_dict, weights_df)
