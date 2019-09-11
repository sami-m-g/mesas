import copy
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


class Weighted:

    def __init__(self, sas_fun_dict, weights_df):
        self.weights_df = weights_df
        self.N = len(self.weights_df)
        self.components = OrderedDict(
            (label, Component(label, sas_fun_dict[label], weights_df[label]))
            for label in sas_fun_dict.keys())
        self._blend()

    def subdivided_copy(self, label, segment):
        new_blend = copy.deepcopy(self)
        new_blend.components[label].sas_fun = self.components[label].sas_fun.subdivided_copy(segment)
        return new_blend

    def _blend(self):
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
            self._interp1d_inv[i] = interp1d(self.P[:, i], self.ST[:, i],
                                             kind='linear', copy=False,
                                             bounds_error=True, assume_sorted=True)
            self._interp1d[i] = interp1d(self.ST[:, i], self.P[:, i],
                                         kind='linear', copy=False,
                                         bounds_error=True, assume_sorted=True)

    def __call__(self, ST, i):
        return self._interp1d[i](ST)

    def inv(self, P, i):
        return self._interp1d_inv[i](P)

    def __repr__(self):
        repr = '\n'
        for label, component in self.components.items():
            repr += component.__repr__()
            repr += ''+'\n'
        return repr

    def update(self, P_list, segment_list):
        starti = 0
        assert len(P_list) == len(segment_list)
        for label, component in self.components.items():
            nparams = nP = len(component.sas_fun.P) - 3
            component.sas_fun.P[2: -1] = P_list[starti: starti + nP]
            component.sas_fun.segment_list = segment_list[starti: starti + nparams]
            starti += nP
        self._blend()

    def get_P_list(self):
        return np.concatenate([component.sas_fun.P for label, component in self.components.items()])

    def update_from_P_list(self, P_list):
        starti = 0
        for label, component in self.components.items():
            nP = len(component.sas_fun.P)-3
            component.sas_fun.P[2: -1] = P_list[starti: starti+nP]
            starti += nP
        self._blend()

    def get_segment_list(self):
        return np.concatenate([component.sas_fun.segment_list for label, component in self.components.items()])

    def update_from_segment_list(self, segment_list):
        starti = 0
        for label, component in self.components.items():
            nparams = len(component.sas_fun.segment_list)
            component.sas_fun.segment_list = segment_list[starti: starti + nparams]
            starti += nparams
        self._blend()

    def plot(self, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        componentlines = dict([(label, component.plot(ax=ax, **kwargs))
                               for label, component in self.components.items()])
        # timelines = [ax.plot(self.ST[:, i], self.P[:, i], color='0.3', alpha=0.5, zorder=5)
        #             for i in range(self.N)]
        ax.set_xlabel('$S_T$')
        ax.set_ylabel('$\Omega(S_T)$')
        ax.legend(frameon=False)
        # return componentlines, timelines
        return componentlines

    def get_jacobian(self, *args, **kwargs):
        cat_me = [component.sas_fun.get_jacobian(*args, **kwargs).T * component.weights.values
                  for label, component in self.components.items()]
        return np.concatenate(cat_me, axis=0).T

    def get_residual_parts(self, *args, **kwargs):
        cat_me = [component.sas_fun.get_residual_parts(*args, **kwargs).T * component.weights.values
                  for label, component in self.components.items()]
        return np.concatenate(cat_me, axis=0).T


class Component:
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

    def __init__(self, sas_fun, N):
        weights_df = pd.DataFrame(index=range(N))
        weights_df['Fixed'] = 1.
        Weighted.__init__(self, {'Fixed': sas_fun}, weights_df)


class StateSwitch(Weighted):

    def __init__(self, sas_fun_dict, state_ts):
        weights_df = pd.DataFrame(index=range(len(state_ts)))
        for label in sas_fun_dict:
            weights_df[label] = 1. * (state_ts == label)
        Weighted.__init__(self, sas_fun_dict, weights_df)
