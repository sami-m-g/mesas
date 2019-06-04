# Copyright [2019] [Esther Fei Xu]
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from __future__ import division
import numpy as np


def lookup_init(N, npiece=1):
    '''
    This is used for finding n piece of look-up table with the last piece
    utilizing the regularization with an exponential shape

    :param N: the number of time steps
    :param ST_max: initial guess for ST_max (list, one entry per flux)
    :param npiece: the partition of the CDF (P_list)

    :return: lookup table and P_list
             lookup[:,0,0] - the corresponding x value for the cdf
             lookup[0,:,0] - timestep
             lookup[0,0,:] - the lookup table for certain flux
    '''
    npiece = np.reshape(npiece, -1)
    numflux = len(npiece)
    lookup_list = []
    P_list = []
    params_list = []
    for q in range(numflux):
        P = np.linspace(0, 1, npiece[q]+1, endpoint=True)
        lookup = np.zeros((len(P), N))
        lookup[npiece[q], :] = 1
        for i in range(npiece[q] - 1):
            lookup[npiece[q] - i - 1, :] = np.random.uniform(
                0, lookup[npiece[q] - i, 0], 1)
        lookup_list.append(lookup)
        P_list.append(P)
        params_list.append(np.diff(lookup[:, 0]))
    return lookup_list, P_list, params_list


def reg(ST, P_list):
    '''
    This is using the Taylor expansion method

    :param ST: desired interval
    :param P_list: probability range
    :return: The ST values correspond to the probability
    '''

    npiece = len(ST)
    P95 = P_list[-2] + 0.95 * (P_list[-1] - P_list[-2])
    ST95 = ST[-2] + 0.95 * (ST[-1] - ST[-2])
    P_list2 = np.linspace(P95, 1, 101, endpoint=True)
    P_list_reg = np.concatenate((P_list[:-1], P_list2), axis=None)
    ST_reg = np.zeros_like(P_list_reg)
    ST_reg[:npiece-1] = ST[:-1]
    ST_reg[npiece-1:] = ST95 - (ST[-1]-ST95) * np.log(
        (1-P_list_reg[npiece-1:])/(1-P95))

    return ST_reg, P_list_reg
