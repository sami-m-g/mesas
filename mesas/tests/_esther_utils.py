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
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def lkup_init(N, numflux = 1,npiece = 1):
    '''
    This is used for finding n piece of look-up table with the last piece
    utilizing the regularization with an exponential shape

    :param N: the number of time steps
    :param numflux: the number of fluxes
    :param npiece: the partition of the CDF (P_list)

    :return: lookup table and P_list
             lookup[:,0,0] - the corresponding x value for the cdf
             lookup[0,:,0] - timestep
             lookup[0,0,:] - the lookup table for certain flux
    '''

    P_list = np.linspace(0, 1, npiece+1, endpoint=True)

    if numflux == 1:
        rSAS_lookup = np.zeros((len(P_list), N))
    else:
        rSAS_lookup = np.zeros((len(P_list), N, numflux))

    return rSAS_lookup,P_list

def exp_fit(ST, P_list, npiece):

    '''
    This is a little bit complicated function for fitting.

    :param ST: desired interval
    :param P_list: probability range
    :param npiece: number of pieces for SAS
    :return: The ST values correspond to the probability
    '''
    # This part is used to fit the exponential
    Y1 = P_list[(npiece - 1):]
    Y2 = np.repeat(1, 1000, axis=0)
    Y = np.concatenate((Y1, Y2), axis=None)
    # Given ST max = 0.95
    p = 95
    X = np.linspace(ST[npiece-1], ((ST[npiece] - ST[npiece-1])/ p) * len(Y), len(Y), endpoint=False)
    X = X/1000

    def func(x, a):
        return (npiece-1)/npiece + (1/npiece)*(1 - np.exp(-a*(x - ST[npiece - 1] / 1000)))

    popt, pcov = curve_fit(func, X, Y)

    return (-1 / popt * np.log(npiece * (P_list[(npiece - 1):] - 1)) + ST[npiece - 1] / 1000) * 1000


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
    ST_reg[npiece-1:] = ST95 - (ST[-1]-ST95) * np.log((1-P_list_reg[npiece-1:])/(1-P95))

    return ST_reg, P_list_reg
