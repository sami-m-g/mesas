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

# TODO: separate the functions into two
# TODO: pieces for Q and E
# TODO: run fmin
# TODO: set the partition condition

from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rsas
from scipy.optimize import fmin
from _esther_utils import lkup_init, reg
from scipy.interpolate import interp1d

'''
    This part import data and save variable names
'''
csvfile = '../data/lower_hafren.csv'

data = pd.read_csv(csvfile, index_col=0, header=0, parse_dates=True)[:2000]
J = data['J']
Q = data['Q']
ET = data['Q_2']
C_Q = data['C_Q']  # the measured concentration, observation
C_J = data['C_J']  # input concentration
C_Q = C_Q.replace(0, np.nan)  # replace the error measurements with 0
C_old = [7.7]


'''
    This part can be refined in the future for more flexable use
'''
# numflux = 2  # number of fluxes involved
# npieces = np.ones(numflux)  # number of SAS pieces for both fluxes
# npiece = max(npieces) # takes the maximum pieces
# N = len(C_Q) # The number of observations
# '''
# Find the maximum pieces from npieceQ and npieceE and use it to make the lookup table
# '''
# [lkup, P_list] = lkup_init(N, numflux, npiece)

'''
    This part makes the lookup take and the CDF partition
'''
numflux = 2  # number of fluxes involved
npieceQ = 3  # number of SAS pieces for both fluxes
npieceE = 1
N = len(C_Q)  # The number of observations
'''
Find the maximum pieces from npieceQ and npieceE and use it to make the lookup table
'''
[lkupQ, P_listQ] = lkup_init(N, 1, npieceQ)
[lkupE, P_listE] = lkup_init(N, 1, npieceE)
lkup = [lkupQ, lkupE]
plist = [P_listQ, P_listE]

# '''
#     This part defines the rSAS function, Let's start with the values given by Harman, 2015
# '''
# # find ST, initial guess
# ST = np.zeros((numflux, npiece+1))
# ST[0][npieces[0]] = 3533
# ST[1][npieces[1]] = 2267
# for i in range(npiece - 1):
#     ST[0][npieces[0] - i - 1] = np.random.uniform(0, ST[0][npieces[0] - i], 1)
#     ST[1][npieces[1] - i - 1] = np.random.uniform(0, ST[0][npieces[1] - i], 1)

'''
    This part defines the rSAS function, Let's start with the values given by Harman, 2015
'''
np.random.seed(1)
# find SQ and SE, initial guess
SQ = np.zeros(npieceQ + 1)
SQ[npieceQ] = 3533
for i in range(npieceQ - 1):
    SQ[npieceQ - i - 1] = np.random.uniform(0, SQ[npieceQ - i], 1)

SE = np.zeros(npieceE + 1)
SE[npieceE] = 2267
for i in range(npieceE - 1):
    SE[npieceE - i - 1] = np.random.uniform(0, SE[npieceE - i], 1)

# the initial guess for params
params0 = np.r_[np.log(np.diff(SQ)), np.log(np.diff(SE))]


def _run(params, lkup, plist, J, Q, ET, C_J, C_old, C_Q):
    '''
    Given ST and initialized lkup table and corresponding P lists, along with other inputs, run the model

    :param lkup: initialized lkup table, to be updated
    :param plist:
    :param J: precip
    :param Q: discharge
    :param ET: evapotranspiration
    :param C_J: concentration of precip
    :param C_old: concentration of previous water
    :return:  opt - output from the model, lkupn - updated lkup table based on given ST
    '''

    # params to SQ and SE

    SQ = np.r_[0, np.cumsum(np.exp(params[:npieceQ]))]
    SE = np.r_[0, np.cumsum(np.exp(params[npieceQ:npieceQ+npieceE]))]

    N = len(C_J)
    P_listQ = plist[0]
    P_listE = plist[1]
    # corresponding lkup table with regularization applied
    lkupQ = np.zeros((npieceQ+101, N))
    lkupE = np.zeros((npieceE+101, N))
    for i in range(N):  # for all timesteps
        lkupQ[:, i], P_listQ_reg = reg(SQ, P_listQ)
        lkupE[:, i], P_listE_reg = reg(SE, P_listE)
    lkupn = [lkupQ, lkupE]
    P_listn = [P_listQ_reg, P_listE_reg]
    P_list_union = np.union1d(P_listQ_reg, P_listE_reg)
    lkup_union = np.zeros((len(P_list_union), N, numflux))
    for i in range(N):  # for all timesteps
        for q in range(numflux):
            lkup_union[:, i, q] = interp1d(P_listn[q], lkupn[q][:, i])(P_list_union)
    # todo: need to be fixed, "P_list must be a 1-D array"
    outputs = rsas.solve(J, [Q, ET], lkup_union, P_list=P_list_union,
                         C_J=C_J, verbose=False, C_old=C_old, n_substeps=1)
    opt = outputs['C_Q'][:, 0, 0]

    return opt


def obj(params, *args):
    # todo: how could I make sure that these params are nested?
    opt = _run(params, *args)

    # define observations
    obs_index = np.logical_not(np.isnan(C_Q))  # find out where the observation has values
    opt = opt[obs_index]
    obs = C_Q[obs_index]
    err = np.sqrt(np.mean((obs - opt)**2))
    print err, params
    return err


params_opt = fmin(obj, params0, args=(lkup, plist, J, Q, ET, C_J, C_old, C_Q))
C_Q_pred = _run(params_opt, lkup, plist, J, Q, ET, C_J, C_old, C_Q)


'''
    We define the lkup table based on the given SQ and SE
'''


''' Define the objective function as minimum SSE between outputs and observations'''


# notice that these parameters should be listed in sequences


# # Initialization
# NSE = np.zeros(N-1) # NSE metric
#
# # Observation
# obs = np.array(C_Q)
# obs_index = obs > 0 # find out where the observation has values
# obs = obs[obs_index]
# obs_mean = obs.mean()


# for n in range(N-1): # since we have N of observations, we can ultimately divide N-1 pieces
#     # Step 1: define the lookup table
#     if n == 0:
#         for i in range(npiece+1):
#             lkup[i, :, 0] = S_Q * i/(npiece)
#             lkup[i, :, 1] = S_E * i/(npiece)
#     else:
#
#         for i in range(n-1):
#             j = randint(0, npiece)
#
#
#
#
#
#
#     # Step 2: stochastical model run
#     outputs = rsas.solve(J, [Q, ET], lkup, P_list= P_list, C_J=C_J, verbose=True, C_old=C_old)
#
#     # Step 3: using three methods to find model performance
#     opt = outputs['C_Q'][:,0,0][obs_index]
#
#
#     # Step 4: find the SSE
#     SSE =
#
#
#
#     # Step 5: Decide if end the loop
#
#
#     print ("This is loop {}".format(n+1))
#
#
#
# outfile = '/home/esther/PycharmProjects/rsas/rsas/lower_hafren_2000 runs_w_3_params.csv'
#
# d = {'Q1': para_a1,'Q2': para_a2,'Q3': para_a3, 'ET1': para_b1,'ET2': para_b2,'ET3': para_b3, 'C_old': para_c, 'NSE': NSE, 'MAE': MAE, 'RMSE': RMSE}
# df = pd.DataFrame(data=d)
# df.to_csv(outfile)
#
# print("done with the process!!! file saved to {}".format(outfile))
#
#
