# -*- coding: utf-8 -*-
from timeit import Timer
import mesas.sas as sas
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
makeplot = False
# Initializes the random number generator so we always get the same result
np.random.seed(0)
# =====================================
# Load the input data
# =====================================
# data = pd.read_csv('../data/lower_hafren.csv', index_col=0, parse_dates=[1])[0:1000]  # 9375]
data = pd.read_csv('../data/lower_hafren.csv', index_col=0, parse_dates=False)[0:1000]  # 9375]
# length of the dataset
N = len(data)
# The individual timeseries can be pulled out of the dataframe
dS_old = data['Delta S'].values
dS = data['Delta S new'].values
J = data['J'].values
Q = data['Q'].values
E = data['Q_2'].values
C_J = data['C_J'].values
C_obs = data['C_Qo_1'].values
# =========================
# Parameters needed by sas
# =========================
# The concentration of water older than the start of observations
C_old = 7.7
ST_init = np.zeros(N + 1)

# =========================
# Create the sas functions
# =========================
# Parameters for the rSAS function
E_rSAS_fun_type = 'uniform'
ST_min_E = data['ST_min_2'].values
ST_max_E = data['ST_max_2'].values
E_rSAS_fun_parameters = np.c_[ST_min_E*1., ST_max_E*1.]
rSAS_fun_E = sas.create_function(E_rSAS_fun_type, E_rSAS_fun_parameters)
Q_rSAS_fun_type = 'gamma'
ST_min = data['ST_min'].values
ST_max = data['ST_max'].values
scale = data['scale'].values
shape = data['shape'].values
Q_rSAS_fun_parameters = np.c_[ST_min, ST_max, scale, shape]
rSAS_fun_Q_gamma = sas.create_function(Q_rSAS_fun_type, Q_rSAS_fun_parameters)
rSAS_fun = [rSAS_fun_Q_gamma, rSAS_fun_E]
# outputs = sas.solve(J, [Q, E], rSAS_fun, alpha=[[1], [0]], ST_init=ST_init,
# dt=1., n_substeps=5, C_J=C_J, C_old=[C_old], verbose=False, debug=False)


def fun(): return sas.solve(J, [Q, E], rSAS_fun, alpha=[[1], [0]], ST_init=ST_init,
                            dt=1., n_substeps=5, C_J=C_J, C_old=[C_old], verbose=False, debug=False)


t = Timer(fun)
print(t.timeit(20))


if makeplot:
    C_Q1m1 = outputs['C_Q'][:, 0, 0]
    data['C_Q_gamma'] = outputs['C_Q'][:, 0, 0]
    fig = plt.figure(2)
    plt.clf()
    plt.step(data.index, data['C_Q_gamma'], 'g', ls='-', label='gamma', lw=2, where='post')
    gooddata = ~np.isnan(C_obs)
    plt.plot(data.index[gooddata], C_obs[gooddata], 'ro',
             label='observed', lw=1, markerfacecolor='auto', markeredgecolor='r')
    plt.step(data.index[gooddata], C_obs[gooddata], 'r-', label='observed', lw=1, where='post')
    plt.legend(loc=0)
    plt.ylabel('Concentration [-]')
    plt.xlabel('time')
    plt.title('Outflow concentration')
    plt.show()
