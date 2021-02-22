# Todo: incorporate pytest-benchmark
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import pytest

# 2. A 'blender' that assumes that SAS function is fixed in time
from mesas.sas.specs import Fixed, Weighted
# classes we will use to construct the SAS model
# 1. A piecewise constant SAS function
from mesas.sas.functions import Piecewise
# 3. The sas model class
from mesas.sas.model import Model

dt = 0.001
J_max = 1.0 / dt  # <-- steady-state flow rate
C_J = 1000.
C_std = 500.
C_old = C_J
N = 500
burn = 0 #int(N*0.5)
ST_min1 = 0.
ST_max1 = 5.
ST_min21 = 3.
ST_max21 = 10.
ST_min22 = 0.
ST_max22 = 20.
eps = 0.0000001
n_substeps = 2

n_segment = 2
fQ = 0.3
fc = 0.5


def create_true(iq=None, ic=None, j=None):
    np.random.seed(1234)
    data_df = pd.DataFrame()
    data_df['J'] = (np.random.rand(N+burn)**2) * J_max
    J_mean = data_df['J'].mean()
    UH1 = np.exp(-np.arange(N+burn)/20.); UH1 = UH1/UH1.sum()
    UH2 = np.exp(-np.arange(N+burn)/5.); UH2 = UH2/UH2.sum()
    data_df['Q1'] = (signal.convolve(data_df['J']-J_mean, UH1, mode='full')[:N+burn]+J_mean) * fQ
    data_df['Q2'] = (signal.convolve(data_df['J']-J_mean, UH2, mode='full')[:N+burn]+J_mean) * (1-fQ)
    data_df['c21'] = data_df['Q2']/data_df['Q2'].max() * fc
    data_df['c22'] = 1 - data_df['c21']
    data_df['Ca'] = (np.random.rand(N+burn)-0.5) * C_std + C_J
    data_df['Cb'] = (np.random.rand(N+burn)-0.5) * C_std + C_J
    sas_fun1 = Piecewise(nsegment=n_segment,  ST_min=ST_min1, ST_max=ST_max1, auto='random')
    sas_fun21 = Piecewise(nsegment=n_segment, ST_min=ST_min21, ST_max=ST_max21, auto='random')
    sas_fun22 = Piecewise(nsegment=n_segment, ST_min=ST_min22, ST_max=ST_max22, auto='random')
    if j is not None:
        if iq == 0:
            sas_fun_pert = sas_fun1
        elif ic == 0:
            sas_fun_pert = sas_fun21
        elif ic == 1:
            sas_fun_pert = sas_fun22
        STp = sas_fun_pert.ST
        STp[j] = STp[j] + eps
        sas_fun_pert.ST = STp
    sas_blends = {'Q1': Fixed(sas_fun1),
                  'Q2': Weighted({'c21': sas_fun21, 'c22': sas_fun22})}
    solute_parameters = {
        'Ca': {
            'C_old': C_old,
            'observations': ['Q1', 'Q2']
        }
    }
    model = Model(data_df, sas_blends, solute_parameters, debug=False, verbose=True, dt=dt, n_substeps=n_substeps, jacobian=False)
    model.run()
    data_df['Ca Q1'] = model.result['C_Q'][:,0,0]
    data_df['Ca Q2'] = model.result['C_Q'][:,1,0]
    data_df = data_df[burn:]
    return model, data_df


true_model, data_df = create_true()

axP = plt.subplot2grid((1,2), (0,0))
true_model.sas_blends['Q1'].plot(ax=axP)
true_model.sas_blends['Q2'].plot(ax=axP)

axC = plt.subplot2grid((1,2), (0,1))
plt.plot(data_df['Ca'], label='Ca', alpha=0.5)
plt.plot(data_df['Ca Q1'], label='Ca Q1')
plt.plot(data_df['Ca Q2'], label='Ca Q2')

sas_fun1 = Piecewise(nsegment=1,  ST_min=0.01, ST_max=120.)
sas_fun21 = Piecewise(nsegment=1, ST_min=0.01, ST_max=120.)
sas_fun22 = Piecewise(nsegment=1, ST_min=0.01, ST_max=120.)
sas_blends = {
    'Q1': Fixed(sas_fun1),
    'Q2': Weighted({
                    'c21': sas_fun21,
                    'c22': sas_fun22
                    })}
solute_parameters = {
    'Ca': {
        'C_old': C_old,
        'observations': {'Q1':'Ca Q1', 'Q2':'Ca Q2'}
    }
}
inv_model = Model(data_df, sas_blends, solute_parameters,
                  debug=False, verbose=False,
                  dt=dt, n_substeps=n_substeps, ST_largest_segment=1000., jacobian=True)

# Run the recursive_split algorithm
from mesas.me.recursive_split import run as recursive_split
can_model = recursive_split(inv_model, n_splits=5, search_mode='scanning')

axC.plot(data_df.index, can_model.result['C_Q'][:,0,0], '--', label='Ca Q1 pred')
axC.plot(data_df.index, can_model.result['C_Q'][:,1,0], '--', label='Ca Q2 pred')

can_model.sas_blends['Q1'].plot(ax=axP, ls='--')
can_model.sas_blends['Q2'].plot(ax=axP, ls='--')
