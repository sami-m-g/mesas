from test_sas import steady_run
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from mesas.sas.blender import Fixed
from mesas.sas.functions import Piecewise
from mesas.sas.model import Model
import mesas.utils.vis as vis

timeseries_length =  50
max_age = timeseries_length
dt = 1/timeseries_length
Q_0 = 3.0/timeseries_length / dt  # <-- steady-state flow rate
C_J = 1000.
C_old = 0
S_0 = 1
S_m = 0.1
eps = 0.0001/timeseries_length
n_substeps = 20
n_segment = 1
fQ=0.3
fc=0.1
jacobian = False

def steady_run(timeseries_length, dt, Q_0, S_0, C_J, j=None, ST_min=0., debug=True, verbose=False, n_substeps=n_substeps, jacobian=jacobian, max_age=max_age, C_old=C_old, n_segment=n_segment):
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(timeseries_length) * Q_0
    data_df['J'] = np.ones(timeseries_length) * Q_0
    data_df['Ca'] = np.ones(timeseries_length) * C_J
    data_df['Cb'] = np.ones(timeseries_length) * 0
    data_df['Cb'].iloc[int(0.05 * timeseries_length):int(0.1 * timeseries_length)] = C_J
    ST = np.linspace(ST_min, S_0, n_segment+1)
    if j is not None:
        ST[j] = ST[j] + eps
    sas_fun1 = Piecewise(ST=ST)
    sas_blends = {'Q1': Fixed(sas_fun1, N=len(data_df))}
    solute_parameters = {'Ca': {'C_old': C_old, 'observations': ['Q1']}, 'Cb': {'C_old': C_old, 'observations': ['Q1']}}
    model = Model(data_df, sas_blends, solute_parameters, debug=True, verbose=False, dt=dt, n_substeps=n_substeps, jacobian=jacobian, max_age=max_age)
    model.run()
    return model


model = steady_run(timeseries_length, dt, Q_0, S_0, C_J, j=None, ST_min=S_m, debug=True, verbose=False, n_substeps=n_substeps, jacobian=jacobian, max_age=max_age, C_old=C_old)

flux = 'Q1'
sol = 'Cb'
i = 40

from collections import OrderedDict
artists = OrderedDict()
vis.plot_SAS_cumulative(model, flux, i=i, artists_dict=artists)

vis.plot_transport_column_with_timeseries(model, flux, sol, i=i, nST=20, cmap='cividis_r', TC_frac=0.3, vrange=[0, C_J], ST_max=S_0)

ani = vis.make_transport_column_animation(model, flux, sol, nST=20, cmap='cividis_r', TC_frac=0.3, vrange=[0, C_J], ST_max=S_0)
ani.save(f'test_make_transport_column_animation.gif', writer='imagemagick')

