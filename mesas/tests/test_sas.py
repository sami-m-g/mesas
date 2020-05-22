# Todo: incorporate pytest-benchmark
import numpy as np
import pandas as pd
import pytest

# 2. A 'blender' that assumes that SAS function is fixed in time
from mesas.sas.blender import Fixed
# classes we will use to construct the SAS model
# 1. A piecewise constant SAS function
from mesas.sas.functions import Piecewise
# 3. The sas model class
from mesas.sas.model import Model

dt = 24.
Q_0 = 1.0 / dt  # <-- steady-state flow rate
C_J = 100.0
N = 5
S_0 = 3.
eps = 0.0000001


def steady_run(N, dt, Q_0, S_0, C_J, eps):
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(N) * Q_0 / 2
    data_df['Q2'] = np.ones(N) * Q_0 / 2
    data_df['J'] = np.ones(N) * Q_0
    data_df['Ca'] = np.ones(N) * C_J
    data_df['Cb'] = np.zeros(N)
    data_df['Cb'].loc[0] = C_J
    sas_fun1 = Piecewise(nsegment=1, ST_max=S_0 + eps)
    sas_fun2 = Piecewise(nsegment=1, ST_max=S_0)
    sas_blends = {'Q1': Fixed(sas_fun1, N=len(data_df)), 'Q2': Fixed(sas_fun2, N=len(data_df))}
    solute_parameters = {'Ca': {}, 'Cb': {}}
    model = Model(data_df, sas_blends, solute_parameters, debug=False, verbose=True, dt=dt, n_substeps=2)
    print('running test_steady')
    model.run()
    return model


def test_steady_vs_analytical():
    model = steady_run(N, dt, Q_0, S_0, C_J, 0)
    rdf = model.result
    T1 = S_0 / Q_0
    CDF = 1 - np.exp(-dt * np.arange(N + 1) / T1)
    F = 1 - CDF * T1 / dt
    sT_analytical = S_0 * np.diff(CDF) / dt
    pQ_analytical = np.zeros_like(sT_analytical)
    pQ_analytical[0] = F[1] / dt
    pQ_analytical[1:] = np.diff(F, n=2) / dt
    sT_analytical = np.tril(np.tile(sT_analytical, [N, 1])).T
    pQ_analytical = np.tril(np.tile(pQ_analytical, [N, 1])).T
    print('Expected:')
    print(sT_analytical[:, -1])
    print('Got:')
    print(rdf['sT'][:, -1])
    print('Difference/expected:')
    print(((sT_analytical - rdf['sT'][:, 1:]) / sT_analytical)[:, -1])
    print('Expected:')
    print(pQ_analytical[:, -1])
    print('Got:')
    for iq in range(1):
        print(f'iq = {iq}')
        print(rdf['pQ'][:, -1, iq])
    print('Difference/expected:')
    for iq in range(1):
        print(f'iq = {iq}')
        print(((pQ_analytical - rdf['pQ'][:, :, iq]) / pQ_analytical)[:, -1])
    print('Water Balance:')
    print(rdf['WaterBalance'][:, -1])
    print('Solute Balance:')
    for s in range(1):
        print(rdf['SoluteBalance'][:, -1, s])
    print('T compare')
    print('Expected:')
    print(sT_analytical[:, -1] / pQ_analytical[:, -1] / Q_0 / T1)
    print('Got:')
    print(rdf['sT'][:, -1] / rdf['pQ'][:, -1, 0] / Q_0 / T1)
    #
    assert np.nanmax(np.abs((sT_analytical - rdf['sT'][:, 1:]) / sT_analytical)) < 1.0E-3
    for iq in range(1):
        assert np.nanmax(np.abs((pQ_analytical - rdf['pQ'][:, :, iq]) / pQ_analytical)) < 1.0E-3
    assert np.abs(rdf['WaterBalance']).max() < 1.0E-6
    assert np.abs(rdf['SoluteBalance']).max() < 1.0E-6


def test_sensitivity():
    model = steady_run(N, dt, Q_0, S_0, C_J, 0)
    model2 = steady_run(N, dt, Q_0, S_0, C_J, eps)
    rdf = model.result
    rdf2 = model2.result
    SAS_lookup, _, _, _ = model._create_sas_lookup()
    SAS_lookup2, _, _, _ = model2._create_sas_lookup()
    j = 1
    dSj = SAS_lookup2[j, N-1] - SAS_lookup[j, N-1]
    dsTdSj = ((rdf2['sT'][:, 1:]-rdf['sT'][:, 1:])/dSj)
    print('dsTdSj: actual')
    print(rdf['dsTdSj'][:, 1:, j])
    print('dsTdSj: expected')
    print(dsTdSj)
    print('dsTdSj: Difference/expected:')
    err = (rdf['dsTdSj'][:, 1:, j] - dsTdSj) / dsTdSj
    print(err)
    assert np.nanmax(np.abs(err))<1.0E-3
    for s in range(2):
        print(s)
        dmTdSj = ((rdf2['mT'][:, 1:, s]-rdf['mT'][:, 1:, s])/dSj)
        print('dmTdSj: actual')
        print(rdf['dmTdSj'][:, 1:, j, s])
        print('dmTdSj: expected')
        print((rdf2['mT'][:, 1:, s]-rdf['mT'][:, 1:, s])/dSj)
        print('dmTdSj: Difference/expected:')
        err = (rdf['dmTdSj'][:, 1:, j, s] - dmTdSj) / dmTdSj
        print(err)
        assert np.nanmax(np.abs(err))<1.0E-3
