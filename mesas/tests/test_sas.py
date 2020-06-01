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

dt = 25
Q_0 = 1.0 / dt  # <-- steady-state flow rate
C_J = 1000.0
N = 5
S_0 = 5.
eps = 0.0000000001


def steady_run(N, dt, Q_0, S_0, C_J, eps1=(0, 0), eps2=(0, 0), ST_min=0.):
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(N) * Q_0# / 2
    #data_df['Q2'] = np.ones(N) * Q_0 / 2
    data_df['J'] = np.ones(N) * Q_0
    data_df['Ca'] = np.ones(N) * C_J
    data_df['Cb'] = np.zeros(N)
    data_df['Cb'].loc[0] = C_J
    sas_fun1 = Piecewise(nsegment=1, ST_min=ST_min+eps1[0], ST_max=S_0 + eps1[1])
    #sas_fun2 = Piecewise(nsegment=1, ST_min=ST_min+eps2[0], ST_max=S_0 + eps2[1])
    sas_blends = {'Q1': Fixed(sas_fun1, N=len(data_df))}
    #sas_blends = {'Q1': Fixed(sas_fun1, N=len(data_df)), 'Q2': Fixed(sas_fun2, N=len(data_df))}
    #solute_parameters = {'Ca': {'C_old': -00, 'observations': ['Q1']}}
    solute_parameters = {'Ca': {'C_old': -0, 'observations': ['Q1', 'Q2']},
                         'Cb': {'C_old': -0, 'observations': ['Q1', 'Q2']}}
    model = Model(data_df, sas_blends, solute_parameters, debug=False, verbose=True, dt=dt, n_substeps=50)
    print('running test_steady')
    model.run()
    return model


#def test_steady_vs_analytical():
#    model = steady_run(N, dt, Q_0, S_0, C_J)
#    rdf = model.result
#    T1 = S_0 / Q_0
#    CDF = 1 - np.exp(-dt * np.arange(N + 1) / T1)
#    F = 1 - CDF * T1 / dt
#    sT_analytical = S_0 * np.diff(CDF) / dt
#    pQ_analytical = np.zeros_like(sT_analytical)
#    pQ_analytical[0] = F[1] / dt
#    pQ_analytical[1:] = np.diff(F, n=2) / dt
#    sT_analytical = np.tril(np.tile(sT_analytical, [N, 1])).T
#    pQ_analytical = np.tril(np.tile(pQ_analytical, [N, 1])).T
#    print('Expected:')
#    print(dt * sT_analytical[:, -3:])
#    print('Got:')
#    print(dt * rdf['sT'][:, -3:])
#    print('Difference/expected:')
#    print(((sT_analytical - rdf['sT'][:, 1:]) / sT_analytical)[:, -3:])
#    print('Expected:')
#    print(dt * pQ_analytical[:, -3:])
#    print('Got:')
#    for iq in range(1):
#        print(f'iq = {iq}')
#        print(dt * rdf['pQ'][:, -3:, iq])
#    print('Difference/expected:')
#    for iq in range(1):
#        print(f'iq = {iq}')
#        print(dt * ((pQ_analytical - rdf['pQ'][:, :, iq]) / pQ_analytical)[:, -3:])
#    print('Water Balance:')
#    print(dt * rdf['WaterBalance'][:, -3:])
#    print('Solute Balance:')
#    for s in range(1):
#        print(dt * rdf['SoluteBalance'][:, -3:, s])
#    print('T compare')
#    print('Expected:')
#    print(dt * sT_analytical[:, -3:] / pQ_analytical[:, -3:] / Q_0 / T1)
#    print('Got:')
#    print(dt * rdf['sT'][:, -3:] / rdf['pQ'][:, -3:, 0] / Q_0 / T1)
#    #
#    assert np.nanmax(np.abs((sT_analytical - rdf['sT'][:, 1:]) / sT_analytical)) < 1.0E-3
#    for iq in range(1):
#        assert np.nanmax(np.abs((pQ_analytical - rdf['pQ'][:, :, iq]) / pQ_analytical)) < 1.0E-3
#    assert np.abs(rdf['WaterBalance']).max() < 1.0E-6
#    assert np.abs(rdf['SoluteBalance']).max() < 1.0E-6


def test_sensitivity():
    model = steady_run(N, dt, Q_0, S_0, C_J)
    model2 = steady_run(N, dt, Q_0, S_0, C_J, eps1=[0, eps])
    rdf = model.result
    rdf2 = model2.result
    SAS_lookup, _, _, _, _, _, _ = model._create_sas_lookup()
    SAS_lookup2, _, _, _, _, _, _ = model2._create_sas_lookup()
    j = 1
    dSj = SAS_lookup2[j, N - 1] - SAS_lookup[j, N - 1]
    dsTdSj = ((rdf2['sT'][:, 1:] - rdf['sT'][:, 1:]) / dSj)
    print('dsTdSj: actual')
    print(rdf['dsTdSj'][:, 1:, j])
    print('dsTdSj: expected')
    print(dsTdSj)
    print('dsTdSj: Difference/expected:')
    err = (rdf['dsTdSj'][:, 1:, j] - dsTdSj) / dsTdSj
    print(err)
    assert np.nanmax(np.abs(err)) < 1.0E-3
    for s in range(2):
        print(s)
        dmTdSj = ((rdf2['mT'][:, 1:, s] - rdf['mT'][:, 1:, s]) / dSj)
        print('dmTdSj: actual')
        print(rdf['dmTdSj'][:, 1:, j, s])
        print('dmTdSj: expected')
        print((rdf2['mT'][:, 1:, s] - rdf['mT'][:, 1:, s]) / dSj)
        print('dmTdSj: Difference/expected:')
        err = (rdf['dmTdSj'][:, 1:, j, s] - dmTdSj) / dmTdSj
        print(err)
        assert np.nanmax(np.abs(err)) < 1.0E-3


def test_Jacobian():
    ST_min=0.1
    model = steady_run(N, dt, Q_0, S_0, C_J, ST_min=ST_min)
    rdf = model.result
    modeps = [[0, 0], [0, 0]]
    rdfeps = [[0, 0], [0, 0]]
    modeps[0][0] = steady_run(N, dt, Q_0, S_0, C_J, eps1=[eps, 0], ST_min=ST_min)
    modeps[0][1] = steady_run(N, dt, Q_0, S_0, C_J, eps1=[0, eps], ST_min=ST_min)
    #modeps[1][0] = steady_run(N, dt, Q_0, S_0, C_J, eps2=[eps, 0], ST_min=0.1)
    #modeps[1][1] = steady_run(N, dt, Q_0, S_0, C_J, eps2=[0, eps], ST_min=0.1)
    for i in range(1):
        for j in range(2):
            rdfeps[i][j] = modeps[i][j].result
    Jac = model.get_jacobian(mode='endpoint', logtransform=False)
    print('Jacobian: calculated')
    print(Jac[:,:2]/C_J)
    print('Jacobian: should be')
    J = None
    for isol, sol in enumerate(model._solorder):
        if 'observations' in model.solute_parameters[sol]:
            for isolflux, solflux in enumerate(model._fluxorder):
                if solflux in model.solute_parameters[sol]['observations']:
                    J_seg = None
                    for iflux, flux in enumerate(model._comp2learn_fluxorder):
                        for ilabel, label in enumerate(model.sas_blends[flux]._comp2learn_componentorder):
                            i = iflux
                            J_seg_this = np.column_stack([((rdfeps[i][j]['C_Q'][:, isolflux, isol] - rdf['C_Q'][:, isolflux, isol]) / eps) for j in range(2)])
                            if J_seg is None:
                                J_seg = J_seg_this
                            else:
                                J_seg = np.c_[J_seg, J_seg_this]
                    if J is None:
                        J = J_seg
                    else:
                        J = np.concatenate((J, J_seg), axis=0)
    print(J/C_J)
