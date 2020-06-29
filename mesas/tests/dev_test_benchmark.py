# Todo: incorporate pytest-benchmark
import numpy as np
import pandas as pd
import pytest

# 2. A 'blender' that assumes that SAS function is fixed in time
from mesas.sas.blender import Fixed, Weighted
# classes we will use to construct the SAS model
# 1. A piecewise constant SAS function
from mesas.sas.functions import Piecewise
# 3. The sas model class
from mesas.sas.model import Model

import os
TASK_ID = os.environ('SLURM_ARRAY_TASK_ID')


timeseries_length =  500
max_age = timeseries_length
dt = 0.01
Q_0 = 1.0/timeseries_length / dt  # <-- steady-state flow rate
C_J = 1000.
C_old = 2000.
S_0 = 5
S_m = 0.601/timeseries_length
eps = 0.0001/timeseries_length
n_substeps = 20
n_segment = 5
fQ=0.3
fc=0.1
jacobian = True

print(f'timeseries_length = {timeseries_length}')
print(f'dt = {dt}')
print(f'Q_0 = {Q_0}')
print(f'C_J = {C_J}')
print(f'C_old = {C_old}')
print(f'S_0 = {S_0}')
print(f'S_m = {S_m}')
print(f'eps = {eps}')
print(f'n_substeps = {n_substeps}')
print(f'n_segment = {n_segment}')
print(f'fQ = {fQ}')
print(f'fc = {fc}')
print(f'jacobian = {jacobian}')


def steady_run(timeseries_length, dt, Q_0, S_0, C_J, j=None, ST_min=0., n_substeps=n_substeps, n_segment=n_segment, max_age=max_age):
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(timeseries_length) * Q_0
    data_df['J'] = np.ones(timeseries_length) * Q_0
    data_df['Ca'] = np.ones(timeseries_length) * C_J
    ST = np.linspace(ST_min, S_0, n_segment+1)
    if j is not None:
        ST[j] = ST[j] + eps
    sas_fun1 = Piecewise(ST=ST)
    sas_blends = {'Q1': Fixed(sas_fun1, N=len(data_df))}
    solute_parameters = {'Ca': {'C_old': C_old, 'observations': ['Q1']}}
    model = Model(data_df, sas_blends, solute_parameters, debug=True, verbose=False, dt=dt, n_substeps=n_substeps, jacobian=jacobian, max_age=max_age)
    model.run()
    return model

def test_steady_uniform(benchmark):
#def test_steady_uniform():
    print('')
    print('running test_steady_uniform')

    n = np.arange(timeseries_length)
    T_0 = S_0 / Q_0
    Delta = dt / T_0
    Kappa = np.exp(-Delta)
    Eta = Kappa ** n

    sTdisc = -((Q_0 * Eta * (-1 + Kappa)) / Delta)
    pQdisc = (Q_0 * Eta * (-1 + Kappa) ** 2) / (S_0 * Delta ** 2 * Kappa)
    pQdisc[0] = (Q_0 * (-1 + Delta + n[0] * Delta + Eta[0] * Kappa)) / (S_0 * Delta ** 2)
    mQdisc = pQdisc * C_J * Q_0
    mTdisc = -((C_J * Q_0 * Eta * (-1 + Kappa)) / Delta)
    CQdisc = (C_J * (Delta + Eta * (-1 + Kappa)) - C_old *Eta * (-1 + Kappa)) / Delta

    CQdisc = np.array([[CQdisc]]).T
    sTdisc = np.tril(np.tile(sTdisc, [timeseries_length, 1])).T
    sTdisc = np.c_[np.zeros(timeseries_length), sTdisc]
    pQdisc = np.tril(np.tile(pQdisc, [timeseries_length, 1]))
    pQdisc = np.array([pQdisc]).T
    mQdisc = np.tril(np.tile(mQdisc, [timeseries_length, 1]))
    mQdisc = np.array([[mQdisc]]).T
    mTdisc = np.tril(np.tile(mTdisc, [timeseries_length, 1])).T
    mTdisc = np.c_[np.zeros(timeseries_length), mTdisc]
    mTdisc = np.array([mTdisc.T]).T

    model = benchmark(steady_run, timeseries_length, dt, Q_0, S_0, C_J, n_segment=1, max_age=max_age)
    #model = steady_run(timeseries_length, dt, Q_0, S_0, C_J, n_segment=1)
    rdf = model.result

    def printcheck(rdf, varstr, analy):
        err = (analy[:max_age, ...] - rdf[varstr][:max_age, ...]) / analy[:max_age, ...]
        try:
            assert np.nanmax(np.abs(err)) < 1.0E-4
        except AssertionError:
            print(f'{varstr} Expected:')
            print(analy[:max_age, ...].T)
            print(f'{varstr} Got:')
            print(rdf[varstr][:max_age, ...].T)
            print(f'{varstr} Difference/expected:')
            print(err[..., :].T)
            print('')
            raise

    printcheck(rdf, 'pQ', pQdisc)
    printcheck(rdf, 'sT', sTdisc)
    printcheck(rdf, 'mQ', mQdisc)
    printcheck(rdf, 'mT', mTdisc)
    printcheck(rdf, 'C_Q', CQdisc)

    try:
        assert np.abs(rdf['WaterBalance'] / Q_0).max() < 1.0E-6
    except AssertionError:
        print('Water Balance:')
        print(rdf['WaterBalance'][:, -3:] / Q_0)
        raise

    for s in range(1):
        try:
            assert np.abs(rdf['SoluteBalance'] / (Q_0 * C_J)).max() < 1.0E-6
        except AssertionError:
            print(f'Solute Balance {s}:')
            print(rdf['SoluteBalance'][:, -3:, s] / (Q_0 * C_J))
            raise


    if jacobian:

        dsTdSjdisc = -((Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta))
        dmTdSjdisc = -((C_J * Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta))
        dCQdSjdisc = ((C_J - C_old) *Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 *Delta)

        dsTdSjdisc = np.tril(np.tile(dsTdSjdisc, [timeseries_length, 1])).T
        dsTdSjdisc = np.c_[np.zeros(timeseries_length), dsTdSjdisc]
        dmTdSjdisc = np.tril(np.tile(dmTdSjdisc, [timeseries_length, 1])).T
        dmTdSjdisc = np.c_[np.zeros(timeseries_length), dmTdSjdisc]
        dmTdSjdisc = np.array([dmTdSjdisc.T]).T
        dCQdSjdisc = np.array([[dCQdSjdisc]]).T

        model2 = steady_run(timeseries_length, dt, Q_0, S_0, C_J, j=1, n_segment=1)
        rdf2 = model2.result
        SAS_lookup, _, _, _, _, _, _ = model._create_sas_lookup()
        SAS_lookup2, _, _, _, _, _, _ = model2._create_sas_lookup()
        j = 1
        dSj = SAS_lookup2[j, timeseries_length - 1] - SAS_lookup[j, timeseries_length - 1]

        def printcheck(rdfi, rdfp, varstr, ostr, analy):
            if varstr=='dCdSj':
                var = rdfi[varstr][:,1,...]
            else:
                var = rdfi[varstr][:,:,1,...]
            err = (analy[:max_age, ...] - var[:max_age, ...]) / analy[:max_age, ...]
            check = (((rdfp[ostr] - rdfi[ostr]) / dSj))
            errcheck = (analy[:max_age, ...] - check[:max_age, ...]) / analy[:max_age, ...]
            try:
                assert np.nanmax(np.abs(err)) < 1.0E-3
                assert np.nanmax(np.abs(errcheck)) < 1.0E-2
            except AssertionError:
                print(f'{varstr} Expected:')
                print(analy[:max_age, ...].T)
                print(f'{varstr} eps check:')
                print(check[:max_age, ...].T)
                print(f'{varstr} Got:')
                print(var[:max_age, ...].T)
                print(f'{varstr} Difference/expected:')
                print(err[..., -3:].T)
                print('')
                print(f'{varstr} Difference/expected check:')
                print(errcheck[..., -3:].T)
                print('')
                raise

        printcheck(rdf, rdf2, 'dsTdSj', 'sT', dsTdSjdisc)
        printcheck(rdf, rdf2, 'dmTdSj', 'mT', dmTdSjdisc)
        printcheck(rdf, rdf2, 'dCdSj', 'C_Q', dCQdSjdisc)

