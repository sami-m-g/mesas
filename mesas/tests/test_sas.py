# Todo: incorporate pytest-benchmark
import numpy as np
import pandas as pd
import pytest

from mesas.sas.model import Model

timeseries_length = 50
max_age = timeseries_length
dt = 0.01
Q_0 = 1.0 / timeseries_length / dt  # <-- steady-state flow rate
C_J = 1000.
C_old = 2000.
S_0 = 5
S_m = 0.601 / timeseries_length
eps = 0.0001 / timeseries_length
n_substeps = 20
n_segment = 5
fQ = 0.3
fc = 0.1
jacobian = True
debug = False
verbose = False

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

sas_specs = {
    'Q1': {
        'ST': [S_m, S_0]
    }
}

sas_specs = {
    'Q1': {
        'ST': [0, 'S_0']
    }
}

from scipy.stats import gamma

sas_specs = {
    'Q1': {
        'scipy.stats': gamma,
        'args': {
            'a': 2.,
            'scale': 'S_0',
            'loc': 'S_m'
        },
        'nsegments': 100
    }
}


sas_specs = {
    'Q1': {
        'scipy.stats': gamma,
        'args': {
            'a': 2.,
            'scale': 'S_0',
            'loc': 'S_m'
        },
        'P': np.arange(0,1,100)
    }
}


sas_specs = {
    'Q1': {
        'scipy.stats': gamma,
        'args': {
            'a': 2.,
            'scale': 'S_0',
            'loc': 'S_m'
        },
        'ST': ['S_m', 1000.]
    }
}


def steady_run(timeseries_length, dt, Q_0, S_0, C_J, j=None, ST_min=0., debug=debug, verbose=verbose,
               n_substeps=n_substeps, jacobian=jacobian, max_age=max_age, C_old=C_old, n_segment=n_segment):
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(timeseries_length) * Q_0
    data_df['J'] = np.ones(timeseries_length) * Q_0
    data_df['Ca'] = np.ones(timeseries_length) * C_J
    ST = np.linspace(ST_min, S_0, n_segment + 1)
    if j is not None:
        ST[j] = ST[j] + eps
    sas_specs = {'Q1':
                     {'steady_run':
                          {'ST': ST }}}
    solute_parameters = {'Ca': {'C_old': C_old, 'observations': ['Q1']}}
    model = Model(data_df, sas_specs, solute_parameters, debug=debug, verbose=verbose, dt=dt, n_substeps=n_substeps,
                  jacobian=jacobian, max_age=max_age)
    model.run()
    return model


from scipy.stats import gamma


def steady_run_continuous(timeseries_length, dt, Q_0, S_0, C_J, a, j=None, ST_min=0., debug=debug, verbose=verbose,
                          n_substeps=n_substeps, jacobian=jacobian, max_age=max_age, C_old=C_old, n_segment=n_segment):
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(timeseries_length) * Q_0
    data_df['J'] = np.ones(timeseries_length) * Q_0
    data_df['Ca'] = np.ones(timeseries_length) * C_J
    data_df['S_0'] = np.ones(timeseries_length) * S_0
    data_df['S_m'] = np.ones(timeseries_length) * S_m
    kwargs = {'a': a, 'scale': S_0, 'loc': ST_min}
    if j is not None:
        kwargs[j] = kwargs[j] + eps
    sas_specs = {'Q1':
                     {'steady_run_continuous':
                          {'scipy.stats': gamma,
                           'args': { 'a': 2.,
                                     'scale': 'S_0',
                                     'loc': 'S_m' },
                           'nsegments': 200}}}
    solute_parameters = {'Ca': {'C_old': C_old, 'observations': ['Q1']}}
    model = Model(data_df, sas_specs, solute_parameters, debug=debug, verbose=verbose, dt=dt, n_substeps=n_substeps,
                  jacobian=jacobian, max_age=max_age)
    model.run()
    return model


def steady_run_multiple(timeseries_length, dt, Q_0, S_0, C_J, iq=None, ic=None, j=None, ST_min=0.,
                        n_substeps=n_substeps, fQ=fQ, fc=fc, n_segment=n_segment):
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(timeseries_length) * Q_0 * fQ
    data_df['Q2'] = np.ones(timeseries_length) * Q_0 * (1 - fQ)
    data_df['c21'] = np.ones(timeseries_length) * fc
    data_df['c22'] = np.ones(timeseries_length) * (1 - fc)
    data_df['J'] = np.ones(timeseries_length) * Q_0
    data_df['Ca'] = np.ones(timeseries_length) * C_J
    data_df['Cb'] = np.ones(timeseries_length) * C_J
    ST = np.linspace(ST_min, S_0, n_segment + 1)
    ST1 = ST.copy()
    ST21 = ST.copy()
    ST22 = ST.copy()
    if j is not None:
        STp = ST.copy()
        STp[j] = STp[j] + eps
        if iq == 0:
            ST1 = STp
        elif ic == 0:
            ST21 = STp
        elif ic == 1:
            ST22 = STp
    sas_spec = {'Q1':
                    {'Fixed':
                           {'ST': ST1}},
                'Q2':
                    {'c21':
                         {'ST': ST21},
                     'c22':
                         {'ST': ST22}}}
    solute_parameters = {'Ca': {'C_old': C_old, 'observations': ['Q1', 'Q2']},
                         'Cb': {'C_old': C_old, 'observations': ['Q1', 'Q2']}}
    model = Model(data_df, sas_spec, solute_parameters, debug=False, verbose=False, dt=dt, n_substeps=n_substeps,
                  jacobian=jacobian)
    model.run()
    return model


def test_steady_uniform(benchmark):
    # def test_steady_uniform():
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
    CQdisc = (C_J * (Delta + Eta * (-1 + Kappa)) - C_old * Eta * (-1 + Kappa)) / Delta

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

    # model = benchmark(steady_run, timeseries_length, dt, Q_0, S_0, C_J, n_segment=1, max_age=max_age)
    model = steady_run(timeseries_length, dt, Q_0, S_0, C_J, n_segment=1, max_age=max_age)
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
        dCQdSjdisc = ((C_J - C_old) * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta)

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
        dSj = SAS_lookup2[timeseries_length - 1, j] - SAS_lookup[timeseries_length - 1, j]

        def printcheck(rdfi, rdfp, varstr, ostr, analy):
            if varstr == 'dCdSj':
                var = rdfi[varstr][:, 1, ...]
            else:
                var = rdfi[varstr][:, :, 1, ...]
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


def test_steady_gamma():
    print('')
    print('running test_steady_gamma')

    model = steady_run_continuous(timeseries_length, dt, Q_0, S_0, C_J, a=0.5, max_age=max_age)
    rdf = model.result
    print(model.sas_blends['Q1'])


def test_steady_piston_uniform():
    print('running test_steady_piston_uniform')

    n = np.arange(timeseries_length)
    T_0 = S_0 / Q_0
    Delta = dt / T_0
    Kappa = np.exp(-Delta)
    Eta = Kappa ** n
    T_m = S_m / Q_0
    m = T_m / dt
    HeavisideTheta = lambda x: np.heaviside(x, 1)

    sTdisc = Q_0 + (Q_0 * Kappa ** (m / (-1 + m * Delta)) * (
            -((-1 + Delta + n * Delta) * Kappa ** (m / (1 - m * Delta))) + (-1 + m * Delta) * Kappa ** (
            (1 + n) / (1 - m * Delta)) + (
                    (-1 + n * Delta) * Kappa ** (m / (1 - m * Delta)) + (1 - m * Delta) * Kappa ** (
                    n / (1 - m * Delta))) * HeavisideTheta(-m + n)) * HeavisideTheta(1 - m + n)) / Delta
    dsTdSjdisc = (Q_0 * Kappa ** ((1 + n) / (1 - m * Delta)) * (
            (1 + (1 - 2 * m + n) * Delta) * Kappa ** (m / (-1 + m * Delta)) + (-1 + m * Delta) * Kappa ** (
            (1 + n) / (-1 + m * Delta)) + Kappa ** (1 / (-1 + m * Delta)) * (
                    (-1 + 2 * m * Delta - n * Delta) * Kappa ** (m / (-1 + m * Delta)) + (
                    1 - m * Delta) * Kappa ** (n / (-1 + m * Delta))) * HeavisideTheta(
        -m + n)) * HeavisideTheta(1 - m + n)) / (S_0 * Delta * (-1 + m * Delta))
    pQdisc = -((Q_0 * (
            (1 + Delta - n * Delta + (-1 + m * Delta) * Kappa ** ((1 + m - n) / (-1 + m * Delta))) * HeavisideTheta(
        ((-1 - m + n) * S_0 * Delta) / Q_0) + 2 * (
                    -1 + n * Delta + (1 - m * Delta) * Kappa ** ((m - n) / (-1 + m * Delta))) * HeavisideTheta(
        ((-m + n) * S_0 * Delta) / Q_0) + (1 - (1 + n) * Delta + (-1 + m * Delta) * Kappa ** (
            (1 - m + n) / (1 - m * Delta))) * HeavisideTheta(((1 - m + n) * S_0 * Delta) / Q_0))) / (
                       S_0 * Delta ** 2))
    pQdisc[0] = (Q_0 * (-1 + Delta + n[0] * Delta + (1 - m * Delta) * Kappa ** (
            (1 - m + n[0]) / (1 - m * Delta))) * HeavisideTheta(((1 - m + n[0]) * S_0 * Delta) / Q_0)) / (
                        S_0 * Delta ** 2)
    mQdisc = pQdisc * Q_0 * C_J
    dpQdSjdisc = -((Q_0 * Kappa ** ((1 + m - n) / (-1 + m * Delta)) * ((1 + (-1 - 2 * m + n) * Delta + (
            -1 + m * Delta) * Kappa ** ((1 + m - n) / (1 - m * Delta))) * HeavisideTheta(
        ((-1 - m + n) * S_0 * Delta) / Q_0) + 2 * ((-1 + 2 * m * Delta - n * Delta) * Kappa ** (1 / (1 - m * Delta)) + (
            1 - m * Delta) * Kappa ** ((1 + m - n) / (1 - m * Delta))) * HeavisideTheta(
        ((-m + n) * S_0 * Delta) / Q_0) - ((1 - m * Delta) * Kappa ** ((1 + m - n) / (1 - m * Delta)) + (
            -1 + (-1 + 2 * m - n) * Delta) / Kappa ** (2 / (-1 + m * Delta))) * HeavisideTheta(
        ((1 - m + n) * S_0 * Delta) / Q_0))) / (S_0 ** 2 * Delta ** 2 * (-1 + m * Delta)))
    dpQdSjdisc[0] = (Q_0 * (1 - m * Delta + (-1 + (-1 + 2 * m - n[0]) * Delta) * Kappa ** (
            (1 - m + n[0]) / (1 - m * Delta))) * HeavisideTheta(((1 - m + n[0]) * S_0 * Delta) / Q_0)) / (
                            S_0 ** 2 * Delta ** 2 * (-1 + m * Delta))
    mTdisc = (C_J * Q_0 * (Delta + Kappa ** (m / (-1 + m * Delta)) * (
            -((-1 + Delta + n * Delta) * Kappa ** (m / (1 - m * Delta))) + (-1 + m * Delta) * Kappa ** (
            (1 + n) / (1 - m * Delta)) + (
                    (-1 + n * Delta) * Kappa ** (m / (1 - m * Delta)) + (1 - m * Delta) * Kappa ** (
                    n / (1 - m * Delta))) * HeavisideTheta(-m + n)) * HeavisideTheta(1 - m + n))) / Delta
    dmTdSjdisc = (C_J * Q_0 * Kappa ** ((1 + n) / (1 - m * Delta)) * (
            (1 + (1 - 2 * m + n) * Delta) * Kappa ** (m / (-1 + m * Delta)) + (-1 + m * Delta) * Kappa ** (
            (1 + n) / (-1 + m * Delta)) + Kappa ** (1 / (-1 + m * Delta)) * (
                    (-1 + 2 * m * Delta - n * Delta) * Kappa ** (m / (-1 + m * Delta)) + (
                    1 - m * Delta) * Kappa ** (n / (-1 + m * Delta))) * HeavisideTheta(
        -m + n)) * HeavisideTheta(1 - m + n)) / (S_0 * Delta * (-1 + m * Delta))
    CQdisc = C_old - ((C_J - C_old) * Kappa ** (m / (-1 + m * Delta)) * (
            -((-1 + Delta + n * Delta) * Kappa ** (m / (1 - m * Delta))) + (-1 + m * Delta) * Kappa ** (
            (1 + n) / (1 - m * Delta)) + (
                    (-1 + n * Delta) * Kappa ** (m / (1 - m * Delta)) + (1 - m * Delta) * Kappa ** (
                    n / (1 - m * Delta))) * HeavisideTheta(-m + n)) * HeavisideTheta(1 - m + n)) / Delta
    dCQdSjdisc = ((C_J - C_old) * (
            1 - m * Delta + (-1 + (-1 + 2 * m - n) * Delta) * Kappa ** ((1 - m + n) / (1 - m * Delta)) + (
            -1 + m * Delta + (1 - 2 * m * Delta + n * Delta) * Kappa ** (
            (m - n) / (-1 + m * Delta))) * HeavisideTheta(-m + n)) * HeavisideTheta(1 - m + n)) / (
                         S_0 * Delta * (-1 + m * Delta))

    dsTdSmdisc = -((Q_0 * Kappa ** ((1 - m + n) / (1 - m * Delta)) * (
            1 - m + n + (m - n) * Kappa ** (1 / (-1 + m * Delta)) * HeavisideTheta(-m + n)) * HeavisideTheta(
        1 - m + n)) / (S_0 * (-1 + m * Delta)))
    dpQdSmdisc = -((Q_0 * Kappa ** ((1 - m + n) / (1 - m * Delta)) * (
            (1 + m - n) * Kappa ** (2 / (-1 + m * Delta)) * HeavisideTheta(-1 - m + n) - 2 * (m - n) * Kappa ** (
            1 / (-1 + m * Delta)) * HeavisideTheta(-m + n) + (-1 + m - n) * HeavisideTheta(1 - m + n))) / (
                           S_0 ** 2 * Delta * (-1 + m * Delta)))
    dpQdSmdisc[0] = -(
            ((-1 + m - n[0]) * Q_0 * Kappa ** ((1 - m + n[0]) / (1 - m * Delta)) * HeavisideTheta(1 - m + n[0])) / (
            S_0 ** 2 * Delta * (-1 + m * Delta)))
    dmTdSmdisc = -((C_J * Q_0 * Kappa ** ((1 - m + n) / (1 - m * Delta)) * (
            1 - m + n + (m - n) * Kappa ** (1 / (-1 + m * Delta)) * HeavisideTheta(-m + n)) * HeavisideTheta(
        1 - m + n)) / (S_0 * (-1 + m * Delta)))
    dCQdSmdisc = ((C_J - C_old) * Kappa ** ((1 - m + n) / (1 - m * Delta)) * (
            1 - m + n + (m - n) * Kappa ** (1 / (-1 + m * Delta)) * HeavisideTheta(-m + n)) * HeavisideTheta(
        1 - m + n)) / (S_0 * (-1 + m * Delta))

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

    dsTdSjdisc = np.tril(np.tile(dsTdSjdisc, [timeseries_length, 1])).T
    dsTdSjdisc = np.c_[np.zeros(timeseries_length), dsTdSjdisc]
    dmTdSjdisc = np.tril(np.tile(dmTdSjdisc, [timeseries_length, 1])).T
    dmTdSjdisc = np.c_[np.zeros(timeseries_length), dmTdSjdisc]
    dmTdSjdisc = np.array([dmTdSjdisc.T]).T
    dCQdSjdisc = np.array([[dCQdSjdisc]]).T

    dsTdSmdisc = np.tril(np.tile(dsTdSmdisc, [timeseries_length, 1])).T
    dsTdSmdisc = np.c_[np.zeros(timeseries_length), dsTdSmdisc]
    dmTdSmdisc = np.tril(np.tile(dmTdSmdisc, [timeseries_length, 1])).T
    dmTdSmdisc = np.c_[np.zeros(timeseries_length), dmTdSmdisc]
    dmTdSmdisc = np.array([dmTdSmdisc.T]).T
    dCQdSmdisc = np.array([[dCQdSmdisc]]).T

    model = steady_run(timeseries_length, dt, Q_0, S_0, C_J, ST_min=S_m, n_segment=1)
    rdf = model.result

    def printcheck(rdfi, varstr, analy):
        err = (analy[:max_age, ...] - rdfi[varstr][:max_age, ...]) / analy[:max_age, ...]
        try:
            assert np.nanmax(np.abs(err)) < 5.0E-2
        except AssertionError:
            print(f'{varstr} Expected:')
            print(analy[:max_age, ...].T)
            print(f'{varstr} Got:')
            print(rdfi[varstr][:max_age, ...].T)
            print(f'{varstr} Difference/expected:')
            print(err[..., -3:].T)
            print('')
            raise

    printcheck(rdf, 'sT', sTdisc)
    printcheck(rdf, 'pQ', pQdisc)
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

        def printcheck(rdfi, rdfp, varstr, ostr, analy, ip):
            if varstr == 'dCdSj':
                var = rdfi[varstr][:, ip, ...]
            else:
                var = rdfi[varstr][:, :, ip, ...]
            err = (analy[:max_age, ...] - var[:max_age, ...]) / analy[:max_age, ...]
            check = (((rdfp[ostr] - rdfi[ostr]) / dSj))
            errcheck = (analy[:max_age, ...] - check[:max_age, ...]) / analy[:max_age, ...]
            try:
                assert np.nanmax(np.abs(err)) < 2.0E-1
                assert np.nanmax(np.abs(errcheck)) < 2.0E-1
            except AssertionError:
                print(f'{varstr} j={j}  Expected:')
                print(analy.T)
                print(f'{varstr} j={j}  eps check:')
                print(check.T)
                print(f'{varstr} j={j}  Got:')
                print(var.T)
                print(f'{varstr} j={j}  Difference/expected:')
                print(err[..., -3:].T)
                print('')
                raise

        model0 = steady_run(timeseries_length, dt, Q_0, S_0, C_J, j=1, ST_min=S_m, n_segment=1)
        rdf0 = model0.result
        SAS_lookup, _, _, _, _, _, _ = model._create_sas_lookup()
        SAS_lookup0, _, _, _, _, _, _ = model0._create_sas_lookup()
        j = 1
        dSj = SAS_lookup0[timeseries_length - 1, j] - SAS_lookup[timeseries_length - 1, j]
        printcheck(rdf, rdf0, 'dsTdSj', 'sT', dsTdSjdisc, j)
        printcheck(rdf, rdf0, 'dmTdSj', 'mT', dmTdSjdisc, j)
        printcheck(rdf, rdf0, 'dCdSj', 'C_Q', dCQdSjdisc, j)

        modelm = steady_run(timeseries_length, dt, Q_0, S_0, C_J, j=0, ST_min=S_m, n_segment=1)
        rdfm = modelm.result
        SAS_lookup, _, _, _, _, _, _ = model._create_sas_lookup()
        SAS_lookupm, _, _, _, _, _, _ = modelm._create_sas_lookup()
        j = 0
        dSj = SAS_lookupm[timeseries_length - 1, j] - SAS_lookup[timeseries_length - 1, j]
        printcheck(rdf, rdfm, 'dmTdSj', 'mT', dmTdSmdisc, j)
        printcheck(rdf, rdfm, 'dsTdSj', 'sT', dsTdSmdisc, j)
        printcheck(rdf, rdfm, 'dCdSj', 'C_Q', dCQdSmdisc, j)


def test_multiple():
    print('running test_multiple')

    n = np.arange(timeseries_length)
    T_0 = S_0 / Q_0
    Delta = dt / T_0
    Kappa = np.exp(-Delta)
    Eta = Kappa ** n
    T_m = S_m / Q_0
    m = T_m / dt
    HeavisideTheta = lambda x: np.heaviside(x, 1)

    sTdisc = Q_0 + (Q_0 * Kappa ** (m / (-1 + m * Delta)) * (
            -((-1 + Delta + n * Delta) * Kappa ** (m / (1 - m * Delta))) + (-1 + m * Delta) * Kappa ** (
            (1 + n) / (1 - m * Delta)) + (
                    (-1 + n * Delta) * Kappa ** (m / (1 - m * Delta)) + (1 - m * Delta) * Kappa ** (
                    n / (1 - m * Delta))) * HeavisideTheta(-m + n)) * HeavisideTheta(1 - m + n)) / Delta
    pQdisc = -((Q_0 * (
            (1 + Delta - n * Delta + (-1 + m * Delta) * Kappa ** ((1 + m - n) / (-1 + m * Delta))) * HeavisideTheta(
        ((-1 - m + n) * S_0 * Delta) / Q_0) + 2 * (
                    -1 + n * Delta + (1 - m * Delta) * Kappa ** ((m - n) / (-1 + m * Delta))) * HeavisideTheta(
        ((-m + n) * S_0 * Delta) / Q_0) + (1 - (1 + n) * Delta + (-1 + m * Delta) * Kappa ** (
            (1 - m + n) / (1 - m * Delta))) * HeavisideTheta(((1 - m + n) * S_0 * Delta) / Q_0))) / (
                       S_0 * Delta ** 2))
    pQdisc[0] = (Q_0 * (-1 + Delta + n[0] * Delta + (1 - m * Delta) * Kappa ** (
            (1 - m + n[0]) / (1 - m * Delta))) * HeavisideTheta(((1 - m + n[0]) * S_0 * Delta) / Q_0)) / (
                        S_0 * Delta ** 2)
    mQdisc = pQdisc * Q_0 * C_J
    mTdisc = (C_J * Q_0 * (Delta + Kappa ** (m / (-1 + m * Delta)) * (
            -((-1 + Delta + n * Delta) * Kappa ** (m / (1 - m * Delta))) + (-1 + m * Delta) * Kappa ** (
            (1 + n) / (1 - m * Delta)) + (
                    (-1 + n * Delta) * Kappa ** (m / (1 - m * Delta)) + (1 - m * Delta) * Kappa ** (
                    n / (1 - m * Delta))) * HeavisideTheta(-m + n)) * HeavisideTheta(1 - m + n))) / Delta
    CQdisc = C_old - ((C_J - C_old) * Kappa ** (m / (-1 + m * Delta)) * (
            -((-1 + Delta + n * Delta) * Kappa ** (m / (1 - m * Delta))) + (-1 + m * Delta) * Kappa ** (
            (1 + n) / (1 - m * Delta)) + (
                    (-1 + n * Delta) * Kappa ** (m / (1 - m * Delta)) + (1 - m * Delta) * Kappa ** (
                    n / (1 - m * Delta))) * HeavisideTheta(-m + n)) * HeavisideTheta(1 - m + n)) / Delta

    CQdisc = np.array([[CQdisc]]).T
    sTdisc = np.tril(np.tile(sTdisc, [timeseries_length, 1])).T
    sTdisc = np.c_[np.zeros(timeseries_length), sTdisc]
    pQdisc = np.tril(np.tile(pQdisc, [timeseries_length, 1]))
    pQdisc = np.array([pQdisc]).T
    mQdisc = np.tril(np.tile(mQdisc, [timeseries_length, 1])).T
    mTdisc = np.tril(np.tile(mTdisc, [timeseries_length, 1])).T
    mTdisc = np.c_[np.zeros(timeseries_length), mTdisc]
    mTdisc = np.array([mTdisc.T]).T

    model = steady_run_multiple(timeseries_length, dt, Q_0, S_0, C_J, iq=None, ic=None, j=None, ST_min=S_m,
                                n_substeps=n_substeps, fQ=fQ, fc=fc, n_segment=n_segment)
    rdf = model.result

    def printcheck(rdfi, varstr, analy):
        if varstr == 'mQ':
            var = rdfi[varstr][:, :, 0, 0]
            analy = analy * fQ
        else:
            var = rdfi[varstr]
        err = (analy[:max_age, ...] - var[:max_age, ...]) / analy[:max_age, ...]
        try:
            assert np.nanmax(np.abs(err)) < 5.0E-2
        except AssertionError:
            print(f'{varstr} Expected:')
            print(analy[:max_age, ...].T)
            print(f'{varstr} Got:')
            print(var[:max_age, ...].T)
            print(f'{varstr} Difference/expected:')
            print(err[..., -3:].T)
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

        def printcheck(rdfi, rdfp, varstr, ostr, norm, ip):
            var = rdfi[varstr][:, :, ip, ...]
            check = (rdfp[ostr] - rdfi[ostr]) / dSj
            err = (check[:max_age, ...] - var[:max_age, ...]) / norm
            try:
                assert np.nanmax(np.abs(err)) < 5.0E-2
            except AssertionError:
                print(f'{varstr} ip={ip} eps check:')
                print(check[:max_age, ...].T)
                print(f'{varstr} ip={ip} Got:')
                print(var[:max_age, ...].T)
                print(f'{varstr} ip={ip} Difference/norm')
                print(err[..., -3:].T)
                print('')
                raise

        def printcheckC(rdfi, rdfp, varstr, ostr, ip, iq, s):
            var = rdfi[varstr][:, ip, iq, s]
            check = (rdfp[ostr][:, iq, s] - rdfi[ostr][:, iq, s]) / dSj
            err = (check - var) / C_J
            try:
                assert np.nanmax(np.abs(err)) < 5.0E-2
            except:
                print(f'{varstr} ip={ip} eps check:')
                print(check.T)
                print(f'{varstr} ip={ip} Got:')
                print(var.T)
                print(f'{varstr} ip={ip} Difference/CJ:')
                print(err[..., :].T)
                print('')
                raise

        SAS_lookup, _, _, _, _, _, _ = model._create_sas_lookup()
        for iq, ic, ip0 in [(0, 0, 0 * (n_segment + 1)), (1, 0, 1 * (n_segment + 1)), (1, 1, 2 * (n_segment + 1))]:
            for j in range(n_segment + 1):
                ip = ip0 + j
                print(iq, ic, j, ip)
                modelp = steady_run_multiple(timeseries_length, dt, Q_0, S_0, C_J, iq=iq, ic=ic, j=j, ST_min=S_m,
                                             n_substeps=n_substeps, fQ=fQ, fc=fc, n_segment=n_segment)
                rdfp = modelp.result
                SAS_lookupp, _, _, _, _, _, _ = modelp._create_sas_lookup()
                dSj = SAS_lookupp[timeseries_length - 1, ip] - SAS_lookup[timeseries_length - 1, ip]

                printcheck(rdf, rdfp, 'dsTdSj', 'sT', Q_0, ip)
                printcheck(rdf, rdfp, 'dmTdSj', 'mT', Q_0 * C_J, ip)
                for iqq in range(2):
                    for s in range(2):
                        # print(iq, ic, j, ip, iqq, s)
                        printcheckC(rdf, rdfp, 'dCdSj', 'C_Q', ip, iqq, s)

# def test_Jacobian():
#    n = np.arange(timeseries_length)
#    T_0 = S_0 / Q_0
#    Delta = dt / T_0
#    Kappa = np.exp(-Delta)
#    Eta = Kappa ** n
#
#    dsTdSjdisc = -((Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta))
#    dSTdSjdisc = (Q_0 - Q_0 * (1 + Delta + n * Delta) * Eta * Kappa) / (S_0 * Delta)
#    dpQdSjdisc = (Q_0 * (-1 + (1 + Delta + n[0] * Delta) * Eta * Kappa)) / (S_0 ** 2 * Delta ** 2)
#    dpQdSjdisc[0] = (Q_0 * Eta[0] * (-1 + Kappa) * (-1 + Kappa + Delta * (1 + n[0] * (-1 + Kappa) + Kappa))) / (
#            S_0 ** 2 * Delta ** 2 * Kappa)
#    dPQdSjdisc = (Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 ** 2 * Delta ** 2)
#    dPQdSjdisc[0] = (Q_0 * (-1 + Kappa + Delta * Kappa)) / (S_0 ** 2 * Delta ** 2)
#    dmQdSjdisc = dpQdSjdisc * C_J * Q_0
#    dmTdSjdisc = -((C_J * Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta))
#    dCQdSjdisc = ((C_J + C_old) * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta)
#
#    dsTdSjdisc = np.tril(np.tile(dsTdSjdisc, [timeseries_length, 1])).T
#    dsTdSjdisc = np.c_[np.zeros(timeseries_length), dsTdSjdisc]
#    dCQdSjdisc = np.array([[dCQdSjdisc]]).T
#    dpQdSjdisc = np.tril(np.tile(dpQdSjdisc, [timeseries_length, 1]))
#    dpQdSjdisc = np.array([dpQdSjdisc]).T
#    dmQdSjdisc = np.tril(np.tile(dmQdSjdisc, [timeseries_length, 1]))
#    dmQdSjdisc = np.array([[dmQdSjdisc]]).T
#    dmTdSjdisc = np.tril(np.tile(dmTdSjdisc, [timeseries_length, 1])).T
#    dmTdSjdisc = np.c_[np.zeros(timeseries_length), dmTdSjdisc]
#    dmTdSjdisc = np.array([dmTdSjdisc.T]).T
#
#    ST_min = 1.5
#    model = steady_run(timeseries_length, dt, Q_0, S_0, C_J, ST_min=ST_min)
#    rdf = model.result
#    modeps = [[0, 0], [0, 0]]
#    rdfeps = [[0, 0], [0, 0]]
#    modeps[0][0] = steady_run(timeseries_length, dt, Q_0, S_0, C_J, eps1=[eps, 0], ST_min=ST_min)
#    modeps[0][1] = steady_run(timeseries_length, dt, Q_0, S_0, C_J, eps1=[0, eps], ST_min=ST_min)
#    # modeps[1][0] = steady_run(timeseries_length, dt, Q_0, S_0, C_J, eps2=[eps, 0], ST_min=ST_min)
#    # modeps[1][1] = steady_run(timeseries_length, dt, Q_0, S_0, C_J, eps2=[0, eps], ST_min=ST_min)
#    for i in range(1):
#        for j in range(2):
#            rdfeps[i][j] = modeps[i][j].result
#    Jac = model.get_jacobian(mode='endpoint', logtransform=False)
#    print('Jacobian: calculated')
#    print(Jac[:, :2])
#    print('Jacobian: should be')
#    J = None
#    for isol, sol in enumerate(model._solorder):
#        if 'observations' in model.solute_parameters[sol]:
#            for isolflux, solflux in enumerate(model._fluxorder):
#                if solflux in model.solute_parameters[sol]['observations']:
#                    J_seg = None
#                    for iflux, flux in enumerate(model._comp2learn_fluxorder):
#                        for ilabel, label in enumerate(model.sas_blends[flux]._comp2learn_componentorder):
#                            i = iflux
#                            J_seg_this = np.column_stack(
#                                [((rdfeps[i][j]['C_Q'][:, isolflux, isol] - rdf['C_Q'][:, isolflux, isol]) / eps) for j
#                                 in range(2)])
#                            if J_seg is None:
#                                J_seg = J_seg_this
#                            else:
#                                J_seg = np.c_[J_seg, J_seg_this]
#                    if J is None:
#                        J = J_seg
#                    else:
#                        J = np.concatenate((J, J_seg), axis=0)
#    print(J)
