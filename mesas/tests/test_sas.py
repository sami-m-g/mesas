#Todo: incorporate pytest-benchmark
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


#@pytest.fixture
#def steady_df():
    #return data_df


def test_steady():
    #data_df = steady_df
    # Steady-state flow in and out for N timesteps
    dt = 24.
    Q_0 = 1.0/dt  # <-- steady-state flow rate
    C_J = 100.0
    N = 3
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(N) * Q_0/2
    data_df['Q2'] = np.ones(N) * Q_0/2
    data_df['J'] = np.ones(N) * Q_0
    data_df['Ca'] = np.ones(N) * C_J
    data_df['Cb'] = np.zeros(N)
    data_df['Cb'].loc[0] = C_J
    N = len(data_df)
    S_0 = 3.
    sas_fun1 = Piecewise(nsegment=1, ST_max=S_0)
    sas_fun2 = Piecewise(nsegment=1, ST_max=S_0)
    sas_blends = {'Q1': Fixed(sas_fun1, N=len(data_df)), 'Q2': Fixed(sas_fun1, N=len(data_df))}
    solute_parameters = { 'Ca': { }, 'Cb': { }}
    model = Model(data_df, sas_blends, solute_parameters, debug=False, verbose=True, dt=dt, n_substeps=10)
    print('running test_steady')
    model.run()
    rdf = model.result
    T1 = S_0 / Q_0
    CDF = 1 - np.exp(-dt * np.arange(N + 1) / T1)
    F = 1 - CDF * T1 / dt
    sT_analytical = S_0 * np.diff(CDF)/dt
    pQ_analytical = np.zeros_like(sT_analytical)
    pQ_analytical[0] = F[1]/dt
    pQ_analytical[1:] = np.diff(F, n=2)/dt
    sT_analytical = np.tril(np.tile(sT_analytical, [N, 1])).T
    pQ_analytical = np.tril(np.tile(pQ_analytical, [N, 1])).T
    print('Expected:')
    print(sT_analytical)
    print('Got:')
    print(rdf['sT'][:,1:])
    print('Difference:')
    print(sT_analytical-rdf['sT'][:,1:])
    print('Expected:')
    print(pQ_analytical)
    print('Got:')
    for iq in range(1):
        print(f'iq = {iq}')
        print(rdf['pQ'][:, :, iq])
    print('Difference:')
    for iq in range(1):
        print(f'iq = {iq}')
        print(pQ_analytical - rdf['pQ'][:, :, iq])
    print('Water Balance:')
    print(rdf['WaterBalance'])
    print('Solute Balance:')
    for s in range(1):
        print(rdf['SoluteBalance'][:, :, s])
    assert np.abs(rdf['WaterBalance']).max()<1.0E-6
    assert np.abs(rdf['SoluteBalance']).max()<1.0E-6
    assert np.abs(sT_analytical-rdf['sT'][:,1:]).max()<1.0E-6
    assert np.abs(rdf['SoluteBalance']).max()<1.0E-6
    for iq in range(1):
        assert np.abs(pQ_analytical - rdf['pQ'][:, :, iq]).max()<1.0E-6




