#Todo: incorporate pytest-benchmark
import numpy as np
import pandas as pd

# 2. A 'blender' that assumes that SAS function is fixed in time
from mesas.sas.blender import Fixed
# classes we will use to construct the SAS model
# 1. A piecewise constant SAS function
from mesas.sas.functions import Piecewise
# 3. The sas model class
from mesas.sas.model import Model


#@pytest.fixture
def steady_df():
    # Steady-state flow in and out for N timesteps
    Q_0 = 1.0  # <-- steady-state flow rate
    C_J = 100.0
    N = 10
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(N) * Q_0
    data_df['Q2'] = np.ones(N) * Q_0
    data_df['J'] = np.ones(N) * Q_0 * 2
    data_df['C_J'] = np.ones(N) * C_J
    return data_df


def test_steady(steady_df):
    data_df = steady_df
    S_0 = 4.
    sas_fun1 = Piecewise(nsegment=1, ST_max=S_0)
    sas_fun2 = Piecewise(nsegment=1, ST_max=S_0)
    sas_blends = {'Q1': Fixed(sas_fun1, N=len(data_df)), 'Q2': Fixed(sas_fun2, N=len(data_df))}
    solute_parameters = { 'C_J': { }}
    model = Model(data_df, sas_blends, solute_parameters, debug=True, verbose=True)
    print('running test_steady')
    model.run()


test_steady(steady_df())


