import pytest
import pandas as pd
import numpy as np
from mesas.sas.model import Model

def test_testtest():
    N = 1000 # number of timesteps to run
    t = np.arange(N)

    J = 2.5 # inflow rate
    k = 1/20 # hydraulic response rate
    Q = J * (1 - np.exp(-k * t))
    S = Q / k
    S[0] = S[1]/1000

    Cin = np.zeros(N)
    Cin[10] = 100./J

    my_dataframe = pd.DataFrame(index = t, data={'J':J, 'Q':Q, 'S':S, 'Cin':Cin})

    my_sas_specs = {'Q':
                    {'a uniform distribution over total storage':
                         {'ST': [0, 'S'] }}}

    my_solute_parameters = {'Cin':{}}

    my_model = Model(data_df=my_dataframe, sas_specs=my_sas_specs, solute_parameters=my_solute_parameters)

    my_model.run()

    print(('Cin --> Q' in my_model.data_df.columns) & (not np.any(np.isnan(my_model.data_df['Cin --> Q']))), my_model.data_df['Cin --> Q'].sum())

