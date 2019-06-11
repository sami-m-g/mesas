from copy import deepcopy
from sas_model import Model, SASTimeseries
from sas_functions import Piecewise
import pandas as pd
import numpy as np
from scipy.optimize import fmin
import matplotlib.pyplot as plt
plt.ion()
np.random.seed(1)

'''
    Import data
'''
csvfile = '../data/lower_hafren.csv'
data_df = pd.read_csv(csvfile, index_col=0, header=0, parse_dates=False)
N = len(data_df)  # The number of observations

'''
    Initialize the model
'''
# create the SAS functions
# initially these are a single piece -- i.e. a uniform distribution
sas_fun_Q_0 = Piecewise(npiece=1, ST_max=3533.)
sas_fun_E = Piecewise(npiece=1, ST_max=2267.)
# Associate the SAS functions with fluxes, and create timeseries of SAS functions
sas_ts = {
    'Q': SASTimeseries(sas_fun_Q_0, N=N),
    'ET': SASTimeseries(sas_fun_E, N=N)
}
# Specify parameters of the solutes to be transported
solute_parameters = {
    'Cl mg/l': {'C_old': 7.11, 'alpha': {'Q': 1., 'ET': 0.}}
}
# Create the model
mymodel = Model(
    data_df=data_df,
    sas_ts=sas_ts,
    solute_parameters=solute_parameters,
    n_substeps=3,
    verbose=True,
    debug=False,
    full_outputs=False
)


def run(params, model):
    """Helper function to run the model with a new set of parameters

    """
    # update the SAS function timeseries with the new parameters
    model.sas_ts['Q'].update_from_paramlist(params)
    # run the model
    model.run()
    # extract the timeseries of predictions
    pred = model.result['C_Q'][:, 0, 0]
    return pred


def objective_function(params, *args):
    """Evaluates the objective function for a given set of parameters

    """
    # run the helper function
    pred = run(params, *args)
    # define observations
    # find out where the observation has values
    obs_index = np.logical_not(np.isnan(data_df['Q Cl mg/l']))
    # select only those values
    pred = pred[obs_index]
    obs = data_df['Q Cl mg/l'].loc[obs_index]
    # calculate the RMSE
    err = np.sqrt(np.mean((obs - pred)**2))
    print(err, params)
    return err


'''
    Model 0 -- one piece
'''
if False:
    # Retrieve the current parameters as a list
    params_0 = mymodel.sas_ts['Q'].get_paramlist()
    # Use the built in optimizer to minimize the objective function
    #params_0 = fmin(objective_function, params_0, args=(mymodel, ))
    # get the prediction using the optimized parameters
    C_Q_pred = run(params_0, mymodel)

    # plot the result
    plt.figure(0)
    plt.clf()
    plt.plot(data_df['Q Cl mg/l'], 'r.')
    plt.plot(C_Q_pred, 'b')
    plt.show()

'''
    Model 1 -- two pieces
'''
if False:
    mymodel.results = None
    sas_fun_Q_1 = Piecewise(npiece=2, ST_max=3533.)
    mymodel.set_sas_ts('Q', SASTimeseries(sas_fun_Q_1, N=N))

    # Retrieve the current parameters as a list
    params_1 = mymodel.sas_ts['Q'].get_paramlist()
    # Use the built in optimizer to minimize the objective function
    #params_1 = fmin(objective_function, params_1, args=(mymodel, ))
    # get the prediction using the optimized parameters
    C_Q_pred = run(params_1, mymodel)

    # plot the result
    plt.figure(1)
    plt.clf()
    plt.plot(data_df['Q Cl mg/l'], 'ro')
    plt.plot(C_Q_pred, 'b')
    plt.show()

'''
    Model 2 -- two pieces, two states
'''
if True:
    mymodel.results = None
    # Create two SAS functions, one for high flows, one for low
    sas_fun_Q_3_0 = Piecewise(npiece=2, ST_max=3533.)
    sas_fun_Q_3_1 = Piecewise(npiece=2, ST_max=3533.)
    # Define a scalar that corresponds with the 'state' of the system
    state_ts = np.where(data_df['Q'] > data_df['Q'].median(), 1., 0.)
    # set the SAS function to switch between these SAS functions
    mymodel.set_sas_ts(
        'Q', SASTimeseries({0.: sas_fun_Q_3_0,
                            1.: sas_fun_Q_3_1},
                           state_ts=state_ts))

    # Retrieve the current parameters as a list
    params_2 = mymodel.sas_ts['Q'].get_paramlist()
    # Use the built in optimizer to minimize the objective function
    params_2 = fmin(objective_function, params_2, args=(mymodel, ))
    # get the prediction using the optimized parameters
    C_Q_pred = run(params_2, mymodel)

    # plot the result
    plt.figure(2)
    plt.clf()
    plt.plot(data_df['Q Cl mg/l'], 'ro')
    plt.plot(C_Q_pred, 'b')
    plt.show()
