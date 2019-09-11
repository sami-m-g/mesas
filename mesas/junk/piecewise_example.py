rom
sas_model
#from scipy.optimize import fmin
# from scipy.optimize import fmin
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sas_blender import Weighted, Fixed
from sas_functions import Piecewise
from sas_model import Model

plt.ion()
np.random.seed(1)

'''
    Define some functions we will use to optimize SAS function parameters
'''


def run(params, model):
    """Helper function to run the old_model with a new set of parameters

    """
    # update the SAS function timeseries with the new parameters
    model.sas_blends['Q'].update_from_segment_list(params)
    # run the old_model
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
    Import some data

    The column labels in the dataframe are important, and must correspond with
    relevant labels provided when initializing the old_model
'''
csvfile = '../data/lower_hafren.csv'
data_df = pd.read_csv(csvfile, index_col=0, header=0, parse_dates=False)
# Increase the input flux to account for occult precipitation
data_df['Cl mg/l'] = data_df['Cl mg/l'] * 1.1
N = len(data_df)  # The number of observations

'''
    Initialize a simple version of the old_model
'''
# create two SAS functions, one for each flux
# initially these are a single piece -- i.e. a uniform distribution
sas_fun_Q_1 = Piecewise(npiece=1, ST_max=3533.)
sas_fun_E = Piecewise(npiece=1, ST_max=2267.)
# Before being used in the old_model, we need to say how the sas function changes in time
# This is handled by the blender classes defined in blender.py
# for now, we will use the 'Fixed' blender
sas_blends = {
    'Q': Fixed(sas_fun_Q_1, N=N),
    'ET': Fixed(sas_fun_E, N=N)
}
# Specify parameters of the solutes to be transported
solute_parameters = {
    'Cl mg/l': {
        'C_old': 7.11/1.1,                  # Concentration for water of unknown age
        'alpha': {'Q': 1., 'ET': 0.}}       # partitioning coeff. Accounts for lack of Cl in ET
}
# Create the old_model
mymodel = Model(
    data_df=data_df,
    sas_blends=sas_blends,
    solute_parameters=solute_parameters,
    n_substeps=1,                           # substeps increase numerical accuracy
    verbose=True,                           # print information about calculation progress
    debug=False,                            # print (lots of) information
    full_outputs=False,                      # calculate and return mass balance information (slow)
    influx='J'                              # label of the water influx data in data_df
)

'''
    Below we will modify and run this base old_model for a few different cases.
    We will only modify and optimize the SAS function for the discharge flux 'Q'
'''

'''
    Case 1 -- Fixed, single piece for Q
'''
if False:
    # Retrieve the current parameters as a list.
    # The method .get_paramlist() collates the parameters of all the component
    # sas functions in the blend into a single list. A sister method
    # .update_from_paramlist(segment_list) can update all the component functions
    # from a similar list
    params_1 = mymodel.sas_blends['Q'].get_segment_list()
    #
    # Use the built in optimizer to minimize the objective function
    # (uncomment this line to run the optimization)
    #params_1 = fmin(objective_function, params_1, args=(mymodel, ))
    #
    # get the prediction using the optimized parameters
    C_Q_pred = run(params_1, mymodel)

    # plot the result
    plt.figure(0, figsize=[10, 6])
    plt.clf()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    mymodel.sas_blends['Q'].plot(ax=ax1)
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    plt.plot(data_df['Q Cl mg/l'], 'r.')
    plt.plot(C_Q_pred, 'b')
    plt.show()

'''
    Case 2 -- Fixed, two pieces for Q

    Here we assume the SAS function for Q (discharge) has two piecewise linear segments
'''
if False:
    # Clear the previous results
    mymodel.results = None
    #
    # Create a new Piecewise function with two parts.
    # The location of the divisions between the parts (i.e. the initial parameters)
    # are selected randomly, starting from the largest before ST_max
    sas_fun_Q_2 = Piecewise(npiece=2, ST_max=3533.)
    #
    # Add this new SAS function to the old_model using the .set_sas_blend() method.
    # This leaves the component function and blend for ET unmodified.
    mymodel.set_sas_blend('Q', Fixed(sas_fun_Q_2, N=N))
    #
    # Retrieve the current parameters as a list
    params_2 = mymodel.sas_blends['Q'].get_segment_list()
    #
    # Use the built in optimizer to minimize the objective function
    #params_2 = fmin(objective_function, params_2, args=(mymodel, ))
    # Or use these parameters I prepared earlier:
    params_2 = [819., 9902.]
    #
    # get the prediction using the optimized parameters
    C_Q_pred = run(params_2, mymodel)

    # plot the result
    plt.figure(1, figsize=[10, 6])
    plt.clf()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    mymodel.sas_blends['Q'].plot(ax=ax1)
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    plt.plot(data_df['Q Cl mg/l'], 'r.')
    plt.plot(C_Q_pred, 'b')
    plt.show()

'''
    Case 3 -- StateSwitch, two pieces for Q

    Here we assume the SAS function is different for low flows versus high flows.
    We can use the StateSwitch blender to switch between them, based on the current
    value of Q. Each component function for Q (discharge) has two piecewise linear
    segments, as before.
'''
if False:
    mymodel.results = None
    # Create two SAS functions, one for high flows, one for low
    sas_fun_Q_3_low = Piecewise(npiece=2, ST_max=3533.)
    sas_fun_Q_3_high = Piecewise(npiece=2, ST_max=3533.)
    # Define a category timeseries that gives the 'state' of the system
    state_ts = np.where(data_df['Q'] > data_df['Q'].median(), 'high', 'low')
    # set the SAS function to switch between these SAS functions
    mymodel.set_sas_blend(
        'Q', StateSwitch({'low': sas_fun_Q_3_low,
                          'high': sas_fun_Q_3_high},
                         state_ts=state_ts))
    #
    # Retrieve the current parameters as a list
    params_3 = mymodel.sas_blends['Q'].get_segment_list()
    #
    # Use the built in optimizer to minimize the objective function
    #params_3 = fmin(objective_function, params_3, args=(mymodel, ))
    # Or use these parameters I prepared earlier:
    params_3 = [1579., 22684.,   594.,  2506., ]
    #
    # get the prediction using the optimized parameters
    C_Q_pred = run(params_3, mymodel)

    # plot the result
    plt.figure(2, figsize=[10, 6])
    plt.clf()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    mymodel.sas_blends['Q'].plot(ax=ax1)
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    plt.plot(data_df['Q Cl mg/l'], 'r.')
    plt.plot(C_Q_pred, 'b')
    plt.show()

'''
    Case 4 -- Weighted, two pieces for Q
'''
if True:
    mymodel.results = None
    # Create two SAS functions, one for high flows, one for low
    sas_fun_Q_4_min = Piecewise(npiece=2, ST_max=3533.)
    sas_fun_Q_4_max = Piecewise(npiece=2, ST_max=3533.)
    # Define a scalar that corresponds with the 'state' of the system
    weights_df = pd.DataFrame(index=data_df.index)
    weights_df['max'] = (data_df['Q'].rank(pct=True))
    weights_df['min'] = 1 - weights_df['max']
    # set the SAS function to switch between these SAS functions
    mymodel.set_sas_blend(
        'Q', Weighted({'min': sas_fun_Q_4_min,
                       'max': sas_fun_Q_4_max},
                      weights_df=weights_df))
    #
    # Retrieve the current parameters as a list
    params_4 = mymodel.sas_blends['Q'].get_segment_list()
    #
    # Use the built in optimizer to minimize the objective function
    # params_4 = fmin(objective_function, params_4, args=(mymodel, ))
    # Or use these parameters I prepared earlier
    params_4 = [1.52787214e+03, 6.11900453e+06, 4.39158178e+02, 1.93276464e+03]
    #
    # get the prediction using the optimized parameters
    C_Q_pred = run(params_4, mymodel)

    # plot the result

    plt.figure(3, figsize=[10, 6])
    plt.clf()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    plt.xscale('log')
    mymodel.sas_blends['Q'].plot(ax=ax1)
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    plt.plot(data_df['Q Cl mg/l'], 'r.')
    plt.plot(C_Q_pred, 'b')
    plt.show()
