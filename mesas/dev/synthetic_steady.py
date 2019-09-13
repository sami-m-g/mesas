#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
   Example 1 -- Synthetic steady

   In this example, inflow and one outflow from a control volume are held steady, and a timeseries of
   random inputs are generated and run through the SAS model with a fixed SAS function.

   The recursive_split algorithm is then used to reconstruct the assumed SAS function from the
   output concentrations.

"""
# First we import some things

# library imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# a library of plotting functions that can be used to visualize the output and internal steps
import plots
# the mesas algorithm used
from recursive_split import run as recursive_split

# 2. A 'blender' that assumes that SAS function is fixed in time
from mesas.sas.blender import Fixed
# classes we will use to construct the SAS model
# 1. A piecewise constant SAS function
from mesas.sas.functions import Piecewise
# 3. The sas model class
from mesas.sas.model import Model

# this makes the plots we create interactive
plt.ion()
# Set the random seed
np.random.seed(2)

# These are the column headers of important input timeseries
influx_name = 'J'  # volume/time input of water
sol_name = 'C_J'  # mass/volume of input solute
outflux_name = 'Q'  # volume/time output of water
obs_name = 'C_Q obs_name'  # mass/volume of solute outputs <-- this is the data we want to reproduce by 'training' a SAS function

# Parameter used when plotting SAS functions. Sets the x-axis limits
plot_ST_max = 800


def create_synthetic_data(N, influx, gamma_scale, gamma_shape, sampling_interval, C_mean, C_sigma):
    '''
    Creates a synthetic dataset 
    
    Inflows and outflows are assumed to be steady, and a gamma distribution is assumed for the SAS function

    :return: A dataframe containing the synthetic data, and an instance of a sas model used to create it
    '''
    #
    # Create a pandas dataframe to hold the timeseries
    synth_data_df = pd.DataFrame(index=np.arange(N))
    # and populate it with synthetic input data
    synth_data_df[influx_name] = influx
    synth_data_df[outflux_name] = influx
    synth_data_df[sol_name] = C_mean + C_sigma * np.random.randn(N)
    #
    # Create a sas model with a gamma SAS function
    from scipy.special import gammaincinv
    # The segment list is a piecewise-linear approximation of the CDF
    # Each segment is 1% of the probability
    segment_list = np.zeros(100 + 1)
    segment_list[1:] = np.diff(gamma_scale * np.abs(gammaincinv(gamma_shape, np.linspace(0, 1, 100 + 2)[:-1])))
    # the first item in the segment list is the offset from zero (if any)
    segment_list[0] = 0
    # This creates the piecewise constant SAS function from the segment list
    synth_sas_fun_Q = Piecewise(segment_list=segment_list)
    # A sas blend is an object that can be queried to get the sas function at each timestep.
    # Here we assume the sas function is fixed in time.
    # The sas blends are packaged in a dictionary that associates each one with the
    # name of an outflux timeseries in the dataframe
    synth_sas_blends = {outflux_name: Fixed(synth_sas_fun_Q, N=N)}
    # The solute parameters are packaged in a dictionary with keys identical to the
    # name of the input concentration timeseries in the dataframe
    synth_solute_parameters = {
        sol_name: {
            'C_old': C_mean,  # Concentration for water of unknown age
            'alpha': {outflux_name: 1.}}  # partitioning coeff
    }
    # Create the model
    synth_truth_model = Model(
        data_df=synth_data_df,
        sas_blends=synth_sas_blends,
        solute_parameters=synth_solute_parameters,
        n_substeps=4,  # substeps increase numerical accuracy
        verbose=False,  # print information about calculation progress
        debug=False,  # print (lots of) information
        full_outputs=True,  # calculate and return mass balance information (slow)
        influx='J'  # label of the water influx_name data in data_df
    )
    #
    # run the model an generate results
    synth_truth_model.run()
    #
    # Create a timeseries of periodic `samples` spaced by NaNs
    C_train = np.zeros(N) * np.NaN
    C_train[::sampling_interval] = synth_truth_model.result['C_Q'][:, 0, 0][::sampling_interval]
    # add this to the dataframe
    synth_data_df[obs_name] = C_train
    #
    return synth_data_df, synth_truth_model


# Run the function to generate a synthetic dataset
data_df, reference_model = create_synthetic_data(
    N=1 * 365,
    gamma_scale=700,
    gamma_shape=0.5,
    influx=2300 / 365.,
    sampling_interval=1,
    C_mean=10.05,
    C_sigma=1.
)

# Create a new sas model
# This model will be trained to reproduce the observations generated above
# This creates the piecewise constant SAS function with one segment from (0,0) to (1000, 1)
# It is our initial guess of the shape of the sas function
my_sas_fun_Q = Piecewise(segment_list=[0, 1000])

# A sas blend is an object that can be queried to get the sas function at each timestep.
# Here we assume the sas function is fixed in time for as many timesteps as the input timeseries
# The sas blends are packaged in a dictionary that associates each one with the
# name of an outflux timeseries in the dataframe
my_sas_blends = {
    outflux_name: Fixed(my_sas_fun_Q, N=len(data_df))
}

# The solute parameters are packaged in a dictionary with keys identical to the
# name of the input concentration timeseries in the dataframe
my_solute_parameters = {
    sol_name: {
        'C_old': 10.0,  # Concentration for water of unknown age
        'alpha': {outflux_name: 1.},  # partitioning coeff
        'observations': {
            outflux_name: obs_name
        }}}

# Create the model
mymodel = Model(
    data_df=data_df,
    sas_blends=my_sas_blends,
    solute_parameters=my_solute_parameters,
    n_substeps=4,  # substeps increase numerical accuracy
    verbose=False,  # print information about calculation progress
    debug=False,  # print (lots of) information
    full_outputs=True,  # calculate and return calculation arrays, including mass balance (slower)
    influx=influx_name  # column name of the water influx_name data in data_df
)


# Create a function that will be used to make some plots each time the sas function resolution is increased
def incres_plot_fun(old_model, new_model, mse_dict, Ns, segment):
    '''

    :param old_model: a sas model with lower resolution
    :param new_model: a sas model with higher resolution (one segment split in two)
    :param mse_dict: a dictionary showing how the mse changed during convergence to the new model
    :param Ns: the number of segments in the old model
    :param segment: the segment that was split in two
    '''
    fig1 = plt.figure(figsize=[16, 8])
    #
    ax1 = plt.subplot2grid((3, 4), (0, 0), rowspan=3, colspan=2)
    plots.plot_SAS_update(old_model, new_model, reference_model, plot_ST_max, ax1)
    #
    ax2 = plt.subplot2grid((3, 4), (0, 2), colspan=2)
    plots.plot_timeseries_update(old_model, new_model, outflux_name, sol_name, ax2)
    #
    ax3 = plt.subplot2grid((3, 4), (2, 2), colspan=2)
    plots.plot_MSE_improvement(mse_dict, ax3)
    C_train = old_model.data_df[obs_name]
    vartrain = np.nanvar(C_train)
    ax3.set_ylim((vartrain / 100000., vartrain * 3))
    #
    ax4 = plt.subplot2grid((3, 4), (1, 2), colspan=2)
    plots.plot_residuals_timeseries(old_model, new_model, outflux_name, sol_name, ax4)
    ax4.set_ylim((vartrain / 100000., vartrain * 3))
    #
    plt.show()
    plt.tight_layout()
    figname = f'../junk/plots/synth_{Ns}_{segment}.png'
    plt.savefig(figname)


# Run the recursive_split algorithm
rs_model = recursive_split(mymodel,
                           incres_plot_fun=incres_plot_fun,
                           alpha_step=0.5,
                           max_delta=np.log(1.5))
