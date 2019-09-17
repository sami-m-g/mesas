#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
   Example 2 -- Lower Hafren

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
from mesas.sas.blender import Fixed, Weighted
# classes we will use to construct the SAS model
# 1. A piecewise constant SAS function
from mesas.sas.functions import Piecewise
# 3. The sas model class
from mesas.sas.model import Model

# this makes the plots we create interactive
plt.ion()
# Set the random seed
np.random.seed(1)

# These are the column headers of important input timeseries
influx_name = 'J'  # volume/time input of water
outflux_name = 'Q'
sol_name = 'Cl mg/l'  # mass/volume of input solute
obs_name = 'Q Cl mg/l'  # mass/volume of solute outputs <-- this is the data we want to reproduce by 'training' a SAS function
reference_model = None

csvfile = '../data/lower_hafren.csv'
data_df = pd.read_csv(csvfile, index_col=0, header=0, parse_dates=False)  # [:1000]
# Increase the input flux to account for occult precipitation
N = len(data_df)  # The number of observations

# Parameter used when plotting SAS functions. Sets the x-axis limits
plot_ST_max = 3000

# Create a new sas model
# Our initial guess of the shape of the sas function
weights_df = pd.DataFrame(index=data_df.index)
weights_df['max'] = (data_df['Q'].rank(pct=True))
weights_df['min'] = 1 - weights_df['max']
# sas_fun_Q_min = Piecewise(ST=[0., 10000., 20000.], P=[0, 0.5, 1.0])
# sas_fun_Q_max = Piecewise(ST=[0, 117.8, 362.0, 1.612E+3], P=[0, 0.25, 0.5, 1.0])
sas_fun_Q_min = Piecewise(segment_list=[0., 8000.])
sas_fun_Q_max = Piecewise(segment_list=[0., 200.])
sas_fun_E = Piecewise(nsegment=1, ST_max=398.)
my_sas_blends = {
    'Q': Weighted({'max': sas_fun_Q_max,
                   'min': sas_fun_Q_min},
                  weights_df=weights_df),
    'ET': Fixed(sas_fun_E, N=N)
}

# Specify parameters of the solutes to be transported
solute_parameters = {
    'Cl mg/l': {
        'C_old': 7.11 / 1.,  # Concentration for water of unknown age
        'alpha': {'Q': 1., 'ET': 0.},  # partitioning coeff. Accounts for lack of Cl in ET
        'observations': {
            'Q': obs_name
        }}}

# Create the model
mymodel = Model(
    data_df=data_df,
    sas_blends=my_sas_blends,
    solute_parameters=solute_parameters,
    n_substeps=4,  # substeps increase numerical accuracy
    verbose=False,  # print information about calculation progress
    debug=False,  # print (lots of) information
    full_outputs=True,  # calculate and return mass balance information (slow)
    influx='J',  # label of the water influx data in data_df
    max_age=365 * 7  # this limits the memory of the past retained in the model, but reduces computational load
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
    # ax3 = plt.subplot2grid((3, 4), (2, 2), colspan=2)
    # plots.plot_MSE_improvement(mse_dict, ax3)
    C_train = old_model.data_df[obs_name]
    vartrain = np.nanvar(C_train)
    # ax3.set_ylim((vartrain / 100000., vartrain * 3))
    #
    ax4 = plt.subplot2grid((3, 4), (1, 2), colspan=2)
    plots.plot_residuals_timeseries(old_model, new_model, outflux_name, sol_name, ax4)
    ax4.set_ylim((vartrain / 100000., vartrain * 3))
    #
    plt.show()
    plt.tight_layout()
    figname = f'../junk/plots/LH_{Ns}_{segment}.png'
    plt.savefig(figname)


# Run the recursive_split algorithm
rs_model = recursive_split(mymodel,
                           incres_plot_fun=incres_plot_fun,
                           components_to_learn={'Q': ['max', 'min']})
