#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
   Example 2 -- Lower Hafren

"""
# First we import some things

import pickle

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
data_df = pd.read_csv(csvfile, index_col=0, header=0, parse_dates=False)[:5721]
index = np.nonzero(np.isnan(data_df[obs_name]))[0]
index = index[index > 2070]
index = index[index < 5721]

# Increase the input flux to account for occult precipitation
N = len(data_df)  # The number of observations

# Parameter used when plotting SAS functions. Sets the x-axis limits
plot_ST_max = 2000

# Create a new sas model
# Our initial guess of the shape of the sas function
weights_df = pd.DataFrame(index=data_df.index)
weights_df['max'] = (data_df['Q'].rank(pct=True))
weights_df['min'] = 1 - weights_df['max']
sas_fun_Q_min = Piecewise(ST=[0., 1.637e+04], P=[0, 1.0])
sas_fun_Q_max = Piecewise(ST=[0, 48.81, 122.5, 228.2, 373.1, 1.729E+3], P=[0, 0.125, 0.25, 0.375, 0.5, 1.0])
# sas_fun_Q_min = Piecewise(segment_list=[0., 8000.])
# sas_fun_Q_max = Piecewise(segment_list=[0., 200.])
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
    components_to_learn={'Q': ['max', 'min']},
    n_substeps=4,  # substeps increase numerical accuracy
    verbose=False,  # print information about calculation progress
    debug=False,  # print (lots of) information
    full_outputs=True,  # calculate and return mass balance information (slow)
    influx='J',  # label of the water influx data in data_df
    max_age=365 * 1  # this limits the memory of the past retained in the model, but reduces computational load
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
    ax1 = plt.subplot2grid((3, 4), (0, 0), rowspan=2, colspan=4)
    plots.plot_SAS_update(old_model, new_model, reference_model, plot_ST_max, ax1)
    #
    ax2 = plt.subplot2grid((3, 4), (2, 0), colspan=4)
    plots.plot_timeseries_update(old_model, new_model, outflux_name, sol_name, ax2)
    #
    # ax3 = plt.subplot2grid((3, 4), (2, 2), colspan=2)
    # plots.plot_MSE_improvement(mse_dict, ax3)
    # C_train = old_model.data_df[obs_name]
    # vartrain = np.nanvar(C_train)
    # ax3.set_ylim((vartrain / 100000., vartrain * 3))
    #
    # ax4 = plt.subplot2grid((3, 4), (1, 2), colspan=2)
    # plots.plot_residuals_timeseries(old_model, new_model, outflux_name, sol_name, ax4)
    #ax4.set_ylim((vartrain / 100000., vartrain * 3))
    #
    plt.show()
    plt.tight_layout()
    figname = f'../junk/plots/LH_{mode}_{Ns}_{segment}.png'
    plt.savefig(figname)
    with open(f'../junk/LH_{mode}_{Ns}_{segment}.pickle', 'wb') as f:
        pickle.dump(new_model.copy_without_results(), f, pickle.HIGHEST_PROTOCOL)


## Run the recursive_split algorithm
# mode = 'analytical'
# rs1_model = recursive_split(mymodel,
# incres_plot_fun=incres_plot_fun,
# jacobian_mode=mode,
# )
#
# rs1_model.results = None
# with open('LH_rs1.pickle', 'wb') as f:
# pickle.dump(rs1_model, f, pickle.HIGHEST_PROTOCOL)
with open('LH_rs1.pickle', 'rb') as f:
    rs1_model = pickle.load(f)

rs1_model.options = {'max_age': 365 * 10}
rs1_model.components_to_learn = {'Q': ['max', 'min'], 'ET': ['Fixed']}
mode = 'numerical'
rs2_model = recursive_split(rs1_model,
                            incres_plot_fun=incres_plot_fun,
                            jacobian_mode=mode,
                            components_to_learn={'Q': ['max', 'min'], 'ET': ['Fixed']})

with open('LH_rs2.pickle', 'wb') as f:
    pickle.dump(rs2_model, f, pickle.HIGHEST_PROTOCOL)
