import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plots
from recursive_split import run as recursive_split

from mesas.sas.blender import Fixed
from mesas.sas.functions import Piecewise
from mesas.sas.model import Model

plt.ion()
np.random.seed(2)

sol = 'C_J'
obs = 'C_Q obs'
flux = 'Q'
plot_ST_max = 800


def create_synthetic_data():
    N = 1 * 365
    S_scale = 700
    meanflux = 2300 / 365.
    sampling_interval = 1
    synth_data_df = pd.DataFrame(index=np.arange(N))
    synth_data_df['J'] = meanflux
    synth_data_df[flux] = meanflux
    synth_data_df[sol] = np.random.randn(N) + 10.05
    segment_list = np.zeros(100 + 1)
    from scipy.special import gammaincinv
    segment_list[1:] = np.diff(S_scale * np.abs(gammaincinv(0.5, np.linspace(0, 1, 100 + 2)[:-1])))
    segment_list[0] = 0  # 25
    synth_sas_fun_Q = Piecewise(segment_list=segment_list)
    synth_sas_blends = {flux: Fixed(synth_sas_fun_Q, N=N)}
    synth_solute_parameters = {
        sol: {
            'C_old': 10.05,  # Concentration for water of unknown age
            'alpha': {flux: 1.}}  # partitioning coeff
    }
    # Create the old_model
    synth_truth_model = Model(
        data_df=synth_data_df,
        sas_blends=synth_sas_blends,
        solute_parameters=synth_solute_parameters,
        n_substeps=4,  # substeps increase numerical accuracy
        verbose=False,  # print information about calculation progress
        debug=False,  # print (lots of) information
        full_outputs=True,  # calculate and return mass balance information (slow)
        influx='J'  # label of the water influx data in data_df
    )
    synth_truth_model.run()
    C_train = np.zeros(N) * np.NaN
    C_train[::sampling_interval] = synth_truth_model.result['C_Q'][:, 0, 0][::sampling_interval]
    synth_data_df[obs] = C_train
    return synth_data_df, synth_truth_model


data_df, reference_model = create_synthetic_data()

my_sas_fun_Q = Piecewise()
my_sas_blends = {
    flux: Fixed(my_sas_fun_Q, N=len(data_df))
}
my_solute_parameters = {
    sol: {
        'C_old': 10.0,  # Concentration for water of unknown age
        'alpha': {flux: 1.},  # partitioning coeff
        'observations': {
            flux: obs
        }}}
mymodel = Model(
    data_df=data_df,
    sas_blends=my_sas_blends,
    solute_parameters=my_solute_parameters,
    n_substeps=4,  # substeps increase numerical accuracy
    verbose=False,  # print information about calculation progress
    debug=False,  # print (lots of) information
    full_outputs=True,  # calculate and return calculation arrays, including mass balance (slower)
    influx='J'  # column name of the water influx data in data_df
)


def incres_plot_fun(old_model, new_model, mse_dict, Ns, segment):
    fig1 = plt.figure(figsize=[16, 8])
    #
    ax1 = plt.subplot2grid((3, 4), (0, 0), rowspan=3, colspan=2)
    plots.plot_SAS_update(old_model, new_model, reference_model, plot_ST_max, ax1)
    #
    ax2 = plt.subplot2grid((3, 4), (0, 2), colspan=2)
    plots.plot_timeseries_update(old_model, new_model, flux, sol, ax2)
    #
    ax3 = plt.subplot2grid((3, 4), (2, 2), colspan=2)
    plots.plot_MSE_improvement(mse_dict, ax3)
    C_train = old_model.data_df[obs]
    vartrain = np.nanvar(C_train)
    ax3.set_ylim((vartrain / 100000., vartrain * 3))
    #
    ax4 = plt.subplot2grid((3, 4), (1, 2), colspan=2)
    plots.plot_residuals_timeseries(old_model, new_model, flux, sol, ax4)
    ax4.set_ylim((vartrain / 100000., vartrain * 3))
    #
    plt.show()
    plt.tight_layout()
    figname = f'../junk/plots/synth_{Ns}_{segment}.png'
    plt.savefig(figname)


rs_model = recursive_split(mymodel,
                           ST_min=0., ST_max=100.,
                           incres_plot_fun=incres_plot_fun,
                           alpha_step=0.5,
                           maxdelta=np.log(1.5))
