import numpy as np


def plot_SAS_update(old_model, new_model, reference_model, plot_ST_max, ax1):
    for flux in old_model._fluxorder:
        if reference_model:
            reference_model.sas_blends[flux].plot(ax=ax1, zorder=-10, color='r')
        for label in old_model.sas_blends[flux].components.keys():
            old_component, new_component = [m.sas_blends[flux].components[label] for m in [old_model, new_model]]
            old_component.sas_fun.plot(ax=ax1, color='b', label='previous ' + label, alpha=0.5, marker='o')
            new_component.sas_fun.plot(ax=ax1, color='b', ls='-', label='new ' + label, lw=2, marker='o')
            old_ST, new_ST = [c.sas_fun.ST[1:] for c in [old_component, new_component]]
            old_Omega, new_Omega = [c.sas_fun.P[1:] for c in [old_component, new_component]]
            old_ST_interp = np.interp(new_Omega, old_Omega, old_ST)
            old_Omega_interp = new_Omega
            ax1.quiver(old_ST_interp, old_Omega_interp, new_ST - old_ST_interp, new_Omega - old_Omega_interp,
                       scale_units='xy', angles='xy', scale=1, lw=1.5)
    ST = old_model.result['ST']
    ax1.plot([ST[-1, -1], ST[-1, -1]], [0, 1], 'k:')
    ax1.set_ylabel('$\Omega(S_T,t)$')
    ax1.set_xlabel('$S_T$')
    if plot_ST_max in globals():
        ax1.set_xlim((0, plot_ST_max))


def plot_timeseries_update(old_model, new_model, flux, sol, ax2):
    iflux = list(old_model._fluxorder).index(flux)
    isol = list(old_model._solorder).index(sol)
    old_C_pred = old_model.result['C_Q'][:, iflux, isol]
    new_C_pred = new_model.result['C_Q'][:, iflux, isol]
    old_C_old, new_C_old = [m.solute_parameters[sol]['C_old'] for m in [old_model, new_model]]
    ax2.plot(old_C_pred, 'b-', alpha=0.5)
    ax2.plot(new_C_pred, 'b-', lw=2)
    ax2.plot(np.ones_like(old_C_pred) * new_C_old, 'b-', lw=1)
    ax2.plot(np.ones_like(old_C_pred) * old_C_old, 'b-', lw=0.5, alpha=0.3)
    if flux in old_model.solute_parameters[sol]['observations']:
        obs = old_model.solute_parameters[sol]['observations'][flux]
        C_train = old_model.data_df[obs]
        ax2.plot(C_train, 'r.', markersize=4, alpha=0.8)
        if not old_C_old == new_C_old:
            t = np.linspace(0, len(old_C_pred), 50)
            ax2.quiver(t, np.ones_like(t) * old_C_old, np.zeros_like(t), np.ones_like(t) * (new_C_old - old_C_old),
                       scale_units='y', scale=1, lw=1, color='k', alpha=0.3)
    ax2.set_ylabel('Observations')
    ax2.set_xlabel('Time')


def plot_residuals_timeseries(old_model, new_model, flux, sol, ax4):
    iflux = list(old_model._fluxorder).index(flux)
    isol = list(old_model._solorder).index(sol)
    obs = old_model.solute_parameters[sol]['observations'][flux]
    C_train = old_model.data_df[obs]
    old_C_pred = old_model.result['C_Q'][:, iflux, isol]
    new_C_pred = new_model.result['C_Q'][:, iflux, isol]
    old_ri = old_C_pred - C_train
    new_ri = new_C_pred - C_train
    ax4.plot(old_ri ** 2, 'b.-', lw=1, alpha=0.3)
    ax4.plot(new_ri ** 2, 'b.-', lw=1)
    ax4.set_yscale('log')
    vartrain = np.nanvar(C_train)
    ax4.set_ylim((vartrain / 100000., vartrain))
    ax4.set_ylabel('Squared residuals $r_i^2$')
    ax4.set_xlabel('Time')


def plot_MSE_improvement(mse_dict, ax3):
    for label, mse_list in mse_dict.items():
        ax3.plot(np.arange(len(mse_list)) + 1, mse_list, '.-', label=label)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel('$\sum r_i^2/n$')


def plot_expected_vs_actual_change_in_residual(old_model, new_model, flux, sol, alpha_step, ax3):
    old_ri = old_model.get_residuals()
    new_ri = new_model.get_residuals()
    index = np.nonzero(~np.isnan(old_ri))[0]
    J = old_model.jacobian[flux][sol]
    delta = old_model.gn_update_info['delta']
    dri_actual = new_ri - old_ri
    dri_expected = np.zeros_like(new_ri)
    dri_expected[index] = alpha_step * np.dot(J[index, :], delta)
    ax3.plot(dri_expected, dri_actual, 'k.', alpha=0.2)
    ax3.set_xlabel('expected change in residual')
    ax3.set_ylabel('actual change in residual')
    ax3.spines['left'].set_position(('data', 0.0))
    ax3.spines['bottom'].set_position(('data', 0.0))
    ax3.spines['right'].set_color('none')
    ax3.spines['top'].set_color('none')
    ax3.axis('equal')

    # segment_list_old = old_model.sas_blends[flux].get_segment_list()
    # Ns = len(segment_list_old) - 1
    # delta_ip = old_model.gn_update_info['delta_ip']
    # gradient = old_model.gn_update_info['gradient']
    # gradient_ip = old_model.gn_update_info['gradient_ip']
    # J = old_model.jacobian[flux][sol]
    # r_seg, r_old = old_model.get_residual_parts(flux, sol, obs)
    # ri2 = np.sum(r_seg, axis=1) + r_old
    # ri, training_data_index = mymodel.get_residual(flux, sol, obs)
    # S_maxcalc = np.diag(ST)
    # PQ = old_model.result['PQ'][:, :, 0]
    # P_maxcalc = np.diag(PQ)
    # Omegaj = old_model.sas_blends[flux].get_P_list()[1:]
    # J_seg = J[:, :-1]
    ##

    # axJ = plt.subplot2grid((3, 3), (2, 0))
    # plt.pcolormesh(np.arange(len(training_data_index)), np.r_[0, P_list[1:]], J_seg[training_data_index, :Ns + 1].T,
    #               cmap=plt.get_cmap('PiYG'))
    # plt.plot(P_maxcalc[training_data_index], 'r', lw=3)
    # stddelta = iqr(J_seg[training_data_index, :Ns + 1].ravel()) / 2
    # plt.clim((-2 * stddelta, 2 * stddelta))
    # plt.colorbar()
    # plt.title('Jacobian $\partial r_i/\partial \sigma_j$')
    # axr = plt.subplot2grid((3, 3), (2, 1))
    # plt.pcolormesh(np.arange(len(training_data_index)), np.r_[0, P_list[1:]],
    #               gradient_ip[training_data_index, :Ns + 1].T,
    #               cmap=plt.get_cmap('PiYG'))
    # plt.plot(P_maxcalc[training_data_index], 'r', lw=3)
    # stddelta = np.abs(r_seg[training_data_index]).std()
    # plt.clim((-2 * stddelta, 2 * stddelta))
    # plt.colorbar()
    # plt.title(r'Gradient $2r_i\partial r_i/\partial \sigma_j$')
    # axd = plt.subplot2grid((3, 3), (2, 2))
    # plt.pcolormesh(np.arange(len(training_data_index)), np.r_[0, P_list[1:]], delta_ip[training_data_index, :Ns + 1].T,
    #               cmap=plt.get_cmap('PiYG'))
    # plt.plot(P_maxcalc[training_data_index], 'r', lw=3)
    # stddelta = np.abs(delta_ip[training_data_index, :Ns + 1]).std()
    # plt.clim((-2 * stddelta, 2 * stddelta))
    # plt.colorbar()
    # plt.title(r'Change in $\sigma_j$ from $r_i$ (influence $\Delta_ij$)')
    # ax1 = plt.subplot2grid((3, 3), (1, 2))
    # plt.plot(np.arange(len(training_data_index)), delta_ip[training_data_index, -1])
    # plt.title(r'Change in $C_{old}$ from $r_i$')
