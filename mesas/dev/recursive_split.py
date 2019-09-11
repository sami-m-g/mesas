import inspect

import numpy as np
from numpy.linalg import inv

from mesas.sas.functions import Piecewise


def _verbose(string):
    print(string)


def run(model, components_to_learn=None, ST_min=0., ST_max=1., **kwargs):
    """
    Estimates a SAS function that reproduces observations using a piecewise constant pdf

    :param model: a sas old_model instance
    :param components_to_learn: dict of list of str
    :param ST_min:
    :param ST_max:
    :param kwargs:
    :return:
    """
    _verbose(f'Starting {inspect.stack()[0][3]}')
    initial_model = model.copy_without_results()
    if components_to_learn is None:
        components_to_learn = model.get_component_labels()
    new_components = {}
    for flux in components_to_learn.keys():
        new_components[flux] = dict((label, None) for label in components_to_learn[flux])
    for flux in new_components.keys():
        for label in new_components[flux].keys():
            sas_fun = Piecewise([ST_min, ST_max])
            initial_model.set_sas_fun(flux, label, sas_fun)
    models, mse = fit_model(initial_model, **kwargs)
    initial_model = models[-1]
    initial_mse = mse[-1]
    _verbose('Initial old_model is:')
    _verbose(initial_model.sas_blends)
    return increase_resolution(initial_model, initial_mse, new_components, 0, **kwargs)


def increase_resolution(initial_model, initial_mse, new_components, segment, incres_plot_fun=None, **kwargs):
    _verbose(f'\n\n***********************\nStarting {inspect.stack()[0][3]}')
    _verbose(f'Testing increased SAS resolution in segment {segment}')
    new_components_count = 0
    mse_dict = {}
    Ns = len(initial_model.get_segment_list())
    for flux in new_components.keys():
        for label in new_components[flux].keys():
            model_list, mse_list = fit_model(initial_model.subdivided_copy(flux, label, segment), **kwargs)
            new_mse = mse_list[-1]
            mse_dict[f'{flux}, {label}'] = mse_list
            _verbose(f'Initial mse = {initial_mse}')
            _verbose(f'New mse     = {new_mse}')
            if new_mse / initial_mse < 0.9:
                _verbose(f'Subdivision accepted for {flux}, {label}')
                new_components[flux][label] = model_list[-1].sas_blends[flux].components[label]
                new_components_count += 1
            else:
                _verbose(f'Subdivision rejected for {flux}, {label}')
                new_components[flux][label] = None
    if new_components_count > 0:
        if new_components_count > 1:
            new_model = initial_model.copy_without_results()
            for flux in new_components.keys():
                for label in new_components[flux].keys():
                    if new_components[flux][label] is not None:
                        new_model.set_component(flux, new_components[flux][label])
            model_list, mse_list = fit_model(new_model, **kwargs)
        cumulative_new_model = model_list[-1]
        cumulative_new_mse = mse_list[-1]
        _verbose('New model is:')
        _verbose(cumulative_new_model)
        if incres_plot_fun is not None:
            incres_plot_fun(initial_model, cumulative_new_model, mse_dict, Ns, segment)
        return increase_resolution(cumulative_new_model, cumulative_new_mse, new_components, segment,
                                   incres_plot_fun=incres_plot_fun)
    else:
        more_segments = False
        for flux in list(new_components.keys()).copy():
            for label in list(new_components[flux].keys()).copy():
                if segment == initial_model.sas_blends[flux].components[label].sas_fun.nsegment - 1:
                    _verbose(f'No more refinement for {flux}, {label}')
                    del new_components[flux][label]
                else:
                    more_segments = True
        if more_segments:
            _verbose(f'Moving to segment {segment + 1}')
            return increase_resolution(initial_model, initial_mse, new_components, segment + 1)
        else:
            _verbose('Finished increasing old_model resolution')
            return initial_model


def fit_model(model, learn_plot_fun=None, maxiter=20, **kwargs):
    _verbose(f'Starting {inspect.stack()[0][3]}')
    model.run()
    models = [model]
    mse = [np.mean(model.get_residuals() ** 2)]
    for iter in range(maxiter):
        new_model = step(models[-1], **kwargs)
        if iter > 0:
            models[-1].result = None
        models.append(new_model)
        mse.append(np.mean(new_model.get_residuals() ** 2))
    if learn_plot_fun is not None:
        learn_plot_fun(models, mse)
    return models, mse


def step(model, step_plot_fun=None, **kwargs):
    delta, gn_update_info = calc_gauss_newton_update(model, **kwargs)
    new_model = apply_gauss_newton_update(model, delta, **kwargs)
    new_model.run()
    if step_plot_fun is not None:
        step_plot_fun(model, new_model, gn_update_info)
    return new_model


def calc_gauss_newton_update(model, alpha_step=1, max_delta=None, **kwargs):
    # extract the timeseries of predictions
    ri = model.get_residuals()
    index = np.nonzero(~np.isnan(ri))[0]
    J = model.get_jacobian(index=index)
    Nt, Np = J.shape
    #
    # Gauss-Newton Algorithm
    ri[np.isnan(ri)] = 0
    gradient = 2 * np.dot(J.T, ri)
    gradient[np.isnan(gradient)] = 0
    p_not_flat_gradient = gradient != 0
    H = 2 * np.dot(J[:, p_not_flat_gradient].T, J[:, p_not_flat_gradient])
    delta = np.zeros(Np)
    delta[p_not_flat_gradient] = -np.dot(inv(H), gradient[p_not_flat_gradient])
    if max_delta is not None:
        delta = np.where(np.abs(delta) > max_delta, max_delta * np.sign(delta), delta)
    #
    gradient_ip = np.zeros((Nt, Np))
    gradient_ip[:, p_not_flat_gradient] = 2 * (J[:, p_not_flat_gradient].T * ri).T
    delta_ip = np.zeros((Nt, Np))
    delta_ip[:, p_not_flat_gradient] = -np.dot(inv(H), gradient_ip[:, p_not_flat_gradient].T).T
    #
    params_old = model.get_segment_list()
    params_new = params_old * np.exp(alpha_step * delta[:len(params_old)])
    #
    return delta, {
        'H': H,
        'delta': delta,
        'delta_ip': delta_ip,
        'gradient': gradient,
        'gradient_ip': gradient_ip,
        'params_new': params_new
    }


def apply_gauss_newton_update(model, delta, alpha_step=1, include_C_old=True, **kwargs):
    # Update the parameters
    params_old = model.get_segment_list()
    params_new = params_old * np.exp(alpha_step * delta[:-1])
    new_model = model.copy_without_results()
    new_model.update_from_segment_list(params_new)
    if include_C_old:
        for isol, sol in enumerate(model._solorder):
            old_C_old = model.solute_parameters[sol]['C_old']
            new_C_old = old_C_old + alpha_step * delta[-model._numsol + isol]
            new_model.solute_parameters[sol]['C_old'] = new_C_old
    return new_model
