import inspect

import numpy as np
from scipy.optimize import least_squares

VERBOSE = True


def _verbose(string):
    """prints progress updates"""
    if VERBOSE:
        print(string)


def run(model, components_to_learn=None, verbose=True, **kwargs):
    """
    Estimates a SAS function that reproduces observations using a piecewise constant pdf

    :param model: a sas model instance
    :param components_to_learn: dict of list of str
    :param verbose: If True, prints progress updates
    :param kwargs:
    :return:
    """

    # handle verbosity
    global VERBOSE
    VERBOSE = verbose
    _verbose(f'Starting {inspect.stack()[0][3]}')

    # make a copy of the input model
    initial_model = model.copy_without_results()

    _verbose('Initial model is:')
    _verbose(initial_model.sas_blends)

    # Create a dictionary that will hold new components (sas functions) as we find them.
    # It will also serve to keep track of which components to keep subdividing (if there is more than 1
    # NOTE: a component is an object that packages together a sas function and a timeseries of its weight
    # along with a label identifying it. It allows us to specify a time-varying SAS function as a
    # weighted sum of fixed SAS functions
    if components_to_learn is None:
        components_to_learn = model.get_component_labels()
    new_components = {}
    for flux in components_to_learn.keys():
        new_components[flux] = dict((label, None) for label in components_to_learn[flux])

    index = initial_model.get_obs_index()

    # Fit the initial model to the observed data
    better_initial_model, better_initial_mse = fit_model(initial_model, index=index, **kwargs)

    # run the algorithm
    return increase_resolution(better_initial_model, better_initial_mse, new_components, 0, index=index, **kwargs)


def increase_resolution(initial_model, initial_mse, new_components, segment, incres_plot_fun=None, index=None,
                        **kwargs):
    """Increases the resolution of a sas function by splitting a piecewise linear segment

    :param initial_model: The model we want to try increasing the resolution of
    :param initial_mse: How well it currently fits the data
    :param new_components: A dict that keeps track of which components we are trying to refine
    :param segment: The index of the segment we are currently refining
    :param incres_plot_fun: An optional function that can produce some plots along the way
    :param kwargs: Additional arguments to be passed to the fit_model function
    :return:
    """

    _verbose(f'\n***********************\nStarting {inspect.stack()[0][3]}')
    _verbose(f'Testing increased SAS resolution in segment {segment}')

    # This keeps track of how many new components we found
    new_components_count = 0
    # This will keep track of convergence information provided by fit_model in case we want to use it for plotting
    mse_dict = {}
    # total number of segments. Used for plotting (actually for naming the saved figure files)
    Ns = len(initial_model.get_segment_list())

    # loop over the components we want to improve, testing whether subdividing each individually leads to
    # substantial change in the fit
    for flux in new_components.keys():
        for label in new_components[flux].keys():

            new_model, new_mse = fit_model(initial_model.subdivided_copy(flux, label, segment), index=index, **kwargs)

            # The current criteria used for deciding whether to accept the subdivision is simple:
            # Just check whether the model improvement is greater than some threshold
            # This can definitely be improved
            _verbose(f'Initial mse = {initial_mse}')
            _verbose(f'New mse     = {new_mse}')
            if new_mse / initial_mse < 0.99:
                _verbose(f'Subdivision accepted for {flux}, {label}')

                # Add to the record of new components, and increment the counter
                new_components[flux][label] = new_model.sas_blends[flux].components[label]
                new_components_count += 1

            else:
                _verbose(f'Subdivision rejected for {flux}, {label}')

                # Record the rejection by setting the dict entry to None
                new_components[flux][label] = None

            # Keep a record of the mse for plotting
            mse_dict[f'{flux}, {label}'] = new_mse

    # if we accepted refined components we want to add them into the model and then recursively
    # try subdividing them too
    if new_components_count > 0:

        # If we found more than one new component, we need to include all of them, then re-run
        # fit_model so that their interactions can be accounted for
        if new_components_count > 1:
            new_model = initial_model.copy_without_results()
            for flux in new_components.keys():
                for label in new_components[flux].keys():
                    if new_components[flux][label] is not None:
                        new_model.set_component(flux, new_components[flux][label])
            new_model, new_mse = fit_model(new_model, index=index, **kwargs)

        _verbose('New model is:')
        _verbose(new_model)

        # Optionally, make some plots
        if incres_plot_fun is not None:
            incres_plot_fun(initial_model, new_model, mse_dict, Ns, segment)

        # Recursively try subdividing the first of the two new segments
        return increase_resolution(new_model, new_mse, new_components, segment,
                                   incres_plot_fun=incres_plot_fun, index=index, **kwargs)

    # If we didn't accept any new subdivisions, we need to move on
    else:

        # First, we check whether each component even has more segments to subdivide
        more_segments = False
        for flux in list(new_components.keys()).copy():
            for label in list(new_components[flux].keys()).copy():

                if segment == initial_model.sas_blends[flux].components[label].sas_fun.nsegment - 1:

                    # no more segments in this component, so delete it from the dict
                    _verbose(f'No more refinement for {flux}, {label}')
                    del new_components[flux][label]

                else:
                    more_segments = True

        # did any?
        if more_segments:

            # If so, recursively call this function on the next segment
            _verbose(f'Moving to segment {segment + 1}')
            return increase_resolution(initial_model, initial_mse, new_components, segment + 1,
                                       incres_plot_fun=incres_plot_fun, index=index, **kwargs)
        else:

            # If not, we are done. Return the model we've found,
            _verbose('Finished increasing old_model resolution')
            return initial_model


def fit_model(model, include_C_old=True, learn_plot_fun=None, index=None, **kwargs):
    if index is None:
        index = model.get_obs_index()

    new_model = model.copy_without_results()

    segment_list = model.get_segment_list()
    # Assume that any segments of zero length are the first segment
    # TODO: check that zero_segments are indeed ST_min
    non_zero_segments = segment_list > 0
    x0 = np.log(segment_list[non_zero_segments])
    nseg_x0 = len(x0)
    if include_C_old:
        x0 = np.r_[x0, [model.solute_parameters[sol]['C_old'] for sol in model._solorder]]
    keep_jac_columns = np.ones(len(segment_list) + model._numsol) == 1
    keep_jac_columns[:len(segment_list)] = non_zero_segments

    def update_parameters(x):

        new_segment_list = np.zeros_like(segment_list)
        new_segment_list[non_zero_segments] = np.exp(x[:nseg_x0])
        new_model.update_from_segment_list(new_segment_list)

        # optionally update the C_old parameter
        if include_C_old:
            for isol, sol in enumerate(model._solorder):
                new_model.solute_parameters[sol]['C_old'] = x[-model._numsol + isol]

    def f(x):
        update_parameters(x)
        new_model.run()
        # optionally do some plotting
        if learn_plot_fun is not None:
            learn_plot_fun(new_model)
        return new_model.get_residuals()[index]

    def jac(x):
        update_parameters(x)
        new_model.run()
        try:
            J = new_model.get_jacobian(index=index)
        except:
            J = new_model.get_jacobian(index=index)
        return J[:, keep_jac_columns]

    OptimizeResult = least_squares(fun=f,
                                   x0=x0,
                                   jac=jac,
                                   verbose=2)

    xopt = OptimizeResult.x
    update_parameters(xopt)
    new_mse = np.mean(OptimizeResult.fun ** 2)

    return new_model, new_mse

# def fit_model(model, learn_plot_fun=None, maxiter=20, use_scipy=False, **kwargs):
#    """
#    Fits the sas model to the data
#
#    At the moment we are using a Newton-Gauss algorithm that takes advantage of the fact that we can
#    calculate an approximate 'time-localized' Jacobian of the model residuals
#
#    :param model: The model to fit to observations
#    :param learn_plot_fun: An optional function to do some plotting at each iteration of the fit
#    :param maxiter: optional maximum number of iterations
#    :param kwargs: additional keyword arguments to be passed to `step`
#    :return: a list of models and a list of their mean square error (mse)
#    """
#
#    if use_scipy:
#        return fit_model_scipy(model, learn_plot_fun=None, **kwargs):
#
#    _verbose(f'Starting {inspect.stack()[0][3]}')
#
#    # Run the model and initialize the output lists
#    model.run()
#    models = [model]
#    mse = [np.mean(model.get_residuals() ** 2)]
#
#    # loop over the iteration steps
#    for iter in range(maxiter):
#
#        # Call the function that actually does the optimization
#        new_model, new_mse = step(models[-1], mse[-1], **kwargs)
#
#        # Clear out the results of the last model to save memory
#        if iter > 0:
#            models[-1].result = None
#
#        # Add the results to the lists
#        models.append(new_model)
#        mse.append(np.mean(new_model.get_residuals() ** 2))
#
#    # optionally do some plotting
#    if learn_plot_fun is not None:
#        learn_plot_fun(models, mse)
#
#    return models, mse
#
#
# def step(model, mse, alpha_step=1, max_delta=None, LM_lambda=0, step_plot_fun=None, **kwargs):
#    """
#    Takes a step in parameter space
#
#    :param model: Model to start with
#    :param mse: starting mean square error
#    :param step_plot_fun: optional plotting function
#    :param kwargs: additional keyword arguments to be passed to the algorithm
#    :return: a new model with modified parameters and results
#    """
#
#    # Take the step
#    new_model, new_mse, gn_info = apply_step(model,
#                                             alpha_step=alpha_step,
#                                             max_delta=max_delta,
#                                             LM_lambda=LM_lambda,
#                                             **kwargs)
#
#    # Optionally, make some plots
#    if step_plot_fun is not None:
#        step_plot_fun(model, new_model, gn_info)
#
#    return new_model, new_mse
#
#
# def calc_gauss_newton_step(model, LM_lambda=0):
#    """
#    Uses a gauss-newton method to choose a step in parameter space
#
#    :param model: sas model
#    :param LM_lambda: Levenberg–Marquardt algorithm lambda
#    :return: an array of delta values
#    """
#
#    # extract the timeseries of model residuals
#    ri = model.get_residuals()
#    # this timeseries will be NaN except for timesteps with valid observations
#    index = np.nonzero(~np.isnan(ri))[0]
#    # replace the NaNs with zeros
#    ri[np.isnan(ri)] = 0
#
#    # Calculate the jacobian of the residuals with respect to parameters.
#    # Note that the parameters are assumed to be log-transformed `segment_list`s
#    # We pass in the index of valid observations to save time
#    J = model.get_jacobian(index=index)
#
#    Nt, Np = J.shape
#
#    # Gauss-Newton Algorithm
#
#    # Calculate the gradient of the mean square error
#    gradient = 2 * np.dot(J.T, ri)
#    # This should be unnecessary:
#    gradient[np.isnan(gradient)] = 0
#
#    # Find the parameters that have non-zero gradients
#    p_not_flat_gradient = gradient != 0
#
#    # Calculate an approximate Hessian matrix
#    H = 2 * np.dot(J[:, p_not_flat_gradient].T, J[:, p_not_flat_gradient])
#
#    # Calculate the optimum step
#    delta = np.zeros(Np)
#    delta[p_not_flat_gradient] = -np.dot(inv(H + LM_lambda * np.eye(len(H))), gradient[p_not_flat_gradient])
#
#    # Calculate some additional data that can be used for plotting
#    # The gradient term applied to each individual residual:
#    gradient_ip = np.zeros((Nt, Np))
#    gradient_ip[:, p_not_flat_gradient] = 2 * (J[:, p_not_flat_gradient].T * ri).T
#    # The Hessian applied to each individual residual-gradient
#    delta_ip = np.zeros((Nt, Np))
#    delta_ip[:, p_not_flat_gradient] = -np.dot(inv(H), gradient_ip[:, p_not_flat_gradient].T).T
#
#    return delta, {
#        'H': H,
#        'delta': delta,
#        'LM_lambda': LM_lambda,
#        'delta_ip': delta_ip,
#        'gradient': gradient,
#        'gradient_ip': gradient_ip,
#    }
#
#
# def apply_step(model, alpha_step=1, max_delta=None, LM_lambda=0, include_C_old=True, **kwargs):
#    """
#    Applies a delta to a model, returning a new model with modified parameters
#
#    :param model: a sas model
#    :param delta: an array of deltas (changes to log-transformed `segment_list`s)
#    :param max_delta:
#    :param alpha_step:
#    :param include_C_old: if True (default), the value of C_old will be modified. Otherwise it is left as-is
#    :return: a sas model with modified parameters
#    """
#
#    # Calculate how big the step should be
#    delta, gn_info = calc_gauss_newton_step(model, LM_lambda)
#
#    # Apply some limitations on the step size that improve convergence
#    # The fist reduces all steps by a factor
#    delta = alpha_step * delta
#
#    # The second limits the maximum size of step in any direction
#    if max_delta is not None:
#        delta = np.where(np.abs(delta) > max_delta, max_delta * np.sign(delta), delta)
#
#    # Get the old segment list
#    segment_list_old = model.get_segment_list()
#
#    # Apply the delta
#    segment_list_new = segment_list_old * np.exp(delta[:-1])
#
#    # Make a copy of the model without the results
#    new_model = model.copy_without_results()
#
#    # apply the new segment list to the model
#    new_model.update_from_segment_list(segment_list_new)
#
#    # optionally apply the delta to the C_old parameter
#    if include_C_old:
#        for isol, sol in enumerate(model._solorder):
#            old_C_old = model.solute_parameters[sol]['C_old']
#            new_C_old = old_C_old + delta[-model._numsol + isol]
#            new_model.solute_parameters[sol]['C_old'] = new_C_old
#
#    # Run the new model
#    new_model.run()
#    new_mse = np.mean(model.get_residuals() ** 2)
#
#    return new_model, new_mse, gn_info
