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
    """
    Fit the sas function to the data using a least-squares regression optimization

    :param model:
    :param include_C_old:
    :param learn_plot_fun:
    :param index:
    :param kwargs:
    :return:
    """

    if index is None:
        index = model.get_obs_index()

    new_model = model.copy_without_results()

    # Use the current values as the initial estimate
    segment_list = model.get_segment_list()
    # Assume that any segments of zero length are the first segment
    # and exclude these from the estimate of x0
    # TODO: check that zero_segments are indeed ST_min
    non_zero_segments = segment_list > 0
    x0 = np.log(segment_list[non_zero_segments])
    nseg_x0 = len(x0)
    if include_C_old:
        x0 = np.r_[x0, [model.solute_parameters[sol]['C_old'] for sol in model._solorder]]

    # Construct an index into the the columns returned by get_jacobian() that we wish to keep
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
        """return residuals given parameters x"""

        update_parameters(x)
        new_model.run()

        # optionally do some plotting
        if learn_plot_fun is not None:
            learn_plot_fun(new_model)

        return new_model.get_residuals()[index]

    def jac(x):
        """return the jacobian given parameters x"""

        update_parameters(x)
        new_model.run()

        return new_model.get_jacobian(index=index)[:, keep_jac_columns]

    # use the scipy.optimize.least_square package
    OptimizeResult = least_squares(fun=f,
                                   x0=x0,
                                   jac=jac,
                                   verbose=2)

    xopt = OptimizeResult.x
    update_parameters(xopt)
    new_mse = np.mean(OptimizeResult.fun ** 2)

    return new_model, new_mse

