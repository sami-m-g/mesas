import inspect

import numpy as np
from scipy.optimize import least_squares
from sklearn.model_selection import KFold

VERBOSE = True

ITERATION = None

def _verbose(string, **kwargs):
    """prints progress updates"""
    if VERBOSE:
        print(string, **kwargs)


def run(model, verbose=True, search_mode='leftfirst', **kwargs):
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
    _verbose(model)

    # Create a dictionary that will hold new components (sas functions) as we find them.
    # It will also serve to keep track of which components to keep subdividing (if there is more than 1
    # NOTE: a component is an object that packages together a sas function and a timeseries of its weight
    # along with a label identifying it. It allows us to specify a time-varying SAS function as a
    # weighted sum of fixed SAS functions
    new_components = {}
    for flux in initial_model.components_to_learn.keys():
        new_components[flux] = dict((label, None) for label in initial_model.components_to_learn[flux])

    index = initial_model.get_obs_index()

    # Fit the initial model to the observed data
    _verbose('-- Refining the initial model. ', end='')
    better_initial_model, _ = fit_model(initial_model, index=index, **kwargs)

    better_initial_mse = cross_validation_mse(better_initial_model, index=index, **kwargs)

    kwargs['incres_fun'](better_initial_model, 'initial')

    # run the algorithm
    global ITERATION
    ITERATION= 0
    if search_mode=='scanning':
        return increase_resolution_scanning(better_initial_model, better_initial_mse, new_components, segment=0,
                                            index=index, **kwargs)
    elif search_mode=='leftfirst':
        return increase_resolution_leftfirst(better_initial_model, better_initial_mse, new_components, segment=0,
                                             index=index, **kwargs)
    else:
        print('Never be here')



def lookfor_new_components(initial_model, initial_mse, new_components, segment_dict, id_str, incres_fun=None, index=None,
                           **kwargs):
    # This keeps track of how many new components we found
    new_components_count = 0
    # This will keep track of convergence information provided by fit_model in case we want to use it for plotting
    mse_dict = {}

    # These are retained in case they save us having to run the model again
    last_accepted_model, last_accepted_mse = None, None

    # loop over the components we want to improve, testing whether subdividing each individually leads to
    # substantial change in the fit
    for flux in new_components.keys():
        for label in new_components[flux].keys():

            if type(segment_dict) is dict:
                segment = segment_dict[flux][label]
            else:
                segment = segment_dict

            nsegment = initial_model.sas_blends[flux].components[label].sas_fun.nsegment
            if segment < nsegment:
                _verbose(f'Subdividing {flux} {label} segment {segment}')

                new_model = initial_model.subdivided_copy(flux, label, segment)
                segment_list = new_model.get_segment_list()
                if np.any(segment_list[segment_list > 0] < new_model.options['ST_smallest_segment']):
                    _verbose("New segment too small! Try reducing model.options['ST_smallest_segment']")
                    subdivision_accepted = False
                    new_mse = None
                else:
                    _verbose('-- Initial fit of the candidate model ', end='')
                    new_model, new_mse = fit_model(new_model, index=index, **kwargs)

                    # The current criteria used for deciding whether to accept the subdivision is simple:
                    # Just check whether the model improvement is greater than some threshold
                    # This can definitely be improved
                    # subdivision_accepted = (1 - new_mse / initial_mse) > split_tolerance

                    # Use k-fold cross validation to see if the split leads to an improvement in the model
                    new_mse = cross_validation_mse(new_model, index=index, **kwargs)
                    _verbose('Candidate model = ')
                    _verbose(new_model)
                    subdivision_accepted = new_mse < initial_mse

                    _verbose(f'Initial mse = {initial_mse}')
                    _verbose(f'Candidate mse     = {new_mse}')
                    # Optionally, make some plots
                    if incres_fun is not None:
                        incres_fun(new_model, f'{id_str}_{flux}_{label}_{nsegment}_{segment}_{subdivision_accepted}')

                # accept or reject the subdivision
                if subdivision_accepted:
                    _verbose(f'Subdivision accepted for {flux}, {label} segment {segment}')

                    # Add to the record of new components, and increment the counter
                    new_components[flux][label] = new_model.sas_blends[flux].components[label]
                    new_components_count += 1
                    last_accepted_model, last_accepted_mse = new_model, new_mse

                else:
                    _verbose(f'Subdivision rejected for {flux}, {label}')

                    # Record the rejection by setting the dict entry to None
                    new_components[flux][label] = None

                # Keep a record of the mse for plotting
                mse_dict[f'{flux}, {label}'] = new_mse

    return last_accepted_model, last_accepted_mse, new_components, new_components_count, mse_dict


def incorporate_new_components(initial_model, new_components):
    """
    Adds components into a model from a dictionary

    :param initial_model:
    :param new_components:
    :param index:
    :param kwargs:
    :return:
    """
    new_model = initial_model.copy_without_results()

    for flux in new_components.keys():
        for label in new_components[flux].keys():

            if new_components[flux][label] is not None:
                new_model.set_component(flux, new_components[flux][label])

    return new_model


def increase_resolution_scanning(initial_model, initial_mse, new_components, maxscan=100,
                                 index=None, incres_fun=None,
                                 **kwargs):
    """Increases sas function resolution by splitting a piecewise linear segment (scanning method)

    :param initial_model: The model we want to try increasing the resolution of
    :param initial_mse: How well it currently fits the data
    :param new_components: A dict that keeps track of which components we are trying to refine
    :param segment: The index of the segment we are currently refining
    :param kwargs: Additional arguments to be passed to the fit_model function
    :return:
    """

    _verbose(f'\n***********************\nStarting {inspect.stack()[0][3]}')

    any_segments_subdivided = False

    for scan in range(maxscan):
        _verbose(f'\nStarting scan #{scan + 1}')

        # How many segments do we have to scan?
        max_segment = 0
        check_segment_dict = {}
        for flux in new_components.keys():
            check_segment_dict[flux] = {}
            for label in new_components[flux].keys():
                if max_segment < initial_model.sas_blends[flux].components[label].sas_fun.nsegment:
                    max_segment = initial_model.sas_blends[flux].components[label].sas_fun.nsegment
                # start with segment 0
                check_segment_dict[flux][label] = 0

        for segment in range(max_segment):
            _verbose(f'Testing increased SAS resolution in segment {segment}')

            last_accepted_model, last_accepted_mse, new_components, new_components_count, mse_dict = lookfor_new_components(
                initial_model, initial_mse,
                new_components, check_segment_dict,
                id_str=f'scan_{scan}', incres_fun=incres_fun,
                index=None, **kwargs)

            # if we accepted refined components we want to add them into the model
            # before moving on to the next segment or scan
            if new_components_count > 0:

                any_segments_subdivided = True

                # If we found more than one new component, we need to include all of them, then re-run
                # fit_model so that their interactions can be accounted for
                if new_components_count > 1:
                    _verbose('-- Combining and fitting ', end='')
                    new_model, _ = fit_model(incorporate_new_components(initial_model, new_components),
                                             index=index, **kwargs)
                    new_mse = cross_validation_mse(new_model, index=index, **kwargs)

                else:
                    new_model, new_mse = last_accepted_model, last_accepted_mse

                # Optionally, make some plots
                if incres_fun is not None:
                    incres_fun(new_model, f'scan_{scan}_all_components_{max_segment}_{segment}')

                _verbose('New model is:')
                _verbose(new_model)

                initial_model = new_model
                initial_mse = new_mse

            # Increment the segment counter
            for flux in check_segment_dict.keys():
                for label in check_segment_dict[flux].keys():
                    if new_components[flux][label] is not None:
                        check_segment_dict[flux][label] += 1
                    if check_segment_dict[flux][label] < initial_model.sas_blends[flux].components[
                        label].sas_fun.nsegment:
                        check_segment_dict[flux][label] += 1

        if any_segments_subdivided:
            any_segments_subdivided = False
            _verbose(f'Scan {scan + 1} complete')
        else:
            # If not, we are done. Return the model we've found,
            _verbose('Finished increasing model resolution')
            return new_model

    _verbose(f'Maximum number of scans reached ({maxscan})')
    return new_model


def increase_resolution_leftfirst(initial_model, initial_mse, new_components, segment, incres_fun=None, index=None,
                                  **kwargs):
    """Increases sas function resolution by splitting a piecewise linear segment (left-first method)

    :param initial_model: The model we want to try increasing the resolution of
    :param initial_mse: How well it currently fits the data
    :param new_components: A dict that keeps track of which components we are trying to refine
    :param segment: The index of the segment we are currently refining
    :param incres_fun: An optional function that can produce some plots along the way
    :param kwargs: Additional arguments to be passed to the fit_model function
    :return:
    """

    global ITERATION
    ITERATION +=1

    _verbose(f'\n***********************\nStarting {inspect.stack()[0][3]} iteration {ITERATION}')
    _verbose(f'Testing increased SAS resolution in segment {segment}')

    # total number of segments. Used for plotting (actually for naming the saved figure files)
    Ns = len(initial_model.get_segment_list())

    last_accepted_model, last_accepted_mse, new_components, new_components_count, mse_dict = lookfor_new_components(
        initial_model, initial_mse,
        new_components, segment,
        id_str=f'iteration_{ITERATION}', incres_fun=incres_fun,
        index=None, **kwargs)

    # if we accepted refined components we want to add them into the model
    # and then recursively try subdividing them too
    if new_components_count > 0:

        # If we found more than one new component, we need to include all of them, then re-run
        # fit_model so that their interactions can be accounted for
        if new_components_count > 1:
            _verbose('-- Combining and fitting ', end='')
            new_model, _ = fit_model(incorporate_new_components(initial_model, new_components),
                                     index=index, **kwargs)
            new_mse = cross_validation_mse(new_model, index=index, **kwargs)

        else:
            new_model, new_mse = last_accepted_model, last_accepted_mse

        # Optionally, make some plots
        if incres_fun is not None:
            incres_fun(new_model, f'iteration_{ITERATION}_all_components')

        _verbose('New model is:')
        _verbose(new_model)


        # Recursively try subdividing the first of the two new segments
        return increase_resolution_leftfirst(new_model, new_mse, new_components, segment,
                                             incres_fun=incres_fun, index=index, **kwargs)

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
            return increase_resolution_leftfirst(initial_model, initial_mse, new_components, segment + 1,
                                                 incres_fun=incres_fun, index=index, **kwargs)
        else:

            # If not, we are done. Return the model we've found,
            _verbose('Finished increasing model resolution')
            return initial_model


def fit_model(model, include_C_old=True, learn_fun=None, index=None, jacobian_mode='analytical', **kwargs):
    """
    Fit the sas function to the data using a least-squares regression optimization

    :param model:
    :param include_C_old:
    :param learn_fun:
    :param index:
    :param kwargs:
    :return:
    """

    _verbose('Fitting...', end='')

    if index is None:
        index = model.get_obs_index()

    new_model = model.copy_without_results()

    # Use the current values as the initial estimate
    segment_list = new_model.get_segment_list()
    # Assume that any segments of zero length are the first segment
    # and exclude these from the estimate of x0
    # TODO: check that zero_segments are indeed ST_min
    non_zero_segments = segment_list > 0
    x0 = np.log(segment_list[non_zero_segments])
    xmin = np.ones_like(x0) * np.log(new_model.options['ST_smallest_segment'])
    xmax = np.ones_like(x0) * np.log(new_model.options['ST_largest_segment'])
    if any((x0 < xmin) | (x0 > xmax)):
        print(f'A segment size is out of bounds: {np.exp(x0)}')
    nseg_x0 = len(x0)
    if include_C_old:
        x0 = np.r_[x0, [new_model.solute_parameters[sol]['C_old'] for sol in new_model._solorder]]
        xmin = np.r_[xmin, [-np.inf for sol in new_model._solorder]]
        xmax = np.r_[xmax, [np.inf for sol in new_model._solorder]]

    # Construct an index into the the columns returned by get_jacobian() that we wish to keep
    keep_jac_columns = np.ones(len(segment_list) + new_model._numsol) == 1
    keep_jac_columns[:len(segment_list)] = non_zero_segments

    new_model.x_prev = None

    def update_parameters(x):

        if new_model.x_prev is not None and np.all(x == new_model.x_prev):
            return False

        else:
            #_verbose(x)
            new_segment_list = np.zeros_like(segment_list)
            new_segment_list[non_zero_segments] = np.exp(x[:nseg_x0])
            new_model.update_from_segment_list(new_segment_list)

            # optionally update the C_old parameter
            if include_C_old:
                for isol, sol in enumerate(new_model._solorder):
                    new_model.solute_parameters[sol]['C_old'] = x[-new_model._numsol + isol]
            new_model.x_prev = x
            return True

    def f(x):
        """return residuals given parameters x"""

        update_parameters(x)
        new_model.run()

        # optionally do some plotting
        if learn_fun is not None:
            learn_fun(new_model)

        residuals = new_model.get_residuals()[index]
        return residuals

    def jac(x):
        """return the jacobian given parameters x"""

        if update_parameters(x):
            new_model.run()

        jacobian = new_model.get_jacobian(index=index)[:, keep_jac_columns]
        return jacobian

    # use the scipy.optimize.least_square package
    if jacobian_mode == 'analytical':
        OptimizeResult = least_squares(fun=f,
                                       x0=x0,
                                       jac=jac,
                                       verbose=0,
                                       bounds=(xmin, xmax))
    elif jacobian_mode == 'numerical':
        OptimizeResult = least_squares(fun=f,
                                       x0=x0,
                                       verbose=0,
                                       bounds=(xmin, xmax))

    xopt = OptimizeResult.x
    update_parameters(xopt)
    new_mse = np.mean(OptimizeResult.fun ** 2)
    new_model.mse = new_mse

    _verbose('done.')
    return new_model, new_mse


def cross_validation_mse(model, index=None, **kwargs):
    if index is None:
        index = model.get_obs_index()
    n_splits = 5
    kf = KFold(n_splits=n_splits)
    mse_cv_list = []
    _verbose('Getting mse from k-fold cross validation')
    iter = 0
    for train, test in kf.split(index):
        iter +=1
        _verbose(f'-- iteration {iter} of {n_splits} ', end='')
        new_model_cv, _ = fit_model(model, index=index[train], **kwargs)
        mse_cv_list.append(np.mean(new_model_cv.get_residuals()[index[test]]**2))
    _verbose(f'    mse_cv = {mse_cv_list}')
    #mse_cv_list.remove(max(mse_cv_list))
    #mse_cv_list.remove(min(mse_cv_list))
    mse_cv = np.mean(np.sqrt(mse_cv_list))**2
    model.mse_cv = mse_cv
    return mse_cv
