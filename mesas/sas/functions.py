"""
    Module functions
    ================
    This module defines classes representing SAS functions.
    Currently there is only one class specified: :class:`Piecewise`. However this
    is a very flexible way of specifying a function.

    Calling an instance of the class supplied with values of ST (as an array, list, or number) evaluates
    the CDF, and returns corresponding cumulative probabilities. The :func:`Piecewise.inv` method takes probabilities and
    returns values of ST.

    See the class docstring for more information..

"""
import copy

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import rv_continuous


class _SASFunctionBase:
    """
    Base function for SAS functions
    """

    def __init__(self):
        pass

    def __getitem__(self, i):
        return self

    def __setitem__(self, i):
        raise TypeError("You're not allowed to modify the SAS function for a specific index")

    # Next we define a number of properties. These act like attributes, but special functions
    # are called when the attributes are queried or assigned to

    def set_args(self, *args, **kwargs):
        raise NotImplementedError("Method for setting SAS function parameters has not been defined")

    @property
    def has_params(self):
        return self._has_params

    @property
    def ST(self):
        return self._ST

    @ST.setter
    def ST(self, new_ST):
        raise NotImplementedError("Method for setting ST directly has not been defined")

    @property
    def P(self):
        return self._P

    @P.setter
    def P(self, new_P):
        raise NotImplementedError("Method for setting P directly has not been defined")

    @property
    def segment_list(self):
        return self._segment_list

    @segment_list.setter
    def segment_list(self, new_segment_list):
        raise NotImplementedError("Method for setting segments directly has not been defined")

    def _convert_segment_list_to_ST(self, seglst):
        """converts a segment_list array into an array of ST endpoints"""
        return np.cumsum(seglst)

    def _convert_ST_to_segment_list(self, ST):
        """converts a segment_list array into an array of ST endpoints"""
        return np.r_[ST[0], np.diff(ST)]

    def inv(self, P):
        """
        The inverse Cumulative Distribution Function of the SAS function
        """
        raise NotImplementedError("Inverse CDF method not been defined")

    def __call__(self, ST):
        """
        The Cumulative Distribution Function of the SAS function
        """
        raise NotImplementedError("CDF method not been defined")

    def __repr__(self):
        """Return a repr of the SAS function"""
        raise NotImplementedError("__repr__ not been defined")

    def plot(self, ax=None, **kwargs):
        """Return a plot of the SAS function"""
        if ax is None:
            fig, ax = plt.subplots()
        return ax.plot(self.ST, self.P, **kwargs)

    def get_jacobian(self, dCdSj, index=None, mode='segment', logtransform=True):
        """ Calculates a limited jacobian of the model predictions """
        raise NotImplementedError("get_jacobian not been defined")


class Continuous(_SASFunctionBase):
    """
    Base function for SAS functions
    """

    def __init__(self, func, P=None, nsegment=25, ST_max=1.797693134862315e+308):

        self.ST_max = np.float(ST_max)

        assert isinstance(func, rv_continuous)
        self._func = func
        self._has_params = False

        # Make a list of probabilities P
        if P is not None:
            # Use the supplied values
            self.P = np.r_[P]
            self.nsegment = len(P) - 1
        else:
            # Make P equally-spaced
            self.nsegment = nsegment
            self.P = np.linspace(0, 1, self.nsegment + 1, endpoint=True)

    def set_args(self, *args, **kwargs):
        self._frozen_func = self._func(*args, **kwargs)
        self._has_params = True
        self._ST = self.func.ppf(self._P)
        return self

    @property
    def func(self):
        if self._has_params:
            return self._frozen_func
        else:
            raise UnboundLocalError("Parameter values must be set before you can call the underlying function")

    @func.setter
    def func(selfs, new_func):
        raise ValueError("Function type cannot be changed")

    @_SASFunctionBase.ST.setter
    def ST(self, new_ST):
        try:
            assert new_ST[0] >= 0
            assert np.all(np.diff(new_ST) > 0)
            assert len(new_ST) > 1
        except Exception as ex:
            print('Problem with new ST')
            print(f'Attempting to set ST = {new_ST}')
            raise ex
        if new_ST[-1] > self.ST_max:
            self._ST = new_ST
            self.ST_max = self._ST[-1]
        else:
            self._ST = np.r_[new_ST, self.ST_max]
        self._segment_list = self._convert_ST_to_segment_list(self._ST)

    @_SASFunctionBase.P.setter
    def P(self, new_P):
        assert new_P[0] == 0
        assert new_P[-1] == 1
        assert np.all(np.diff(new_P) >= 0)
        self._P = new_P
        if self._has_params:
            self._ST = self.func.ppf(self._P)

    def inv(self, P):
        """
        The inverse Cumulative Distribution Function of the SAS function
        """
        return self.func.ppf(P)

    def __call__(self, ST):
        """
        The Cumulative Distribution Function of the SAS function
        """
        return self.func.cdf(ST)

    def __repr__(self):
        """Return a repr of the SAS function"""
        repr =     f'       {self._func.name} distribution\n'
        if self._has_params:
            repr += f'        parameters: {self._frozen_func.args}\n'
            repr += '        Lookup table version:\n'
            repr += '          ST: '
            for i in range(self.nsegment + 1):
                repr += '{ST:<10.4}  '.format(ST=self.ST[i])
            repr += '\n'
            repr += '          P : '
            for i in range(self.nsegment + 1):
                repr += '{P:<10.4}  '.format(P=self.P[i])
            repr += '\n'
        else:
            repr = f'        (parameters not set)'
        return repr

    def get_jacobian(self, dCdSj, index=None, mode='segment', logtransform=True):
        """ Calculates a limited jacobian of the model predictions """
        raise NotImplementedError("get_jacobian not been defined")


class Piecewise(_SASFunctionBase):
    """
    Class for piecewise-linear SAS functions

    Calling an instance of the class supplied with values of ST (as an array, list, or number) evaluates
    the CDF, and returns corresponding cumulative probabilities. The `inv` method takes probabilities and
    returns values of ST. Additional methods and attributes provide further functionality.

    You can define a piecewise sas function in several different ways. The default is a single segment from 0 to 1:

        >>>sas_fun = Piecewise()
        >>>sas_fun([0.1, 0.3, 0.9])
        array([0.1, 0.3, 0.9])

    The range of the single segment can be modified by providing `ST_min` and `ST_max` keywords

        >>>sas_fun = Piecewise(ST_min=100, ST_max=110)
        >>>sas_fun([99, 101, 103, 109, 111])
        array([0. , 0.1, 0.3, 0.9, 1. ])

    If an integer >1 is passed a `nsegment`, a random distribution will be returned with each of the
    segments representing an equal probability.
    The (ST,P) pairs defining the CDF are given by the `ST` and `P` attributes

        >>>sas_fun.ST
        array([  0.        , 100.        , 103.00391164, 104.17022005,
        110.        ])
        >>>sas_fun.P
        array([0.        , 0.        , 0.33333333, 0.66666667, 1.        ])

    Note that the pair (0,0) is automatically included.
    The segments defining the function are stored in the `segment_list`
    attribute (which is actually a numpy array):

        >>>seglst = sas_fun.segment_list
        array([100.        ,   3.00391164,   1.1663084 ,   5.82977995])

    Note that the first item in the array is ST_min, and each following value is the interval
    along the ST axis to the end of each segment.

    Once created, the values of `segment_list` or `ST` can be modified and the other will update too.
    Values of `P` can also be changed.

    """

    def __init__(self, *args, **kwargs):
        self.set_args(*args, **kwargs)

    def set_args(self, segment_list=None, ST=None, P=None, nsegment=1, ST_max=1., ST_min=0., auto='uniform'):
        """
        Initializes a Piecewise sas function.

        :param segment_list: Optional list of segment lengths in ST. First element is assumed to be ST_min.
        :param ST: Optional list of ST vales at the end of segments. Ignored if `segment_list` is passed
        :param P: Optional list of probabilities at the end of segments. Defaults to evenly spaced intervals
        :param nsegment: If neither `segment_list` nor `ST` are provided, the range between ST_min and ST_max is split into this many randomly, recursively split intervals.
        :param ST_max: The lower bound of the sas function support
        :param ST_min: The upper bound of the sas function support
        """

        self.ST_min = np.float(ST_min)
        self.ST_max = np.float(ST_max)
        self._has_params = False

        # Note that the variables _ST, _P and _segment_list are assigned to below
        # instead of their corresponding properties. This is to avoid triggering the _make_interpolators function
        # until the end

        if segment_list is not None:

            # Use the given segment_list
            self.nsegment = len(segment_list) - 1
            self._segment_list = np.array(segment_list, dtype=np.float)
            self._ST = self._convert_segment_list_to_ST(self._segment_list)

        elif ST is not None:

            # Use the given ST values
            assert ST[0] >= 0
            assert np.all(np.diff(ST) > 0)
            assert len(ST) > 1
            self.nsegment = len(ST) - 1
            self._ST = np.array(ST, dtype=np.float)
            self._segment_list = self._convert_ST_to_segment_list(self._ST)

        else:


            self.nsegment = nsegment

            if auto=='uniform':
                self._ST = np.linspace(self.ST_min, self.ST_max, self.nsegment + 1)
                self._segment_list = self._convert_ST_to_segment_list(self._ST)
            elif auto=='random':
                # Generate a random function
                # Note this should probably be encapsulated in a function
                self._ST = np.zeros(self.nsegment + 1)
                ST_scaled = np.zeros(self.nsegment + 1)
                ST_scaled[-1] = 1.
                for i in range(self.nsegment - 1):
                    ST_scaled[self.nsegment - i - 1] = np.random.uniform(
                        0, ST_scaled[self.nsegment - i], 1)

                self._ST[:] = self.ST_min + (self.ST_max - self.ST_min) * ST_scaled
                self._segment_list = self._convert_ST_to_segment_list(self._ST)

        # Make a list of probabilities P
        if P is not None:

            # Use the supplied values
            assert P[0] == 0
            assert P[-1] == 1
            assert np.all(np.diff(P) >= 0)
            assert len(P) == len(self._ST)
            self._P = np.r_[P]

        else:

            # Make P equally-spaced
            self._P = np.linspace(0, 1, self.nsegment + 1, endpoint=True)

        self._has_params = True
        # call this to make the interpolation functions
        self._make_interpolators()
        return self

    def _make_interpolators(self):
        """
        Create functions that can be called to return the CDF and inverse CDF

        Uses scipy.optimize.interp1d
        """
        self.interp1d_inv = interp1d(self.P, self.ST,
                                     fill_value=(self.ST_min, self.ST_max),
                                     kind='linear', copy=False,
                                     bounds_error=False, assume_sorted=True)
        self.interp1d = interp1d(self.ST, self.P,
                                 fill_value=(0., 1.),
                                 kind='linear', copy=False,
                                 bounds_error=False, assume_sorted=True)

    @_SASFunctionBase.ST.setter
    def ST(self, new_ST):
        try:
            assert new_ST[0] >= 0
            assert np.all(np.diff(new_ST) > 0)
            assert len(new_ST) > 1
        except Exception as ex:
            print('Problem with new ST')
            print(f'Attempting to set ST = {new_ST}')
            if not np.all(np.diff(new_ST) > 0):
                print("   -- if ST values are not distinct, try changing 'ST_largest_segment' and/or 'ST_smallest_segment'")
            raise ex
        self._ST = new_ST
        self.ST_min = self._ST[0]
        self.ST_max = self._ST[-1]
        self._segment_list = self._convert_ST_to_segment_list(self._ST)
        self._make_interpolators()

    @_SASFunctionBase.P.setter
    def P(self, new_P):
        assert new_P[0] == 0
        assert new_P[-1] == 1
        assert len(new_P) == len(self._ST)
        self._P = new_P
        self._make_interpolators()

    @_SASFunctionBase.segment_list.setter
    def segment_list(self, new_segment_list):
        self._segment_list = new_segment_list
        self.ST = self._convert_segment_list_to_ST(self._segment_list)

    def subdivided_copy(self, segment, split_frac=0.5):
        """Returns a new Piecewise object instance with one segment divided into two

        :param segment: The segment to be split
        :param split_frac: (Optional) The fraction of the segment where the split is made. Default is 0.5
        """
        assert segment < self.nsegment
        P1 = self.P[segment]
        P2 = self.P[segment + 1]
        ST1 = self.ST[segment]
        ST2 = self.ST[segment + 1]
        P_new = P1 + split_frac * (P2 - P1)
        ST_new = ST1 + split_frac * (ST2 - ST1)
        P = np.insert(self.P, segment + 1, P_new)
        ST = np.insert(self.ST, segment + 1, ST_new)
        return Piecewise(ST=ST, P=P)

    def inv(self, P):
        """
        The inverse Cumulative Distribution Function of the SAS function

        :param P: Probabilities (can be a number, list or array-like)
        :return: A numpy array of corresponding ST values
        """
        P_arr = np.array(P)
        P_ravel = P_arr.ravel()
        return self.interp1d(P_ravel).reshape(P_arr.shape)

    def __call__(self, ST):
        """
        The Cumulative Distribution Function of the SAS function

        :param ST: Age-ranked Storage values (can be a number, list or array-like)
        :return: A numpy array of corresponding probabilities
        """
        ST_arr = np.array(ST)
        ST_ravel = ST_arr.ravel()
        return self.interp1d(ST_ravel).reshape(ST_arr.shape)

    def __repr__(self):
        """Return a repr of the SAS function"""
        if not hasattr(self, 'nsegment'):
            return 'Not initialized'
        repr = '        ST: '
        for i in range(self.nsegment + 1):
            repr += '{ST:<10.4}  '.format(ST=self.ST[i])
        repr += '\n'
        repr += '        P : '
        for i in range(self.nsegment + 1):
            repr += '{P:<10.4}  '.format(P=self.P[i])
        repr += '\n'
        return repr

    def get_jacobian(self, dCdSj, index=None, mode='segment', logtransform=True):
        """
        Calculates a limited jacobian of the model predictions

        If :math:`\vec{y} = f(\vec{x})` the Jacobian is a matrix where each i, j entry is
        .. math:

            J_{i,j}=\frac{\partial y_i}{\partial x_j}

        Here, the :math:`y_i` values are the predicted concentrations at each timestep,
        and the :math:`x_j` are either the lengths of each piecewise segment along with ST_min (default),
        or the values of ST at the endpoints of each segment (with `mode='endpoint'`).

        The Jacobian is limited as it represents the sensitivity of the model predictions of output
        concentrations at each timestep to variations in the sas function at that timestep, and neglects
        the cumulative effects that changing the sas function would have on the state variables (ST and MS).

        :param dCdSj: jacobian produced by solve.f90
        :param index: index of timesteps to calculate the jacobian for
        :param mode: If 'endpoint' return the derivatives with respect to the ST of the endpoints,
        otherwise if 'segment' return the derivatives wrt the segment lengths
        :param logtransform: If True, return the derivatives with respect to the log transform of the segment lengths.
         Ignored if `mode='endpoint'. Default is True.
        :return: numpy array of the Jacobian
        """

        if index is None:
            index = np.arange(dCdSj.shape[0])

        J_S = dCdSj[index, :]

        if mode == 'endpoint':
            return J_S

        # To get the derivative with respect to the segment length, we add up the derivative w.r.t. the
        # endpoints that would be displaced by varying that segment
        A = np.triu(np.ones(self.nsegment + 1), k=0)
        J_seg = np.dot(A, J_S.T).T

        if logtransform:
            J_seglog = J_seg * self._segment_list
            return J_seglog
        else:
            return J_seg

class PiecewiseASD(Piecewise):
    def __init__(self, STc, QTc, S=None, Q=None):
        self._STc = STc
        self._QTc = QTc
        self._S = None
        self._Q = None
        self._sas_fun_index = None
        ST = STc[-1] - STc[::-1]
        QT = QTc[-1] - QTc[::-1]
        P = QT / QTc[-1]
        super().__init__(ST=ST, P=P)
        self.set_SQ(S, Q)

    def __deepcopy__(self, memo):
        cls = self.__class__
        new_obj = cls.__new__(cls)
        memo[id(self)] = new_obj
        new_obj.__init__(copy.deepcopy(self._STc), copy.deepcopy(self._QTc),
                         copy.deepcopy(self._S), copy.deepcopy(self._Q))
        return new_obj

    def set_SQ(self, S=None, Q=None):
        pass
        # if S is not None:
        #    self._S = S
        # if Q is not None:
        #    self._Q = Q
        # if self._S is not None and self._Q is not None:
        #    self._sas_fun_index = [0] * len(self._S)
        #    for i in range(len(self._S)):
        #        num_orig = ((self._STc < self._S[i]) & (self._QTc < self._Q[i])).sum()
        #        if num_orig > 0:
        #            STc = np.zeros(num_orig + 1)
        #            QTc = np.zeros(num_orig + 1)
        #            STc[:-1] = self._STc[:num_orig]
        #            QTc[:-1] = self._QTc[:num_orig]
        #            STc[-1] = self._S[i]
        #            QTc[-1] = self._Q[i]
        #            index_select = range(1, num_orig+1)
        #            orig_select = range(self.nsegment+1-num_orig, self.nsegment+1)
        #        else:
        #            STc = self._STc[:2]
        #            QTc = np.zeros_like(STc)
        #            orig_select = []
        #            index_select = []
        #        ST = STc[-1] - STc[::-1]
        #        QT = QTc[-1] - QTc[::-1]
        #        P = QT/QT[-1]
        #        self._sas_fun_index[i] = Piecewise(ST=ST, P=P)
        #        self._sas_fun_index[i].orig_select = orig_select
        #        self._sas_fun_index[i].index_select = index_select

    # def __call__(self, *args, **kwargs):
    # warnings.warn('Just FYI, you are calling an ASD function without an index. Is that what you want to do?')
    # return super().__call__(*args, **kwargs)

    # def __getitem__(self, i):
    #    return self._sas_fun_index[i]

    def update_from_Piecewise_SAS(self, sas_fun):
        ST = sas_fun.ST
        P = sas_fun.P
        super().__init__(ST=ST, P=P)
        self._STc = self._STc[-1] - ST[::-1]
        self._QTc = self._QTc[-1] * (1 - P[::-1])
        self.set_SQ()

    @_SASFunctionBase.ST.setter
    def ST(self, new_ST):
        self._STc = self._STc[-1] - new_ST[::-1]
        super(PiecewiseASD, self.__class__).ST.fset(self, new_ST)
        self.set_SQ()

    @_SASFunctionBase.P.setter
    def P(self, new_P):
        self._QTc = self._QTc[-1] * (1 - new_P[::-1])
        super(PiecewiseASD, self.__class__).P.fset(self, new_P)
        self.set_SQ()

    def __repr__(self):
        """Return a repr of the SAS function"""
        repr = '        STc : '
        for i in range(self.nsegment + 1):
            repr += '{STc:<10.4}  '.format(STc=self._STc[i])
        repr += '\n'
        repr += '        QTc : '
        for i in range(self.nsegment + 1):
            repr += '{QTc:<10.4}  '.format(QTc=self._QTc[i])
        repr += '\n'
        return repr

    def plot(self, ax=None, Q_transform=None, ASD=True, **kwargs):
        """Return a repr of the SAS function"""
        if ax is None:
            fig, ax = plt.subplots()
        if ASD:
            if Q_transform is None:
                return ax.plot(self._STc, self._QTc, **kwargs)
            else:
                return ax.plot(self._STc, Q_transform(self._QTc), **kwargs)
        else:
            super().plot(self, ax=ax, **kwargs)
