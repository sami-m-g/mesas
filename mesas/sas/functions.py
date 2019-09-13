"""
    Module functions
    ================
    This module defines classes representing SAS functions.
    Currently there is only one class specified: `Piecewise`. However this
    is a very flexible way of specifying a function.

    Calling an instance of the class supplied with values of ST (as an array, list, or number) evaluates
    the CDF, and returns corresponding cumulative probabilities. The `inv` method takes probabilities and
    returns values of ST.

    See the class docstring for more information.

"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


class Piecewise:
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

    def __init__(self, segment_list=None, ST=None, P=None, nsegment=1, ST_max=1., ST_min=0.):
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
            assert ST[0] == 0
            assert len(ST) > 2
            self._ST = np.array(ST, dtype=np.float)
            self.nsegment = len(ST) - 2
            self._segment_list = self._convert_ST_to_segment_list(self._ST)

        else:

            # Generate a random function
            # Note this should probably be encapsulated in a function
            # and called using an optional keyword, rather than triggered by default

            self.nsegment = nsegment
            self._ST = np.zeros(self.nsegment + 2)

            ST_scaled = np.zeros(self.nsegment + 1)
            ST_scaled[-1] = 1.
            for i in range(self.nsegment - 1):
                ST_scaled[self.nsegment - i - 1] = np.random.uniform(
                    0, ST_scaled[self.nsegment - i], 1)

            self._ST[1:] = self.ST_min + (self.ST_max - self.ST_min) * ST_scaled
            self._segment_list = self._convert_ST_to_segment_list(self._ST)

        # Make a list of probabilities P
        if P is not None:

            # Use the supplied values
            assert P[0] == P[1] == 0
            assert P[-1] == 1
            assert len(P) == len(self._ST)
            self._P = np.r_[P]

        else:

            # Make P equally-spaced
            self._P = np.r_[0, np.linspace(0, 1, self.nsegment + 1, endpoint=True)]

        # call this to make the interpolation functions
        self._make_interpolators()

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

    # Next we define a number of properties. These act like attributes, but special functions
    # are called when the attributes are queried or assigned to

    @property
    def ST(self):
        return self._ST

    @ST.setter
    def ST(self, new_ST):
        assert new_ST[0] == 0
        assert len(new_ST) > 2
        self._ST = new_ST
        self.ST_min = self._ST[1]
        self.ST_max = self._ST[-1]
        self._segment_list = self._convert_ST_to_segment_list(self._ST)
        self._make_interpolators()

    @property
    def P(self):
        return self._P

    @P.setter
    def P(self, new_P):
        assert new_P[0] == new_P[1] == 0
        assert new_P[-1] == 1
        assert len(new_P) == len(self._ST)
        self._P = new_P
        self._make_interpolators()

    @property
    def segment_list(self):
        return self._segment_list

    @segment_list.setter
    def segment_list(self, new_segment_list):
        self._segment_list = new_segment_list
        self.ST = self._convert_segment_list_to_ST(self._segment_list)

    def _convert_segment_list_to_ST(self, seglst):
        """converts a segment_list array into an array of ST endpoints"""
        return np.r_[[0], np.cumsum(seglst)]

    def _convert_ST_to_segment_list(self, ST):
        """converts a segment_list array into an array of ST endpoints"""
        return np.r_[self._ST[1], np.diff(ST[1:])]

    def subdivided_copy(self, segment, split_frac=0.5):
        """Returns a new Piecewise object instance with one segment divided into two

        :param segment: The segment to be split
        :param split_frac: (Optional) The fraction of the segment where the split is made. Default is 0.5
        """
        assert segment < self.nsegment
        P1 = self.P[segment + 1]
        P2 = self.P[segment + 2]
        ST1 = self.ST[segment + 1]
        ST2 = self.ST[segment + 2]
        segment_list = self.segment_list
        P_new = P1 + split_frac * (P2 - P1)
        ST_new = ST1 + split_frac * (ST2 - ST1)
        P = np.insert(self.P, segment + 2, P_new)
        ST = np.insert(self.ST, segment + 2, ST_new)
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
        repr = ''
        repr += 'ST        P' + '\n'
        repr += '--------  --------' + '\n'
        for i in range(self.nsegment + 2):
            repr += '{ST:<8.4}  {P:<8.7} \n'.format(ST=self.ST[i], P=self.P[i])
        return repr

    def plot(self, ax=None, **kwargs):
        """Return a repr of the SAS function"""
        if ax is None:
            fig, ax = plt.subplots()
        return ax.plot(self.ST, self.P, **kwargs)

    def get_jacobian(self, ST, MS, C_old, alpha=1, index=None, mode='segment', logtransform=True):
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

        :param ST: Age-ranked storage array produced by a sas model
        :param MS: Age-ranked solute mass array produced by a sas model
        :param C_old: Old water concentration
        :param alpha: solute partitioning coefficient
        :param index: index of timesteps to calculate the jacobian for
        :param mode: If 'endpoint' return the derivatives with respect to the ST of the endpoints,
        otherwise if 'segment' return the derivatives wrt the segment lengths
        :param logtransform: If True, return the derivatives with respect to the log transform of the segment lengths.
         Ignored if `mode='endpoint'. Default is True.
        :return: numpy array of the Jacobian
        """

        if index is None:
            index = np.arange(ST.shape[1] - 1)

        # Create an array to store the result
        Nt = ST.shape[1] - 1
        Ns = self.nsegment
        J_S = np.zeros((Nt, Ns + 1))

        # pull out the sas function parameters, and get some derivatives and differences
        Omegaj = self.P[1:]
        Sj = self.ST[1:]
        dOmegaj = np.diff(Omegaj)
        dSj = np.diff(Sj)
        omegaj = dOmegaj / dSj

        # The diagonals of these represent the volume and solute mass in the system whose age is known.
        # Everything older is assumed to have concentration C_old
        S_maxcalc = np.diag(ST)
        M_maxcalc = np.diag(MS)

        # Get some differences along the age dimension
        dST = np.diff(ST, axis=0)
        dMS = np.diff(MS, axis=0)

        # This gives the concentration of water of each age
        with np.errstate(divide='ignore', invalid='ignore'):
            CS = np.where(dST > 0, dMS / dST, 0)
        CS[np.isnan(CS)] = 0

        # Loop over the times
        for start_index in index:

            # We are going to do this twice: once for the start and once for the end of each timestep
            # we will average them together at the end
            for stepend in [0, 1]:
                time_index = start_index + stepend

                # Check if we know the concentration in all segments
                if S_maxcalc[time_index] >= Sj[-1]:

                    # If we do set this to the number of segments
                    seg_maxcalc = Ns

                else:

                    # If we don't, set this to the segment that contains the maximum known ST
                    seg_maxcalc = np.argmax(Sj >= S_maxcalc[time_index]) - 1

                # Find the concentration at each segment endpoint
                # Get the age of water at the endpoint
                # NOTE: the two ages Tjl and Tjr are usually the same, but differ only for the case where
                # Sj is exactly equal to ST. In that case, the derivative is actually undefined due to the
                # step-change in concentration. Here we just take the average.
                # Todo: implement upwinding: when Sj=ST, check J_S in both directions and choose steepest
                Tjl = np.digitize(Sj[:seg_maxcalc + 1], ST[:time_index + 1, time_index], right=True) - 1
                Tjr = np.digitize(Sj[:seg_maxcalc + 1], ST[:time_index + 1, time_index], right=False) - 1
                Tjl[Tjl < 0] = 0
                Cpj = (CS[Tjl, time_index] + CS[Tjr, time_index]) / 2.

                # Calculate the bulk-averaged concentration in each segment
                # Age-ranked mass at each segment endpoint
                Mj = np.interp(Sj, ST[:, time_index], MS[:, time_index])
                # Increment of age-ranked mass in each segment
                dMj = np.diff(Mj)
                # Concentration in each segment
                Csj = dMj / dSj

                # Is at least one segment completely full of water of known age?
                if seg_maxcalc > 0:

                    # derivative w.r.t. ST_min
                    J_S[start_index, 0] += alpha * (Csj[0] - Cpj[0]) * omegaj[0]

                    # derivative w.r.t. the endpoints that have completely full segments on both sides
                    if seg_maxcalc > 1:
                        J_S[start_index, 1:seg_maxcalc] += alpha * (
                                Csj[1:seg_maxcalc] * omegaj[1:seg_maxcalc]
                                - Csj[0:seg_maxcalc - 1] * omegaj[0:seg_maxcalc - 1]
                                - Cpj[1:seg_maxcalc] * (omegaj[1:seg_maxcalc] - omegaj[0:seg_maxcalc - 1])
                        )
                # Is the maximum known ST less than ST_max?
                if seg_maxcalc < Ns:

                    # Is it greater than ST_min?
                    if S_maxcalc[time_index] > Sj[0]:
                        # Get the derivative w.r.t the left end of the segment
                        J_S[start_index, seg_maxcalc] += omegaj[seg_maxcalc] * (
                                alpha * (
                                (M_maxcalc[time_index] - Mj[seg_maxcalc]) / dSj[seg_maxcalc] - Cpj[seg_maxcalc])
                                + (C_old * (Sj[1 + seg_maxcalc] - S_maxcalc[time_index])) / dSj[seg_maxcalc]
                        )
                        # Get the derivative w.r.t the right end of the segment
                        J_S[start_index, seg_maxcalc + 1] += omegaj[seg_maxcalc] * (
                                -alpha * (
                                (M_maxcalc[time_index] - Mj[seg_maxcalc]) / dSj[seg_maxcalc])
                                + (C_old * (S_maxcalc[time_index] - Sj[seg_maxcalc])) / dSj[seg_maxcalc]
                        )

                else:

                    # Effect of varying ST_max
                    J_S[start_index, Ns] += -alpha * (Csj[Ns - 1] - Cpj[Ns]) * omegaj[Ns - 1]

        # Average the start and end of the timestep
        J_S = J_S / 2

        if mode == 'endpoint':
            return J_S

        # To get the derivative with respect to the segment length, we add up the derivative w.r.t. the
        # endpoints that would be displaced by varying that segment
        A = np.triu(np.ones(Ns + 1), k=0)
        J_seg = np.dot(A, J_S.T).T

        if logtransform:
            J_seglog = J_seg * np.r_[Sj[0], np.diff(Sj)]
            return J_seglog
        else:
            return J_seg

    def get_residual_parts(self, C_train, ST, MS, alpha=1, index=None):
        """
        Difference between the concentration drawn from each segment and the observed discharge concentration

        :param C_train:
        :param ST: Age-ranked storage array produced by a sas model
        :param MS: Age-ranked solute mass array produced by a sas model
        :param alpha: solute partitioning coefficient
        :param index: index of timesteps to calculate the residuals for
        :return: timeseries of the residuals
        """

        if index is None:
            index = np.arange(ST.shape[1] - 1)

        # Create an array to store the result
        Nt = ST.shape[1] - 1
        Ns = self.nsegment
        r_seg_i = np.zeros((Nt, Ns))

        # pull out the sas function parameters, and get some derivatives and differences
        Omegaj = self.P[1:]
        Sj = self.ST[1:]
        dOmegaj = np.diff(Omegaj)
        dSj = np.diff(Sj)
        omegaj = dOmegaj / dSj

        # The diagonals of these represent the volume and solute mass in the system whose age is known.
        # Everything older is assumed to have concentration C_old
        S_maxcalc = np.diag(ST)
        M_maxcalc = np.diag(MS)

        # Get some differences along the age dimension
        dST = np.diff(ST, axis=0)
        dMS = np.diff(MS, axis=0)

        # This gives the concentration of water of each age
        with np.errstate(divide='ignore', invalid='ignore'):
            CS = np.where(dST > 0, dMS / dST, 0)
        CS[np.isnan(CS)] = 0

        # Loop over the times
        for start_index in index:

            # Pull out the observed concentration we want to compare with
            C_train_this = C_train[start_index]

            # We are going to do this twice: once for the start and once for the end of each timestep
            # we will average them together at the end
            for stepend in [0, 1]:
                time_index = start_index + stepend

                if S_maxcalc[time_index] > Sj[0]:

                    # Check if we know the concentration in all segments
                    if S_maxcalc[time_index] >= Sj[-1]:

                        # If we do set this to the number of segments
                        seg_maxcalc = Ns

                    else:

                        # If we don't, set this to the segment that contains the maximum known ST
                        seg_maxcalc = np.argmax(Sj >= S_maxcalc[time_index]) - 1

                    # Calculate the bulk-averaged concentration in each segment
                    # Age-ranked mass at each segment endpoint
                    Mj = np.interp(Sj, ST[:, time_index], MS[:, time_index])
                    # Increment of age-ranked mass in each segment
                    dMj = np.diff(Mj)
                    # Concentration in each segment
                    Csj = dMj / dSj

                    # Is at least one segment completely full of water of known age?
                    if seg_maxcalc > 0:
                        r_seg_i[start_index, :seg_maxcalc] += (alpha * Csj[:seg_maxcalc] - C_train_this) * dOmegaj[
                                                                                                           :seg_maxcalc]
                    # Is the maximum known ST less than ST_max?
                    if seg_maxcalc < Ns:
                        Omegam = self(S_maxcalc[time_index])
                        r_seg_i[start_index, seg_maxcalc] += (alpha * (M_maxcalc[time_index] - Mj[seg_maxcalc])
                                                              / (S_maxcalc[time_index] - Sj[seg_maxcalc])
                                                              - C_train_this) \
                                                             * (Omegam - Omegaj[seg_maxcalc])

        # Average the start and end of the timestep
        r_seg_i = r_seg_i / 2

        return r_seg_i

