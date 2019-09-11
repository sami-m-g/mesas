import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


class Piecewise:

    def __init__(self, segment_list=None, P=None, nsegment=1, ST_max=1., ST_min=0.):
        self.ST_min = np.float(ST_min)
        self.ST_max = np.float(ST_max)
        self.ST_bound = 10 ** 16
        if segment_list is None:
            self.nsegment = nsegment
            self.nparam = self.nsegment + 1
            if P is None:
                self.P = np.r_[0, np.linspace(0, 1, self.nsegment + 1, endpoint=True)]
            self.ST = np.zeros_like(self.P)
            ST_scaled = np.zeros(self.nparam)
            ST_scaled[-1] = 1.
            for i in range(self.nsegment - 1):
                ST_scaled[self.nsegment - i - 1] = np.random.uniform(
                    0, ST_scaled[self.nsegment - i], 1)
            self.ST[1:] = self.ST_min + (self.ST_max - self.ST_min) * ST_scaled
            self.segment_list = np.diff(self.ST)
        else:
            self.nparam = len(segment_list)
            self.nsegment = self.nparam - 1
            if P is None:
                self.P = np.r_[0, np.linspace(0, 1, self.nsegment + 1, endpoint=True)]
            else:
                self.P = np.r_[P]
            self.segment_list = segment_list

    @property
    def segment_list(self):
        return self._segment_list

    @segment_list.setter
    def segment_list(self, new_segment_list):
        self._segment_list = new_segment_list
        self.ST = self._calc_ST()
        self.ST[self.ST > self.ST_bound] = self.ST_bound
        self.ST_min = self.ST[1]
        self.ST_max = self.ST[-1]
        self.interp1d_inv = interp1d(self.P, self.ST,
                                     kind='linear', copy=False,
                                     bounds_error=True, assume_sorted=True)
        self.interp1d = interp1d(self.ST, self.P,
                                 kind='linear', copy=False,
                                 bounds_error=True, assume_sorted=True)

    def _calc_ST(self):
        return np.r_[[0], np.cumsum(self._segment_list)]

    def subdivided_copy(self, segment):
        assert segment < self.nsegment
        P1 = self.P[segment + 1]
        P2 = self.P[segment + 2]
        ST1 = self.ST[segment + 1]
        ST2 = self.ST[segment + 2]
        segment_list = self.segment_list
        P_new = P1 + (P2 - P1) / 2
        segment_list_new = (ST2 - ST1) / 2
        segment_list[segment + 1] = segment_list_new
        P = np.insert(self.P, segment + 2, P_new)
        segment_list = np.insert(segment_list, segment + 1, segment_list_new)
        return Piecewise(segment_list, P)

    def inv(self, P):
        ST = np.where(P > 0,
                      np.where(P < 1, self.interp1d_inv(P), self.ST_max), self.ST[0])
        return ST

    def __call__(self, ST):
        ST_ravel = ST.ravel()
        return self.interp1d(np.where(ST_ravel <= self.ST_max,
                                      np.where(ST_ravel >= self.ST_min, ST_ravel, 0.), 1.)
                             ).reshape(ST.shape)

    def __repr__(self):
        repr = ''
        repr += 'ST        P' + '\n'
        repr += '--------  --------' + '\n'
        for i in range(self.nsegment + 2):
            repr += '{ST:<8.4}  {P:<8.7} \n'.format(ST=self.ST[i], P=self.P[i])
        return repr

    def plot(self, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        return ax.plot(self.ST, self.P, **kwargs)

    def get_jacobian(self, ST, MS, C_old, alpha=1, index=None):
        if index is None:
            index = np.arange(ST.shape[1] - 1)
        Nt = ST.shape[1] - 1
        Ns = self.nsegment
        Omegaj = self.P[1:]
        dOmegaj = np.diff(Omegaj)
        Sj = self.ST[1:]
        dSj = np.diff(Sj)
        if any(dSj == 0):
            pass
        omegaj = dOmegaj / dSj
        J_S = np.zeros((Nt, Ns + 1))
        S_maxcalc = np.diag(ST)
        M_maxcalc = np.diag(MS)
        dST = np.diff(ST, axis=0)
        dMS = np.diff(MS, axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            CS = np.where(dST > 0, dMS / dST, 0)
        CS[np.isnan(CS)] = 0
        for start_index in index:
            for stepend in [0, 1]:
                time_index = start_index + stepend
                Mj = np.interp(Sj, ST[:, time_index], MS[:, time_index])
                dMj = np.diff(Mj)
                Csj = dMj / dSj
                if S_maxcalc[time_index] >= Sj[-1]:
                    seg_maxcalc = Ns
                else:
                    seg_maxcalc = np.argmax(Sj >= S_maxcalc[time_index]) - 1
                # Todo: when Sj=ST, check J_S in both directions and choose steepest
                Tjl = np.digitize(Sj[:seg_maxcalc + 1], ST[:time_index + 1, time_index], right=True) - 1
                Tjr = np.digitize(Sj[:seg_maxcalc + 1], ST[:time_index + 1, time_index], right=False) - 1
                Tjl[Tjl < 0] = 0
                Cpj = (CS[Tjl, time_index] + CS[Tjr, time_index]) / 2.
                if seg_maxcalc > 0:
                    J_S[start_index, 0] += alpha * (Csj[0] - Cpj[0]) * omegaj[0]
                    if seg_maxcalc > 1:
                        J_S[start_index, 1:seg_maxcalc] += alpha * (
                                Csj[1:seg_maxcalc] * omegaj[1:seg_maxcalc]
                                - Csj[0:seg_maxcalc - 1] * omegaj[0:seg_maxcalc - 1]
                                - Cpj[1:seg_maxcalc] * (omegaj[1:seg_maxcalc] - omegaj[0:seg_maxcalc - 1])
                        )
                if seg_maxcalc < Ns:
                    if S_maxcalc[time_index] > Sj[0]:
                        J_S[start_index, seg_maxcalc] += omegaj[seg_maxcalc] * (
                                alpha * (
                                (M_maxcalc[time_index] - Mj[seg_maxcalc]) / dSj[seg_maxcalc] - Cpj[seg_maxcalc])
                                + (C_old * (Sj[1 + seg_maxcalc] - S_maxcalc[time_index])) / dSj[seg_maxcalc]
                        )
                        J_S[start_index, seg_maxcalc + 1] += omegaj[seg_maxcalc] * (
                                -alpha * (
                                (M_maxcalc[time_index] - Mj[seg_maxcalc]) / dSj[seg_maxcalc])
                                + (C_old * (S_maxcalc[time_index] - Sj[seg_maxcalc])) / dSj[seg_maxcalc]
                        )
                else:
                    J_S[start_index, Ns] += -alpha * (Csj[Ns - 1] - Cpj[Ns]) * omegaj[Ns - 1]
        J_S = J_S / 2
        A = np.triu(np.ones(Ns + 1), k=0)
        J_seg = np.dot(A, J_S.T).T * np.r_[Sj[0], np.diff(Sj)]
        return J_seg

    def get_residual_parts(self, C_train, ST, MS, alpha=1, index=None):
        if index is None:
            index = np.nonzero(~np.isnan(C_train))[0]
        Nt = len(C_train)
        Ns = self.nsegment
        Omegaj = self.P[1:]
        dOmegaj = np.diff(Omegaj)
        Sj = self.ST[1:]
        dSj = np.diff(Sj)
        r_seg_i = np.zeros((Nt, Ns))
        S_maxcalc = np.diag(ST)
        M_maxcalc = np.diag(MS)
        dST = np.diff(ST, axis=0)
        dMS = np.diff(MS, axis=0)
        CS = dMS / dST
        CS[dST == 0] = 0
        for start_index in index:
            C_train_this = C_train[start_index]
            for step_end in [0, 1]:
                time_index = start_index + step_end
                if S_maxcalc[time_index] > Sj[0]:
                    Mj = np.interp(Sj, ST[:, time_index], MS[:, time_index])
                    dMj = np.diff(Mj)
                    Csj = dMj / dSj
                    if S_maxcalc[time_index] >= Sj[-1]:
                        seg_maxcalc = Ns
                    else:
                        seg_maxcalc = np.argmax(Sj >= S_maxcalc[time_index]) - 1
                    if seg_maxcalc > 0:
                        r_seg_i[start_index, :seg_maxcalc] += (alpha * Csj[:seg_maxcalc] - C_train_this) * dOmegaj[
                                                                                                           :seg_maxcalc]
                    if seg_maxcalc < Ns:
                        Omegam = self(S_maxcalc[time_index])
                        r_seg_i[start_index, seg_maxcalc] += (alpha * (M_maxcalc[time_index] - Mj[seg_maxcalc])
                                                              / (S_maxcalc[time_index] - Sj[seg_maxcalc])
                                                              - C_train_this) \
                                                             * (Omegam - Omegaj[seg_maxcalc])
        r_seg_i = r_seg_i / 2
        return r_seg_i

#    def regularize(self):
#        reg_frac = 0.05
#        n_reg = 10
#        P_frac = self.P[-2] + (1-reg_frac) * (self.P[-1] - self.P[-2])
#        ST_frac = self.ST[-2] + (1-reg_frac) * (self.ST[-1] - self.ST[-2])
#        P_temp = np.linspace(P_frac, 1, n_reg+1, endpoint=True)
#        self.P_reg = np.concatenate((self.P[: -1], P_temp), axis=None)
#        self.ST_reg = np.zeros_like(self.P_reg)
#        self.ST_reg[: self.nsegment-1] = self.ST[: -1]
#        self.ST_reg[self.nsegment-1:] = ST_frac - (self.ST[-1]-ST_frac) * np.log(
#            (1-self.P_reg[self.nsegment-1:])/(1-P_frac))
#
