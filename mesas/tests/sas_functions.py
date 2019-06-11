import numpy as np
from scipy.interpolate import interp1d


class Piecewise:
    def __init__(self, params=None, npiece=1, fix_ST_max=False, ST_max=1., fix_ST_min=True, ST_min=0.):
        self.fix_ST_min = fix_ST_min
        self.fix_ST_max = fix_ST_max
        self.ST_min = ST_min
        self.ST_max = ST_max
        if params is None:
            self.npiece = npiece
            self.P = np.r_[0, np.linspace(0, 1, self.npiece+1, endpoint=True)]
            self.ST = np.zeros_like(self.P)
            if self.fix_ST_min:
                self.nparams = self.npiece
            else:
                self.nparams = self.npiece+1
            ST_scaled = np.zeros(self.nparams+1)
            ST_scaled[-1] = 1.
            for i in range(self.nparams-1):
                ST_scaled[self.nparams - i - 1] = np.random.uniform(
                    0, ST_scaled[self.nparams - i], 1)
            if self.fix_ST_min:
                self.ST[1:] = self.ST_min+(self.ST_max-self.ST_min)*ST_scaled
                self.params = np.diff(self.ST)[1:]
            else:
                self.ST[0:] = ST_max*ST_scaled
                self.params = np.diff(self.ST)
        else:
            self.npiece = len(params)
            self.P = np.r_[0, np.linspace(0, 1, self.npiece+1, endpoint=True)]
            self.params = params

    @property
    def params(self):
        return self._params

    @params.setter
    def params(self, new_params):
        self._params = new_params
        self.ST = self._calc_ST()
        self.ST_min = self.ST[1]
        self.ST_max = self.ST[-1]
        self.interp1d_inv = interp1d(self.P, self.ST,
                                     kind='linear', copy=False,
                                     bounds_error=True, assume_sorted=True)
        self.interp1d = interp1d(self.ST, self.P,
                                 kind='linear', copy=False,
                                 bounds_error=True, assume_sorted=True)

    def _calc_ST(self):
        if self.fix_ST_min:
            return np.r_[[0, self.ST_min], np.cumsum(self._params)]
        else:
            return np.r_[[0], np.cumsum(self._params)]

    def inv(self, P):
        ST = np.where(P > 0,
                      np.where(P < 1, self.interp1d_inv(P), self.ST_max), self.ST[0])
        return ST

    def __call__(self, ST):
        return self.interp1d(np.where(ST <= self.ST_max,
                                      np.where(ST >= self.ST_min, ST, 0.), 1.))

    def __repr__(self):
        repr = ''
        repr += 'ST        P'+'\n'
        repr += '--------  --------'+'\n'
        for i in range(self.npiece+2):
            repr += '{ST:<8.4}  {P:<8.7} \n'.format(ST=self.ST[i], P=self.P[i])
        return repr

#    def regularize(self):
#        reg_frac = 0.05
#        n_reg = 10
#        P_frac = self.P[-2] + (1-reg_frac) * (self.P[-1] - self.P[-2])
#        ST_frac = self.ST[-2] + (1-reg_frac) * (self.ST[-1] - self.ST[-2])
#        P_temp = np.linspace(P_frac, 1, n_reg+1, endpoint=True)
#        self.P_reg = np.concatenate((self.P[: -1], P_temp), axis=None)
#        self.ST_reg = np.zeros_like(self.P_reg)
#        self.ST_reg[: self.npiece-1] = self.ST[: -1]
#        self.ST_reg[self.npiece-1:] = ST_frac - (self.ST[-1]-ST_frac) * np.log(
#            (1-self.P_reg[self.npiece-1:])/(1-P_frac))
#
