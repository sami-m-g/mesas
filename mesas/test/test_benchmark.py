#%%
from datetime import time
import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt
import logging

from mesas.sas.model import Model
from scipy.stats import gamma, beta
from scipy.special import lambertw
#%%
logging.basicConfig(filename='test.log', level=logging.INFO)

#%%

"""
"""
def M(delta, i):
    return np.where(i==0, 1, -lambertw(-np.exp(-1 - (i*delta)/2.))*(2 + lambertw(-np.exp(-1 - (i*delta)/2.))))

steady_benchmarks = {
	'Uniform':{
		'spec':{
			'ST': ['S_m', 'S_m0']
        },
        'pQdisc': lambda delta, i: (-1 + np.exp(delta))**2/(np.exp((1 + i)*delta)*delta),
        'pQdisc0':lambda delta: (1 + np.exp(delta)*(-1 + delta))/(np.exp(delta)*delta)
	},
    'Exponential':{
		'spec':{
			"func": "gamma",
			"args": {"a": 1.0, "scale": "S_0", "loc": "S_m"},
        },
        'pQdisc': lambda delta, i: (2*np.log(1 + i*delta) - np.log((1 + (-1 + i)*delta)*(1 + delta + i*delta)))/delta,
        'pQdisc0':lambda delta: (delta + np.log(1/(1 + delta)))/delta
    },
    'Biased old':{
		'spec':{
			"func": "beta",
			"args": {"a": 2.01, "b":1.01, "scale": "S_0", "loc": "S_m"},
        },
        'pQdisc': lambda delta, i: (2*1/np.cosh(delta - i*delta)*1/np.cosh(delta + i*delta)*np.sinh(delta)**2*np.tanh(i*delta))/delta,
        'pQdisc0':lambda delta: 1 - np.tanh(delta)/delta
    },
    'Biased young':{
		'spec':{
			"func": "beta",
			"args": {"a": 1.01, "b":2.01, "scale": "S_0", "loc": "S_m"},
        },
        'pQdisc': lambda delta, i: (2*delta)/((1 + (-1 + i)*delta)*(1 + i*delta)*(1 + delta + i*delta)),
        'pQdisc0':lambda delta: delta/(1 + delta)
    },
    'Partial bypass':{
		'spec':{
			"func": "beta",
			"args": {"a": 1.0/2, "b":1.0, "scale": "S_0", "loc": "S_m"},
        },
        'pQdisc': lambda delta, i: (M(delta, -1 + i) - 2*M(delta, i) + M(delta, 1 + i))/delta,
        'pQdisc0':lambda delta: (-1 + delta + M(delta, 1))/delta
    },
    'Balanced':{
		'spec':{
			"func": "beta",
			"args": {"a": 1.0, "b":1.0/2, "scale": "S_0", "loc": "S_m"},
        },
        'pQdisc': lambda delta, i: delta/2,
        'pQdisc0':lambda delta: delta/4
    }
}
unsteady_benchmarks = {
    'Linear reservoir':{
        'spec':{
            'ST': [0, 'S_0']
        },
        'pQdisc': lambda J, Q, delta, i, j: J[j-i]/Q[j] * (1 - np.exp(-delta)),
        'pQdisc0':lambda J, Q, delta, j: J[j]/Q[j] * (delta - (1 - np.exp(-delta))) / delta
    }
}
#%%
def test_steady(makefigure=True):
#%%
    np.random.seed(2)
    timeseries_length = 100
    max_age = timeseries_length
    dt = 0.1
    N_bm = len(steady_benchmarks)
    Q_0 = 1.0 # <-- steady-state flow rate
    C_J = 1000 + np.random.randn(timeseries_length)*1000.0

    C_old = 2000.0
    J = Q_0
    S_0 = 5 * Q_0
    S_m = 1 * Q_0
    n_substeps = 3
    debug = False
    verbose = False
    jacobian = False

    data_df = pd.DataFrame(index=range(timeseries_length))
    data_df['t'] = data_df.index * dt
    data_df["J"] = J
    data_df["S_0"] = S_0
    data_df["S_m"] = S_m
    data_df["S_m0"] = S_m + S_0
    data_df["C"] = C_J


    if makefigure:
        fig = plt.figure(figsize=[16,9])
        nrow = 2
        ncol = int(len(steady_benchmarks)/nrow)+1
        figcount=0
    for name, bm in steady_benchmarks.items():
        i = np.arange(timeseries_length)
        pQdisc = np.zeros_like(i, dtype=float)
        delta = dt * Q_0 / (S_0)
        Tm = S_m/Q_0
        pQdisc[0] = bm['pQdisc0'](delta) / dt
        pQdisc[1:] = bm['pQdisc'](delta, i[1:]) / dt
        im = int(Tm/dt)
        data_df[f'{name} benchmark C'] = C_old
        data_df[f'{name} benchmark C'][im:] = (np.convolve(C_J, pQdisc, mode='full')[:timeseries_length-im] * dt
        + C_old * (1-np.cumsum(pQdisc)[:timeseries_length-im]*dt))
        data_df[f'{name} benchmark t'] = data_df['t']

        data_df[name] = Q_0
        solute_parameters = {
            "C":{
                "C_old": C_old, 
                "observations": {name:f'{name} benchmark C'},
            }
        }
        sas_specs = {name:{f'{name}_SAS':bm['spec']}}

        model = Model(
            data_df,
            sas_specs,
            solute_parameters,
            debug=debug,
            verbose=verbose,
            dt=dt,
            n_substeps=n_substeps,
            jacobian=jacobian,
            max_age=max_age,
            warning=True
        )
        model.run()
        data_df = model.data_df
        err = (data_df[f'{name} benchmark C'].values-data_df[f'C --> {name}'].values)/data_df[f'{name} benchmark C'].values
        logging.info(f'{name} error = {err.mean()}')
        assert err.mean()<1E-3
        if makefigure:
            icol = int(figcount/nrow)
            irow = figcount - nrow * icol
            ax = plt.subplot2grid((nrow, ncol),(irow, icol))
            ax.plot(data_df['t'], data_df[f'C --> {name}'], alpha=0.3, lw=2)
            ax.plot(data_df[f'{name} benchmark t'], data_df[f'{name} benchmark C'], 'r--', alpha=0.3, lw=2)
            ax.set_title(name)
            ax.set_ylim((500, 2500))
            ax.set_ylabel('Tracer conc.')
            ax.set_xlabel('Time')
            figcount+=1
    if makefigure:
        fig.tight_layout()
        fig.savefig('test_steady.pdf')
# %%

def test_unsteady(makefigure=True):
#%%
    makefigure=True
    np.random.seed(1)
    timeseries_length = 300
    max_age = timeseries_length
    dt = 0.1
    k = 1/5.
    delta = dt * k
    Q_0 = 1.0 # <-- steady-state flow rate
    C_J = 1000 + np.random.randn(timeseries_length)*1000.0
    C_old = 2000.0
    J = Q_0 *np.maximum(Q_0 + np.random.randn(timeseries_length)*Q_0*0.5,0)
    S = np.zeros(timeseries_length+1)
    M = np.zeros(timeseries_length+1)
    S[0] = Q_0 / k
    M[0] = C_old * S[0]
    for j in range(timeseries_length):
        S[j+1] = S[j] * np.exp(-delta) + dt * J[j] * (1 - np.exp(-delta))/delta
        M[j+1] = M[j] * np.exp(-delta) + dt * C_J[j] * J[j] * (1 - np.exp(-delta))/delta
    Q = J - np.diff(S)/dt
    CQ = (C_J * J - np.diff(M)/dt)/Q
    S_0 = (S[:-1] + S[1:])/2
    n_substeps = 5
    debug = False
    verbose = False
    jacobian = False

    data_df = pd.DataFrame(index=range(timeseries_length))
    data_df['t'] = data_df.index * dt
    data_df["J"] = J
    data_df["S_0"] = S_0
    data_df["C"] = C_J


    if makefigure:
        fig = plt.gcf()
    for name, bm in unsteady_benchmarks.items():
        data_df[f'{name} benchmark C'] = CQ
        data_df[f'{name} benchmark t'] = data_df['t']

        data_df[name] = Q_0
        solute_parameters = {
            "C":{
                "C_old": C_old, 
                "observations": {name:f'{name} benchmark C'},
            }
        }
        sas_specs = {name:{f'{name}_SAS':bm['spec']}}

        model = Model(
            data_df,
            sas_specs,
            solute_parameters,
            debug=debug,
            verbose=verbose,
            dt=dt,
            n_substeps=n_substeps,
            jacobian=jacobian,
            max_age=max_age,
            warning=True
        )
        model.run()
        data_df = model.data_df
        err = (data_df[f'{name} benchmark C'].values-data_df[f'C --> {name}'].values)/data_df[f'{name} benchmark C'].values
        logging.info(f'Unsteady {name} error = {err.mean()}')
        assert err.mean()<1E-2
        if makefigure:
            ax = plt.subplot()
            ax.plot(data_df['t'], data_df[f'C --> {name}'], alpha=0.3, lw=2)
            ax.plot(data_df[f'{name} benchmark t'], data_df[f'{name} benchmark C'], 'r--', alpha=0.3, lw=2)
            ax.set_title(name)
            #ax.set_ylim((500, 2500))
            ax.set_ylabel('Tracer conc.')
            ax.set_xlabel('Time')
    if makefigure:
        fig.tight_layout()
        fig.savefig('test_unsteady.pdf')
        fig.show()
# %%

def test_part_multiple(makefigure=True):
#%%
    makefigure=True
    timeseries_length = 300
    max_age = timeseries_length
    dt = 0.1
    Q_0 = 1.0 # <-- steady-state flow rate
    C_J = 1000 * np.ones(timeseries_length)
    C_eq = 00.
    C_old = 000.
    J = Q_0 
    f1, f2, f3 = 0.2, 0.2, 0.6
    S_0 = 10
    n_substeps = 1
    debug = False
    verbose = False
    jacobian = False

    data_df = pd.DataFrame(index=range(timeseries_length))
    data_df['t'] = data_df.index * dt
    data_df["J"] = J
    data_df["Q1"] = J * f1
    data_df["Q2"] = J * f2
    data_df["Q3"] = J * f3
    data_df["S_0"] = S_0
    data_df["C"] = C_J


    if makefigure:
        fig = plt.gcf()
    if True:
        name = 'Uniform'
        bm = steady_benchmarks[name]
        j = np.arange(timeseries_length)
        pQdisc = np.zeros_like(j, dtype=float)
        delta = dt * Q_0 / (S_0)
        pQdisc[0] = bm['pQdisc0'](delta) / dt
        pQdisc[1:] = bm['pQdisc'](delta, j[1:]) / dt
        data_df[f'{name} benchmark C'] = (np.convolve(C_J, pQdisc, mode='full')[:timeseries_length] * dt + C_old * (1-np.cumsum(pQdisc)[:timeseries_length]*dt))
        data_df[f'{name} benchmark t'] = data_df['t']

        solute_parameters = {
            "C":{
                "C_old": C_old,
                'alpha':{'Q1':0.5, 'Q2':1.5, 'Q3':1.}
            },
        }

        sas_specs = {
            'Q1':{
                'Q1_SAS':{
                    'ST': [0, S_0/3, S_0*2/3, S_0]
                }
            },
            'Q2':{
                'Q2_SAS':{
                    'ST': [0, S_0/3, S_0*2/3, S_0]
                }
            },
            'Q3':{
                'Q3_SAS':{
                    'ST': [0, S_0/2, S_0*3/4, S_0],
                    'P': [0, 1/2., 3/4., 1]
                }
            }
        }

        model = Model(
            data_df,
            sas_specs,
            solute_parameters,
            debug=debug,
            verbose=verbose,
            dt=dt,
            n_substeps=n_substeps,
            jacobian=jacobian,
            max_age=max_age,
            warning=True
        )
        model.run()
        data_df = model.data_df
        data_df[f'C --> {name}'] = f1*data_df[f'C --> Q1']+f2*data_df[f'C --> Q2']+f3*data_df[f'C --> Q3']
        err = (data_df[f'{name} benchmark C'].values-data_df[f'C --> {name}'].values)/data_df[f'{name} benchmark C'].values
        logging.info(f'Part/multiple error = {err.mean()}')
        assert err.mean()<1E-2
        if makefigure:
            ax = plt.subplot()
            ax.plot(data_df['t'], data_df[f'C --> {name}'], alpha=0.3, lw=2)
            ax.plot(data_df['t'], data_df[f'C --> Q1'], 'c', alpha=0.3, lw=2)
            ax.plot(data_df['t'], data_df[f'C --> Q2'], 'm', alpha=0.3, lw=2)
            ax.plot(data_df['t'], data_df[f'C --> Q3'], 'k', alpha=0.3, lw=2)
            ax.plot(data_df[f'{name} benchmark t'], data_df[f'{name} benchmark C'], 'k--', alpha=0.3, lw=2)
            ax.set_title(name)
            ax.legend()
            #ax.set_ylim((500, 2500))
            ax.set_ylabel('Tracer conc.')
            ax.set_xlabel('Time')
    if makefigure:
        fig.tight_layout()
        fig.savefig('test_part_multiple.pdf')
        fig.show()
# %%

def test_reaction(makefigure=True):
#%%
    u = 1 # No effect
    v = 1 # No effect
    # TRY LOOKING AT mR matrix?
    makefigure=True
    #timeseries_length = 5000
    #dt = 0.01
    timeseries_length = 500
    dt = 1 * u
    #timeseries_length = 50
    #dt = 0.1
    Q_0 = 1.0 / u * v # <-- steady-state flow rate
    C_J = 0. * np.ones(timeseries_length)
    C_eq = 1000.
    C_old = 1000.
    J = Q_0 
    S_0 = 100. * v
    k1 = 1/200. / u
    max_age = timeseries_length
    n_substeps = 1
    debug = False
    verbose = False
    jacobian = False

    data_df = pd.DataFrame(index=range(timeseries_length))
    data_df['t'] = data_df.index * dt
    data_df["J"] = J
    data_df["Q1"] = J
    data_df["S_0"] = S_0
    data_df["R"] = C_J


    if makefigure:
        fig = plt.gcf()
    if True:
        bm = steady_benchmarks['Uniform']
        j = np.arange(timeseries_length)
        pQdisc = np.zeros_like(j, dtype=float)
        delta = dt * Q_0 / (S_0)
        pQdisc[0] = bm['pQdisc0'](delta) / dt
        pQdisc[1:] = bm['pQdisc'](delta, j[1:]) / dt
        kappa = dt * k1
        pkdisc = np.zeros_like(j, dtype=float)
        pkdisc[0] = bm['pQdisc0'](kappa) / dt
        pkdisc[1:] = bm['pQdisc'](kappa, j[1:]) / dt
        data_df[f'Reaction benchmark R'] = (np.convolve((C_J/k1+C_eq/(Q_0/S_0)), pQdisc*pkdisc, mode='full')[:timeseries_length] * dt 
                                        + C_old * (1-np.cumsum(pkdisc)[:timeseries_length]*dt) * (1-np.cumsum(pQdisc)[:timeseries_length]*dt))
        data_df[f'Reaction benchmark t'] = data_df['t']

        solute_parameters = {
            "R":{
                "C_old": C_old,
                'k1': k1,
                'C_eq': C_eq
            }
        }

        sas_specs = {
            'Q1':{
                'Q1_SAS':{
                    'ST': [0, S_0/3, S_0*2/3, S_0]
                }
            }
        }

        model = Model(
            data_df,
            sas_specs,
            solute_parameters,
            debug=debug,
            verbose=verbose,
            dt=dt,
            n_substeps=n_substeps,
            jacobian=jacobian,
            max_age=max_age,
            warning=True
        )
        model.run()
        data_df = model.data_df
        err = (data_df[f'Reaction benchmark R'].values-data_df[f'R --> Q1'].values)/data_df[f'Reaction benchmark R'].values
        logging.info(f'Reaction error = {err.mean()}')
        assert err.mean()<1E-2
        if makefigure:
            ax = plt.subplot()
            ax.plot(data_df['t'], data_df[f'R --> Q1'], 'r', alpha=0.3, lw=2, label='mesas.py')
            ax.plot(data_df[f'Reaction benchmark t'], data_df[f'Reaction benchmark R'], 'r--', alpha=0.3, lw=2, label='benchmark')
            ax.set_title('First order reaction')
            ax.legend()
            ax.set_ylim((0, 1100))
            ax.set_ylabel('Tracer conc.')
            ax.set_xlabel('Time')
            ax.set_title(data_df[f'R --> Q1'].values[-1])
    if makefigure:
        fig.tight_layout()
        fig.savefig('test_reaction.pdf')
        fig.show()
# %%

if __name__=='__main__':
    test_steady()
    test_unsteady()
    test_part_multiple()
    test_reaction()