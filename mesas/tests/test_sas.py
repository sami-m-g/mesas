# Todo: incorporate pytest-benchmark
import numpy as np
import pandas as pd
import pytest

# 2. A 'blender' that assumes that SAS function is fixed in time
from mesas.sas.blender import Fixed
# classes we will use to construct the SAS model
# 1. A piecewise constant SAS function
from mesas.sas.functions import Piecewise
# 3. The sas model class
from mesas.sas.model import Model

dt = 1.
Q_0 = 1.0 / dt  # <-- steady-state flow rate
C_J = 1000.
C_old = 0.
N = 5
S_0 = 5.
S_m = 0.601
eps = 0.0000001


def steady_run(N, dt, Q_0, S_0, C_J, eps1=(0, 0), ST_min=0., n_substeps=1):
    data_df = pd.DataFrame()
    data_df['Q1'] = np.ones(N) * Q_0
    data_df['J'] = np.ones(N) * Q_0
    data_df['Ca'] = np.ones(N) * C_J
    sas_fun1 = Piecewise(nsegment=1, ST_min=ST_min + eps1[0], ST_max=S_0 + eps1[1])
    sas_blends = {'Q1': Fixed(sas_fun1, N=len(data_df))}
    solute_parameters = {'Ca': {'C_old': C_old, 'observations': ['Q1']}}
    model = Model(data_df, sas_blends, solute_parameters, debug=True, verbose=False, dt=dt, n_substeps=n_substeps)
    print('running test_steady')
    model.run()
    return model

def test_steady_uniform():

    n = np.arange(N)
    T_0 = S_0 / Q_0
    Delta = dt / T_0
    Kappa = np.exp(-Delta)
    Eta = Kappa ** n

    sTdisc = -((Q_0 * Eta * (-1 + Kappa)) / Delta)
    pQdisc = (Q_0 * Eta * (-1 + Kappa) ** 2) / (S_0 * Delta ** 2 * Kappa)
    pQdisc[0] = (Q_0 * (-1 + Delta + n[0] * Delta + Eta[0] * Kappa)) / (S_0 * Delta ** 2)
    mQdisc = pQdisc * C_J * Q_0
    mTdisc = -((C_J * Q_0 * Eta * (-1 + Kappa)) / Delta)
    CQdisc = (C_J * (Delta + Eta * (-1 + Kappa)) - C_old *Eta * (-1 + Kappa)) / Delta

    CQdisc = np.array([[CQdisc]]).T
    sTdisc = np.tril(np.tile(sTdisc, [N, 1])).T
    sTdisc = np.c_[np.zeros(N), sTdisc]
    pQdisc = np.tril(np.tile(pQdisc, [N, 1]))
    pQdisc = np.array([pQdisc]).T
    mQdisc = np.tril(np.tile(mQdisc, [N, 1]))
    mQdisc = np.array([[mQdisc]]).T
    mTdisc = np.tril(np.tile(mTdisc, [N, 1])).T
    mTdisc = np.c_[np.zeros(N), mTdisc]
    mTdisc = np.array([mTdisc.T]).T

    model = steady_run(N, dt, Q_0, S_0, C_J)
    rdf = model.result

    def printcheck(rdf, varstr, analy):
        print(f'{varstr} Expected:')
        print(analy.T)
        print(f'{varstr} Got:')
        print(rdf[varstr].T)
        print(f'{varstr} Difference/expected:')
        err = (analy - rdf[varstr]) / analy
        print(err[..., -3:].T)
        assert np.nanmax(np.abs(err)) < 1.0E-4
        print('')

    printcheck(rdf, 'C_Q', CQdisc)
    printcheck(rdf, 'sT', sTdisc)
    printcheck(rdf, 'pQ', pQdisc)
    printcheck(rdf, 'mT', mTdisc)
    printcheck(rdf, 'mQ', mQdisc)

    print('Water Balance:')
    print(rdf['WaterBalance'][:, -3:] / Q_0)
    assert np.abs(rdf['WaterBalance'] / Q_0).max() < 1.0E-6

    print('Solute Balance:')
    for s in range(1):
        print(rdf['SoluteBalance'][:, -3:, s] / (Q_0 * C_J))
    assert np.abs(rdf['SoluteBalance'] / (Q_0 * C_J)).max() < 1.0E-6

    dsTdSjdisc = -((Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta))
    dmTdSjdisc = -((C_J * Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta))
    dCQdSjdisc = ((C_J - C_old) *Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 *Delta)

    dsTdSjdisc = np.tril(np.tile(dsTdSjdisc, [N, 1])).T
    dsTdSjdisc = np.c_[np.zeros(N), dsTdSjdisc]
    dmTdSjdisc = np.tril(np.tile(dmTdSjdisc, [N, 1])).T
    dmTdSjdisc = np.c_[np.zeros(N), dmTdSjdisc]
    dmTdSjdisc = np.array([dmTdSjdisc.T]).T
    dCQdSjdisc = np.array([[dCQdSjdisc]]).T

    model2 = steady_run(N, dt, Q_0, S_0, C_J, eps1=[0, eps])
    rdf2 = model2.result
    SAS_lookup, _, _, _, _, _, _ = model._create_sas_lookup()
    SAS_lookup2, _, _, _, _, _, _ = model2._create_sas_lookup()
    j = 1
    dSj = SAS_lookup2[j, N - 1] - SAS_lookup[j, N - 1]

    def printcheck(rdfi, rdfp, varstr, ostr, analy):
        if varstr=='dCdSj':
            var = rdfi[varstr][:,1,...]
        else:
            var = rdfi[varstr][:,:,1,...]
        print(f'{varstr} Expected:')
        print(analy.T)
        print(f'{varstr} eps check:')
        print(((rdfp[ostr] - rdfi[ostr]).T / dSj))
        print(f'{varstr} Got:')
        print(var.T)
        print(f'{varstr} Difference/expected:')
        err = (analy - var) / analy
        print(err[..., -3:].T)
        assert np.nanmax(np.abs(err)) < 1.0E-4
        print('')

    printcheck(rdf, rdf2, 'dsTdSj', 'sT', dsTdSjdisc)
    printcheck(rdf, rdf2, 'dmTdSj', 'mT', dmTdSjdisc)
    printcheck(rdf, rdf2, 'dCdSj', 'C_Q', dCQdSjdisc)

def test_steady_piston_uniform():

    n = np.arange(N)
    T_0 = S_0 / Q_0
    Delta = dt / T_0
    Kappa = np.exp(-Delta)
    Eta = Kappa ** n
    T_m = S_m / Q_0
    m = T_m / dt
    HeavisideTheta = lambda x: np.heaviside(x,1)


    sTdisc = Q_0 + (Q_0*Kappa**(m/(-1 + m*Delta))*(-((-1 + Delta + n*Delta)*Kappa**(m/(1 - m*Delta))) + (-1 + m*Delta)*Kappa**((1 + n)/(1 - m*Delta)) + ((-1 + n*Delta)*Kappa**(m/(1 - m*Delta)) + (1 - m*Delta)*Kappa**(n/(1 - m*Delta)))*HeavisideTheta(-m + n))*HeavisideTheta(1 - m + n))/Delta
    dsTdSjdisc = (Q_0*Kappa**((1 + n)/(1 - m*Delta))*((1 + (1 - 2*m + n)*Delta)*Kappa**(m/(-1 + m*Delta)) + (-1 + m*Delta)*Kappa**((1 + n)/(-1 + m*Delta)) + Kappa**(1/(-1 + m*Delta))*((-1 + 2*m*Delta - n*Delta)*Kappa**(m/(-1 + m*Delta)) + (1 - m*Delta)*Kappa**(n/(-1 + m*Delta)))*HeavisideTheta(-m + n))* HeavisideTheta(1 - m + n))/(S_0*Delta*(-1 + m*Delta))
    pQdisc = -((Q_0*((1 + Delta - n*Delta + (-1 + m*Delta)*Kappa**((1 + m - n)/(-1 + m*Delta)))*HeavisideTheta(((-1 - m + n)*S_0*Delta)/Q_0) + 2*(-1 + n*Delta + (1 - m*Delta)*Kappa**((m - n)/(-1 + m*Delta)))*HeavisideTheta(((-m + n)*S_0*Delta)/Q_0) + (1 - (1 + n)*Delta + (-1 + m*Delta)*Kappa**((1 - m + n)/(1 - m*Delta)))*HeavisideTheta(((1 - m + n)*S_0*Delta)/Q_0)))/(S_0*Delta**2))
    pQdisc[0] = (Q_0*(-1 + Delta + n[0]*Delta + (1 - m*Delta)*Kappa**((1 - m + n[0])/(1 - m*Delta)))*HeavisideTheta(((1 - m + n[0])*S_0*Delta)/Q_0))/ (S_0*Delta**2)
    mQdisc = pQdisc * Q_0 * C_J
    dpQdSjdisc = -((Q_0*Kappa**((1 + m - n)/(-1 + m*Delta))*((1 + (-1 - 2*m + n)*Delta + (-1 + m*Delta)*Kappa**((1 + m - n)/(1 - m*Delta)))* HeavisideTheta(((-1 - m + n)*S_0*Delta)/Q_0) + 2*((-1 + 2*m*Delta - n*Delta)*Kappa**(1/(1 - m*Delta)) + (1 - m*Delta)*Kappa**((1 + m - n)/(1 - m*Delta)))*HeavisideTheta(((-m + n)*S_0*Delta)/Q_0) - ((1 - m*Delta)*Kappa**((1 + m - n)/(1 - m*Delta)) + (-1 + (-1 + 2*m - n)*Delta)/Kappa**(2/(-1 + m*Delta)))*HeavisideTheta(((1 - m + n)*S_0*Delta)/Q_0)))/ (S_0**2*Delta**2*(-1 + m*Delta)))
    dpQdSjdisc[0] = (Q_0*(1 - m*Delta + (-1 + (-1 + 2*m - n[0])*Delta)*Kappa**((1 - m + n[0])/(1 - m*Delta)))* HeavisideTheta(((1 - m + n[0])*S_0*Delta)/Q_0))/(S_0**2*Delta**2*(-1 + m*Delta))
    mTdisc = (C_J*Q_0* (Delta + Kappa**(m/(-1 + m*Delta))*(-((-1 + Delta + n*Delta)*Kappa**(m/(1 - m*Delta))) + (-1 + m*Delta)*Kappa**((1 + n)/(1 - m*Delta)) + ((-1 + n*Delta)*Kappa**(m/(1 - m*Delta)) + (1 - m*Delta)*Kappa**(n/(1 - m*Delta)))*HeavisideTheta(-m + n))*HeavisideTheta(1 - m + n)))/Delta
    dmTdSjdisc = (C_J*Q_0*Kappa**((1 + n)/(1 - m*Delta))*((1 + (1 - 2*m + n)*Delta)*Kappa**(m/(-1 + m*Delta)) + (-1 + m*Delta)*Kappa**((1 + n)/(-1 + m*Delta)) + Kappa**(1/(-1 + m*Delta))*((-1 + 2*m*Delta - n*Delta)*Kappa**(m/(-1 + m*Delta)) + (1 - m*Delta)*Kappa**(n/(-1 + m*Delta)))*HeavisideTheta(-m + n))* HeavisideTheta(1 - m + n))/(S_0*Delta*(-1 + m*Delta))
    CQdisc = C_old - ((C_J - C_old)*Kappa**(m/(-1 + m*Delta))* (-((-1 + Delta + n*Delta)*Kappa**(m/(1 - m*Delta))) + (-1 + m*Delta)*Kappa**((1 + n)/(1 - m*Delta)) + ((-1 + n*Delta)*Kappa**(m/(1 - m*Delta)) + (1 - m*Delta)*Kappa**(n/(1 - m*Delta)))*HeavisideTheta(-m + n))*HeavisideTheta(1 - m + n))/Delta
    dCQdSjdisc = ((C_J - C_old)*(1 - m*Delta + (-1 + (-1 + 2*m - n)*Delta)*Kappa**((1 - m + n)/(1 - m*Delta)) + (-1 + m*Delta + (1 - 2*m*Delta + n*Delta)*Kappa**((m - n)/(-1 + m*Delta)))*HeavisideTheta(-m + n))*HeavisideTheta(1 - m + n))/(S_0*Delta*(-1 + m*Delta))

    dsTdSmdisc = -((Q_0*Kappa**((1 - m + n)/(1 - m*Delta))*(1 - m + n + (m - n)*Kappa**(1/(-1 + m*Delta))*HeavisideTheta(-m + n))* HeavisideTheta(1 - m + n))/(S_0*(-1 + m*Delta)))
    dpQdSmdisc = -((Q_0*Kappa**((1 - m + n)/(1 - m*Delta))*((1 + m - n)*Kappa**(2/(-1 + m*Delta))*HeavisideTheta(-1 - m + n) - 2*(m - n)*Kappa**(1/(-1 + m*Delta))*HeavisideTheta(-m + n) + (-1 + m - n)*HeavisideTheta(1 - m + n)))/(S_0**2*Delta*(-1 + m*Delta)))
    dpQdSmdisc[0] = -(((-1 + m - n[0])*Q_0*Kappa**((1 - m + n[0])/(1 - m*Delta))*HeavisideTheta(1 - m + n[0]))/(S_0**2*Delta*(-1 + m*Delta)))
    dmTdSmdisc = -((C_J*Q_0*Kappa**((1 - m + n)/(1 - m*Delta))*(1 - m + n + (m - n)*Kappa**(1/(-1 + m*Delta))*HeavisideTheta(-m + n))* HeavisideTheta(1 - m + n))/(S_0*(-1 + m*Delta)))
    dCQdSmdisc = ((C_J - C_old)*Kappa**((1 - m + n)/(1 - m*Delta))*(1 - m + n + (m - n)*Kappa**(1/(-1 + m*Delta))*HeavisideTheta(-m + n))* HeavisideTheta(1 - m + n))/(S_0*(-1 + m*Delta))

    CQdisc = np.array([[CQdisc]]).T
    sTdisc = np.tril(np.tile(sTdisc, [N, 1])).T
    sTdisc = np.c_[np.zeros(N), sTdisc]
    pQdisc = np.tril(np.tile(pQdisc, [N, 1]))
    pQdisc = np.array([pQdisc]).T
    mQdisc = np.tril(np.tile(mQdisc, [N, 1]))
    mQdisc = np.array([[mQdisc]]).T
    mTdisc = np.tril(np.tile(mTdisc, [N, 1])).T
    mTdisc = np.c_[np.zeros(N), mTdisc]
    mTdisc = np.array([mTdisc.T]).T

    dsTdSjdisc = np.tril(np.tile(dsTdSjdisc, [N, 1])).T
    dsTdSjdisc = np.c_[np.zeros(N), dsTdSjdisc]
    dmTdSjdisc = np.tril(np.tile(dmTdSjdisc, [N, 1])).T
    dmTdSjdisc = np.c_[np.zeros(N), dmTdSjdisc]
    dmTdSjdisc = np.array([dmTdSjdisc.T]).T
    dCQdSjdisc = np.array([[dCQdSjdisc]]).T

    dsTdSmdisc = np.tril(np.tile(dsTdSmdisc, [N, 1])).T
    dsTdSmdisc = np.c_[np.zeros(N), dsTdSmdisc]
    dmTdSmdisc = np.tril(np.tile(dmTdSmdisc, [N, 1])).T
    dmTdSmdisc = np.c_[np.zeros(N), dmTdSmdisc]
    dmTdSmdisc = np.array([dmTdSmdisc.T]).T
    dCQdSmdisc = np.array([[dCQdSmdisc]]).T


    model = steady_run(N, dt, Q_0, S_0, C_J, ST_min=S_m, n_substeps=200)
    rdf = model.result

    def printcheck(rdfi, varstr, analy):
        print(f'{varstr} Expected:')
        print(analy.T)
        print(f'{varstr} Got:')
        print(rdfi[varstr].T)
        print(f'{varstr} Difference/expected:')
        err = (analy - rdfi[varstr]) / analy
        print(err[..., -3:].T)
        assert np.nanmax(np.abs(err)) < 5.0E-2
        print('')

    #printcheck(rdf, 'pQ', pQdisc)
    #printcheck(rdf, 'sT', sTdisc)
    #printcheck(rdf, 'mQ', mQdisc)
    #printcheck(rdf, 'mT', mTdisc)
    #printcheck(rdf, 'C_Q', CQdisc)

    print('Water Balance:')
    print(rdf['WaterBalance'][:, -3:] / Q_0)
    assert np.abs(rdf['WaterBalance'] / Q_0).max() < 1.0E-6

    print('Solute Balance:')
    for s in range(1):
        print(rdf['SoluteBalance'][:, -3:, s] / (Q_0 * C_J))
    assert np.abs(rdf['SoluteBalance'] / (Q_0 * C_J)).max() < 1.0E-6

    def printcheck(rdfi, rdfp, varstr, ostr, analy, ip):
        if varstr=='dCdSj':
            var = rdfi[varstr][:,ip,...]
        else:
            var = rdfi[varstr][:,:,ip,...]
        print(f'{varstr} Expected:')
        print(analy.T)
        print(f'{varstr} eps check:')
        print(((rdfp[ostr] - rdfi[ostr]).T / dSj))
        print(f'{varstr} Got:')
        print(var.T)
        print(f'{varstr} Difference/expected:')
        err = (analy - var) / analy
        print(err[..., -3:].T)
        #assert np.nanmax(np.abs(err)) < 5.0E-2
        print('')

    model0 = steady_run(N, dt, Q_0, S_0, C_J, eps1=[0, eps], ST_min=S_m, n_substeps=200)
    rdf0 = model0.result
    SAS_lookup, _, _, _, _, _, _ = model._create_sas_lookup()
    SAS_lookup0, _, _, _, _, _, _ = model0._create_sas_lookup()
    j = 1
    dSj = SAS_lookup0[j, N - 1] - SAS_lookup[j, N - 1]

    #printcheck(rdf, rdf0, 'dsTdSj', 'sT', dsTdSjdisc, j)
    #printcheck(rdf, rdf0, 'dmTdSj', 'mT', dmTdSjdisc, j)
    #printcheck(rdf, rdf0, 'dCdSj', 'C_Q', dCQdSjdisc, j)

    modelm = steady_run(N, dt, Q_0, S_0, C_J, eps1=[eps, 0], ST_min=S_m, n_substeps=200)
    rdfm = modelm.result
    SAS_lookup, _, _, _, _, _, _ = model._create_sas_lookup()
    SAS_lookupm, _, _, _, _, _, _ = modelm._create_sas_lookup()
    j = 0
    dSj = SAS_lookupm[j, N - 1] - SAS_lookup[j, N - 1]

    printcheck(rdf, rdfm, 'dsTdSj', 'sT', dsTdSmdisc, j)
    printcheck(rdf, rdfm, 'dmTdSj', 'mT', dmTdSmdisc, j)
    printcheck(rdf, rdfm, 'dCdSj', 'C_Q', dCQdSmdisc, j)

def test_Jacobian():
    n = np.arange(N)
    T_0 = S_0 / Q_0
    Delta = dt / T_0
    Kappa = np.exp(-Delta)
    Eta = Kappa ** n

    dsTdSjdisc = -((Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta))
    dSTdSjdisc = (Q_0 - Q_0 * (1 + Delta + n * Delta) * Eta * Kappa) / (S_0 * Delta)
    dpQdSjdisc = (Q_0 * (-1 + (1 + Delta + n[0] * Delta) * Eta * Kappa)) / (S_0 ** 2 * Delta ** 2)
    dpQdSjdisc[0] = (Q_0 * Eta[0] * (-1 + Kappa) * (-1 + Kappa + Delta * (1 + n[0] * (-1 + Kappa) + Kappa))) / (
            S_0 ** 2 * Delta ** 2 * Kappa)
    dPQdSjdisc = (Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 ** 2 * Delta ** 2)
    dPQdSjdisc[0] = (Q_0 * (-1 + Kappa + Delta * Kappa)) / (S_0 ** 2 * Delta ** 2)
    dmQdSjdisc = dpQdSjdisc * C_J * Q_0
    dmTdSjdisc = -((C_J * Q_0 * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta))
    dCQdSjdisc = ((C_J + C_old) * Eta * (-1 + n * Delta * (-1 + Kappa) + Kappa + Delta * Kappa)) / (S_0 * Delta)

    dsTdSjdisc = np.tril(np.tile(dsTdSjdisc, [N, 1])).T
    dsTdSjdisc = np.c_[np.zeros(N), dsTdSjdisc]
    dCQdSjdisc = np.array([[dCQdSjdisc]]).T
    dpQdSjdisc = np.tril(np.tile(dpQdSjdisc, [N, 1]))
    dpQdSjdisc = np.array([dpQdSjdisc]).T
    dmQdSjdisc = np.tril(np.tile(dmQdSjdisc, [N, 1]))
    dmQdSjdisc = np.array([[dmQdSjdisc]]).T
    dmTdSjdisc = np.tril(np.tile(dmTdSjdisc, [N, 1])).T
    dmTdSjdisc = np.c_[np.zeros(N), dmTdSjdisc]
    dmTdSjdisc = np.array([dmTdSjdisc.T]).T

    ST_min = 1.5
    model = steady_run(N, dt, Q_0, S_0, C_J, ST_min=ST_min)
    rdf = model.result
    modeps = [[0, 0], [0, 0]]
    rdfeps = [[0, 0], [0, 0]]
    modeps[0][0] = steady_run(N, dt, Q_0, S_0, C_J, eps1=[eps, 0], ST_min=ST_min)
    modeps[0][1] = steady_run(N, dt, Q_0, S_0, C_J, eps1=[0, eps], ST_min=ST_min)
    # modeps[1][0] = steady_run(N, dt, Q_0, S_0, C_J, eps2=[eps, 0], ST_min=ST_min)
    # modeps[1][1] = steady_run(N, dt, Q_0, S_0, C_J, eps2=[0, eps], ST_min=ST_min)
    for i in range(1):
        for j in range(2):
            rdfeps[i][j] = modeps[i][j].result
    Jac = model.get_jacobian(mode='endpoint', logtransform=False)
    print('Jacobian: calculated')
    print(Jac[:, :2])
    print('Jacobian: should be')
    J = None
    for isol, sol in enumerate(model._solorder):
        if 'observations' in model.solute_parameters[sol]:
            for isolflux, solflux in enumerate(model._fluxorder):
                if solflux in model.solute_parameters[sol]['observations']:
                    J_seg = None
                    for iflux, flux in enumerate(model._comp2learn_fluxorder):
                        for ilabel, label in enumerate(model.sas_blends[flux]._comp2learn_componentorder):
                            i = iflux
                            J_seg_this = np.column_stack(
                                [((rdfeps[i][j]['C_Q'][:, isolflux, isol] - rdf['C_Q'][:, isolflux, isol]) / eps) for j
                                 in range(2)])
                            if J_seg is None:
                                J_seg = J_seg_this
                            else:
                                J_seg = np.c_[J_seg, J_seg_this]
                    if J is None:
                        J = J_seg
                    else:
                        J = np.concatenate((J, J_seg), axis=0)
    print(J)
