from test_sas import steady_run
import matplotlib.pyplot as plt

import mesas.utils.vis as vis

timeseries_length =  50
max_age = timeseries_length
dt = 1/timeseries_length
Q_0 = 3.0/timeseries_length / dt  # <-- steady-state flow rate
C_J = 1000.
C_old = 0.
S_0 = 1
S_m = 0.1
eps = 0.0001/timeseries_length
n_substeps = 20
n_segment = 1
fQ=0.3
fc=0.1
jacobian = False

model = steady_run(timeseries_length, dt, Q_0, S_0, C_J, j=None, ST_min=S_m, debug=True, verbose=False, n_substeps=n_substeps, jacobian=jacobian, max_age=max_age, C_old=C_old)

flux = 'Q1'
sol = 'Cb'
i = 40
fig = plt.figure(figsize=[11.5,  4.])
axTC = plt.subplot2grid((2, 3), (0,1), rowspan=2)
axJ = plt.subplot2grid((2, 3), (0,0))
axQ = plt.subplot2grid((2, 3), (0,2))
axCJ = plt.subplot2grid((2, 3), (1,0))
axCQ = plt.subplot2grid((2, 3), (1,2))
vis.plot_transport_column(model, flux, sol, i=i, ax=axTC, nST=20, cmap='cividis_r', TC_frac=0.3, vrange=[0, C_J], ST_max=S_0)
vis.plot_influx(model, ax=axJ, sharex=axJ, i=i)
vis.plot_outflux(model, flux, ax=axQ, sharex=axJ, i=i)
vis.plot_influx_conc(model, sol, ax=axCJ, sharex=axJ, i=i)
vis.plot_outflux_conc(model, flux, sol, ax=axCQ, sharex=axJ, i=i)
plt.tight_layout()

#ani = vis.make_transport_column_animation(model, flux, sol, nST=20, cmap='cividis_r', TC_frac=0.3, vrange=[0, C_J], ST_max=S_0)
#ani.save(f'test_make_transport_column_animation.gif', writer='imagemagick')
