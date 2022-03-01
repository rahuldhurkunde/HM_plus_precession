import pycbc
import numpy as np
import matplotlib.pyplot as plt
import functions as func
import h5py

def read_FFs(inj_type, nsets):
	FF = np.array([])
	tau0 = np.array([])
	recovered_tau0 = np.array([])
	for k in range(nsets):
		filename = "../%s/output/H1-FF_%s-0-10.hdf" %(inj_type, k)
		f = h5py.File(filename, 'r')

		FF = np.concatenate([FF, np.array(f['FF'])])
		tau0 = np.concatenate([tau0, np.array(f['tau0'])])
		recovered_tau0 = np.concatenate([recovered_tau0, np.array(f['rec_tau'])])
	return FF, tau0, recovered_tau0

def read_tau_file(inj_type, tau_tolerance):
	filename = '/work/rahul.dhurkunde/HM_and_precession/injections/tau_files/tau_crawl_50000_%s.txt' %inj_type
	edges = np.loadtxt(filename)[:,0]
	statistics = np.loadtxt(filename)[:,1]
	tau_func = func.fit_tau_envelope(edges, statistics, tau_tolerance)
	return edges, tau_func

def visual_check(tau0, recovered_tau, edges, tau_func, tau_tolerance, fig, axs, current_inj):
	tau_diff = tau0-recovered_tau
	ax.scatter(tau0, np.abs(tau_diff), label=current_inj)

	x = np.linspace(min(edges), max(edges), 100)
	y = tau_func(x) + tau_tolerance
	ax.set_title(current_inj)
	ax.plot(x,y, color='red')

nsets = 10000
tau_tolerance = 0.35
inj_types = ['aligned', 'onlyHM', 'onlyPrecession', 'HM_prec']
fig, axs = plt.subplots(2,2)

for row in range(2):
	for col in range(2):
		current_inj = inj_types[col + 2*row]
		ax = axs[row, col]
		FF, tau0, recovered_tau0 = read_FFs(current_inj, nsets)
			
		edges, tau_func = read_tau_file(current_inj, tau_tolerance)
		visual_check(tau0, recovered_tau0, edges, tau_func, tau_tolerance, fig, axs, current_inj)

fig.subplots_adjust(right=0.8, hspace = 0.3)
fig.text(0.5, 0.04, 'signal $\\tau_0$', ha='center')
fig.text(0.04, 0.5, '$|\\tau_r - \\tau_0|$', va='center', rotation='vertical')
fig.text(0.5, 0.95, 'Tau recovery for ZDHP', ha='center', fontsize='x-large')
#plt.savefig('tau_recovery_zdhp.png', dpi=600)
plt.show()
