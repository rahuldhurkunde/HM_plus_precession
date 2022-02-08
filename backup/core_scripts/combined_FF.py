import h5py as hf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import functions as func
import scipy
from scipy.interpolate import griddata
import scipy.ndimage as ndimage
import seaborn as sns
import pandas as pd

def nonuniform_imshow(x, y, z, npoints, fig, ax):
	# Create regular grid
	xi, yi = np.linspace(x.min(), x.max(), npoints), np.linspace(y.min(), y.max(), npoints)
	xi, yi = np.meshgrid(xi, yi)

  	# Interpolate missing data
	zi = griddata((x, y), z, (xi, yi), method='linear')
	Z2 = ndimage.gaussian_filter(zi, sigma=0.5, order=0)
	im = plt.contourf(xi, yi, Z2, 15)
	plt.colorbar(im)
	#ax.imshow(xi, yi, Z2) 
	return 

def plot_tau_function(signal_tau, tau_diff):		
	tau_tolerance = 0.1
	bin_edges, statistic = func.compute_tauThreshold_envelope(signal_tau, tau_diff, 50)
	tau_func = func.fit_tau_envelope(bin_edges, statistic, tau_tolerance)
	plt.plot(signal_tau, tau_diff, 'x')
	plt.show()
	
def plot_cdf(FF):
	X = np.sort(FF)
	Y = np.array(range(len(FF)))/float(len(FF))
	plt.yscale('log')
	plt.plot(X,Y)
	plt.grid()
	plt.xlabel('Fitting factors')
	plt.ylabel('Cummulative distribution')
	plt.savefig('aligned_bank_cdf.png', dpi = 600)
	plt.show()

fig1, ax1 = plt.subplots()
FF = np.array([])
tau_diff = np.array([])
m1 = np.array([])
m2 = np.array([])

nsplits = 200 
signal_tau = np.array([]) 
for k in range(nsplits):
	filename = "../%s/output/TF2-FF-0-10.hdf" %k
	f = hf.File(filename, 'r')
	#print(f.keys())
	FF = np.concatenate([FF, np.array(f['FF'])])
	tau_diff = np.concatenate([tau_diff, abs(np.array(f['tau0']) - np.array(f['rec_tau']))])
	m1 = np.concatenate([m1, np.array(f['mass1'])])
	m2 = np.concatenate([m2, np.array(f['mass2'])])
	signal_tau = np.concatenate([signal_tau, np.array(f['tau0'])])

mtotal = m1 + m2
q = np.divide(m1, m2)
#im = plt.scatter(m1, m2, c = tau_diff)
#im = plt.scatter(m1, m2, c = FF)

nonuniform_imshow(mtotal, q, FF, 100, fig1, ax1)
#plt.colorbar(im)
#plot_cdf(FF)
#plt.xlabel('$m_1$')
#plt.ylabel('$m_2$')
#plt.title('Fitting factors (aligned_spin, no HM)')

#plt.savefig('aligned_bankFF.png', dpi=600)
plt.show()
