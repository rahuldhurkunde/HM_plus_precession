import h5py as hf
import numpy as np
import matplotlib.pyplot as plt
import scipy
import injections as inj
import functions as func
from scipy.interpolate import griddata
import scipy.ndimage as ndimage
from scipy.stats import binned_statistic_2d
import lalsimulation as lalsim
import pycbc
from pycbc import waveform, psd, filter, conversions
import lal
import time
import h5py
lal.ClobberDebugLevel(0)
plt.rcParams.update({
    "text.usetex": True})

def read_FFs(filename):
	hf = h5py.File(filename, 'r')
	FF = np.array(f['FF'][:])
	rec_tau = np.array(f['rec_tau'][:])
	sigmasqs = np.array(f['sigmasqs'][:])
	m1 = np.concatenate([m1, np.array(f['mass1'])])
	m2 = np.concatenate([m2, np.array(f['mass2'])])
	return FF, rec_tau, sigmasqs, m1, m2

def binned_plot(x, y, z, xbins, ybins, fig, ax, vmin, vmax):
	x_bins = np.linspace(min(x), max(x), xbins)
	y_bins = np.linspace(min(y), max(y), ybins)

	ret = binned_statistic_2d(x, y, z, statistic = np.mean, bins=[x_bins, y_bins])
	print('After smoothing max value', np.nanmax(ret.statistic), 'min value', np.nanmin(ret.statistic))
	im = ax.imshow(ret.statistic.T, origin='lower', extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]], vmin=vmin, vmax=vmax, aspect='auto', cmap = 'YlGnBu')
	#im = ax.imshow(ret.statistic.T, origin='lower', vmin=vmin, vmax=vmax)
	return im


filename = '50000_FFs.hdf'
vmax = 1.0                  
vmin = 0.76
fref = 100.0
pref = 0.0
f_min = 15.0
levels = np.linspace(vmin, vmax, 8)
xbins = 14
ybins = 12

fig, ax = plt.subplots(1, 1)

FF, rec_tau, sigmasqs, m1, m2 = read_FFs(filename)
im = binned_plot(m1, m2, FF, xbins, ybins, fig, ax, vmin, vmax)

cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.show()
