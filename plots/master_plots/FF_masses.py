import h5py 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import functions as func
import scipy
from scipy.interpolate import griddata
import scipy.ndimage as ndimage
import seaborn as sns
import pandas as pd
import os
import sys
import scipy.ndimage as ndimage
from scipy.stats import binned_statistic_2d
plt.rcParams.update({
    "text.usetex": True})
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

def read_FFs(filename):
	f = h5py.File(filename, 'r')
	FF = np.array(f['FF'][:])
	rec_tau = np.array(f['rec_tau'][:])
	sigmasqs = np.array(f['sigmasqs'][:])
	m1 = np.array(f['mass1'])
	m2 = np.array(f['mass2'])
	return FF, rec_tau, sigmasqs, m1, m2

def binned_plot(x, y, z, xbins, ybins, fig, ax, vmin, vmax):
    x_bins = np.linspace(min(x), max(x), xbins)
    y_bins = np.linspace(min(y), max(y), ybins)

    ret = binned_statistic_2d(x, y, z, statistic = np.mean, bins=[x_bins, y_bins])
    print('After smoothing max value', np.nanmax(ret.statistic), 'min value', np.nanmin(ret.statistic))
    im = ax.imshow(ret.statistic.T, origin='lower', extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]], vmin=vmin, vmax=vmax, aspect='auto', cmap = 'YlGnBu')
    return im

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

xbins = 10
ybins = 8
vmax = 1.0
vmin = 0.76
levels = np.linspace(vmin, vmax, 8)

psd_list = ['ZDHP', 'aplus', 'voyager', 'ce2']
fig_title = ['Advanced LIGO design', 'A+', 'LIGO Voyager', 'Cosmic explorer']
fig, axs = plt.subplots(2,2)
FF = np.array([])

for row in range(2):
	for col in range(2):
		#if (row*2 + col >= 3):
		#	break
		current_psd = psd_list[col + row*2]
		filename = "../%s/HM_prec/50000_FFs.hdf" %(current_psd)
		print(filename)
		FF, rec_tau, sigmasqs, m1, m2 = read_FFs(filename)
		mtotal = m1 + m2
		q = np.divide(m1, m2)
			
		ax = axs[row, col]
		im = binned_plot(m1, m2, FF, xbins, ybins, fig, ax, vmin, vmax)
		ax.set_title(fig_title[col+row*2])

fig.subplots_adjust(right=0.8, hspace = 0.25)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax)
fig.text(0.5, 0.04, '$m_1^{det}$', ha='center', fontsize=13)
fig.text(0.04, 0.5, '$m_2^{det}$', va='center', rotation='vertical', fontsize=13)
fig.text(0.96, 0.43, 'Fitting factors', ha='center', rotation=270, fontsize='large')

plt.savefig('FF_diffPSDs.png', dpi=600)
plt.show()
