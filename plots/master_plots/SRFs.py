import numpy as np
import matplotlib.pyplot as plt
import h5py 
import time
import injections as inj
import functions as func
import pycbc
from scipy.stats import binned_statistic_2d
from pycbc import psd, detector
import lal
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
import os
lal.ClobberDebugLevel(0)
plt.rcParams.update({
    "text.usetex": True})
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

def read_FFs(filename):
	f = h5py.File(filename, 'r')

	FF = np.array(f['FF'])
	m1 = np.array(f['mass1'])
	m2 = np.array(f['mass2'])
	sigmasqs = np.array(f['sigmasqs'])
	return FF, sigmasqs, m1, m2

def get_injInside(binnumber, bin_x, bin_y):
    inj_with_row = np.where(binnumber[0] == bin_x+1)
    inj_with_col = np.where(binnumber[1] == bin_y+1)
    injInd = np.intersect1d(inj_with_row, inj_with_col)
    return injInd

def compute_SRF(param_x, param_y, nbinsx, nbinsy, FF, sigmasqs):
	start = time.time()
	ret = binned_statistic_2d(param_x, param_y, FF, statistic='count', bins=[nbinsx, nbinsy], expand_binnumbers=True)
	end = time.time()
	print('Binned statistics computed in', end-start)

	SRF = []
	indices = []
	ncells = int(nbinsx*nbinsy)
	binnumber = ret.binnumber
	count = ret.statistic
	xedge = ret.x_edge
	yedge = ret.y_edge

	for y in range(nbinsy):
		for x in range(nbinsx):
			numerator = 0.0
			denominator = 0.0
			injInd = np.array([])
			ninjs = int(count[x, y])

			if(ninjs !=0):	
				injInd = get_injInside(binnumber, x, y)
				for n in range(ninjs):
					inj_ind = injInd[n]
					numerator += FF[inj_ind]**3*sigmasqs[inj_ind]**3
					denominator += sigmasqs[inj_ind]**3
				temp_SRF = numerator/denominator
			else:
				temp_SRF = np.nan
				print('No injections found in bin', x, y)

			#print([x,y], 'Inj inside', len(injInd), 'SRF', temp_SRF)
			indices.append(injInd)
			SRF.append(temp_SRF)
	return np.array(SRF), indices, xedge, yedge

def plot_SRF(SRF, xedge, yedge, fig, ax, vmax, vmin):
	#levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
	cmap = plt.cm.YlGnBu
	#norm = BoundaryNorm(levels, ncolors = cmap.N, clip=True)	

	temp_x = xedge[:-1]
	temp_y = yedge[:-1]
	dx = xedge[1]-xedge[0]
	dy = yedge[1]-yedge[0]
	x_new = temp_x + dx/2.0
	y_new = temp_y + dy/2.0
	xx, yy = np.meshgrid(x_new, y_new)

	mtotal = xx + yy
	q = np.divide(xx, yy)
	z = SRF.reshape(len(y_new), len(x_new))
	#z[0,-1] = np.nan 
	#im = ax.pcolormesh(xx, yy, z, cmap=cmap, shading='auto', norm = MidpointNormalize(midpoint=1.0,vmin=vmin, vmax=vmax))
	im = ax.pcolormesh(xx, yy, z, cmap=cmap, vmax=vmax, vmin=vmin, shading='nearest') 
	return im


nbinsx = 8 
nbinsy = 7 
fig, axs = plt.subplots(2,2)

psd_list = ['ZDHP', 'aplus', 'voyager', 'ce2']
fig_title = ['Advanced LIGO design', 'A+', 'LIGO Voyager', 'Cosmic explorer']

for row in range(2):
	for col in range(2):
		#if (row*2 + col >= 3):
		#	break
		current_psd = psd_list[col + row*2]
		filename = '../%s/HM_prec/50000_FFs.hdf' %current_psd
		hf = h5py.File(filename, 'r')
		ax = axs[row, col]
	
		FF, sigmasqs, sg_m1, sg_m2 = read_FFs(filename)
		SRF, indices, xedge, yedge = compute_SRF(sg_m1, sg_m2, nbinsx, nbinsy, FF, sigmasqs)

		vmax = np.nanmax(SRF)
		vmin = np.nanmin(SRF)
		print('Vmax', vmax, 'Vmin', vmin)
		im = plot_SRF(SRF, xedge, yedge, fig, ax, 1.0, 0.6)
		ax.set_title(fig_title[col + row*2])

fig.subplots_adjust(right=0.8, hspace = 0.3)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax)
fig.text(0.5, 0.04, '$m_1^{det}$', ha='center', fontsize=13)
fig.text(0.04, 0.5, '$m_2^{det}$', va='center', rotation='vertical', fontsize=13)
fig.text(0.96, 0.4, 'Signal recovery fractions', ha='center', rotation=270, fontsize='large')

plt.savefig('SRF_diffPSDs.jpg', dpi=600)
plt.show()							
