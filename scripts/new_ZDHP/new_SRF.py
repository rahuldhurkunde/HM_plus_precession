import numpy as np
import matplotlib.pyplot as plt
import h5py 
import time
import injections as inj
import functions as func
import pycbc
from scipy.interpolate import griddata
from scipy.stats import binned_statistic_2d
from pycbc import psd, detector
import lal
import seaborn as sns
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import pandas as pd
import argparse as ap
from matplotlib import colors
from matplotlib.patches import Rectangle
import os
lal.ClobberDebugLevel(0)

class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def read_FFs(inj_type):
	filename = "%s/50000_FFs.hdf" %(inj_type)
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
	levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
	cmap = plt.cm.YlGnBu
	norm = BoundaryNorm(levels, ncolors = cmap.N, clip=True)	

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
	#im = ax.pcolormesh(xx, yy, z, cmap=cmap, shading='auto', norm = MidpointNormalize(midpoint=1.0,vmin=vmin, vmax=vmax))
	im = ax.pcolormesh(xx, yy, z, cmap=cmap, vmax=vmax, vmin=vmin, shading='auto') 
	return im

nbinsx = 10 
nbinsy = 8 
inj_types = ['aligned', 'onlyHM', 'onlyPrecession', 'HM_prec']

#Command line parser
parser = ap.ArgumentParser()
parser.add_argument('--all',
			type=int, default = 0,
			help = 'Specify 1 to compute SRF (default option will only plot the SRF)')
parser.add_argument('--inj_class',
			type = int,
			help = 'Injection class (integer) 0-aligned, 1-onlyHM, 2-onlyPrecession, 3-HM_prec')
	
args, remaining_args = parser.parse_known_args()

if args.all == 1:
	fig, axs = plt.subplots(2,2)

	print("Computing SRF for bins", nbinsx, 'x', nbinsy)
	for row in range(2): 
		for col in range(2):
			current_inj = inj_types[row*2 + col]
			print(current_inj)
			ax = axs[row, col]	

			start = time.time()
			FF, sigmasqs, sg_m1, sg_m2 = read_FFs(current_inj)
			SRF, indices, xedge, yedge = compute_SRF(sg_m1, sg_m2, nbinsx, nbinsy, FF, sigmasqs)
			end = time.time()
			print("Time taken for SRF computation", end-start)

			vmax = np.nanmax(SRF)
			vmin = np.nanmin(SRF)
			#vmin = 0.1
			print('Vmax', vmax, 'Vmin', vmin)

			plot_SRF(SRF, xedge, yedge, fig, ax, 1.0, 0.6)
			ax.set_title(current_inj)
	plt.show()
